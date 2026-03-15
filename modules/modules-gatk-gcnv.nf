#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL FILE INSTANTIATION
// =====================================================================================
params {
    outdir = './output' // Change to your desired output path
    bin_length = 1000 // Example parameter for bin length
    padding = 100 // Example parameter for padding
    scatter_count = 10 // Example parameter for scatter count
    is_wgs = false // Set to false indicating the analysis is exome sequencing
}

// Define output directory
outdir = file(params.outdir, type: 'dir')

// Process to generate ploidy priors
process GENERATE_PLOIDY_PRIORS {
    tag "Ploidy_Priors"
    publishDir "${outdir}/out_GCNV/priors", mode: 'copy', overwrite: true

    input: 
    path fai
    
    output: 
    path "ploidy_priors.tsv", emit: priors
    
    script:
    """
    # Create header with 4 ploidy states (0, 1, 2, 3)
    echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" > ploidy_priors.tsv

    # Default for all contigs: high priority for Ploidy 2 (Autosomes)
    awk -v OFS="\t" '{print \$1, "0.01", "0.01", "0.97", "0.01"}' $fai >> ploidy_priors.tsv

    # Adjust Sex Chromosomes (X and Y) for mixed cohort
    sed -i 's/^chrX\t0.01\t0.01\t0.97\t0.01/chrX\t0.01\t0.49\t0.49\t0.01/' ploidy_priors.tsv
    sed -i 's/^chrY\t0.01\t0.01\t0.97\t0.01/chrY\t0.50\t0.50\t0.00\t0.00/' ploidy_priors.tsv
    """
}

// Process to preprocess intervals
process PREPROCESS_INTERVALS {
    tag "Preprocess_Intervals"
    label 'gatk'
    
    input: 
    path fasta
    path fai
    path dict
    path targets
    
    output: 
    path "preprocessed.interval_list", emit: interval_list
    
    script:
    """
    gatk --java-options "-Xmx4g" PreprocessIntervals \\
        -R $fasta \\
        -L $targets \\
        --bin-length ${params.bin_length} \\
        --padding ${params.padding} \\
        -imr OVERLAPPING_ONLY \\
        -O preprocessed.interval_list
    """
}

// Process to annotate intervals
process ANNOTATE_INTERVALS {
    tag "Annotate_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    path fasta
    path fai
    path dict
    
    output: 
    path "annotated.tsv", emit: annotated_intervals
    
    script:
    """
    gatk --java-options "-Xmx4g" AnnotateIntervals \\
        -L $intervals -R $fasta -imr OVERLAPPING_ONLY -O annotated.tsv
    """
}

// Process to collect read counts
process COLLECT_READ_COUNTS {
    tag "$sample_id"
    label 'gatk'
    
    input: 
    tuple val(sample_id), path(reads), path(index)
    path intervals
    path fasta
    path fai
    path dict

    output: 
    tuple val(sample_id), path("${sample_id}.counts.hdf5"), emit: counts

    script:
    """
    gatk --java-options "-Xmx4g" CollectReadCounts \\
        -L $intervals -R $fasta -I $reads \\
        --format HDF5 -imr OVERLAPPING_ONLY \\
        -O ${sample_id}.counts.hdf5
    """
}

// Process to filter intervals
process FILTER_INTERVALS {
    tag "Filter_Intervals"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path annotated
    path counts
    
    output: 
    path "cohort.filtered.interval_list", emit: filtered_intervals
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" FilterIntervals \\
        -L $intervals \\
        --annotated-intervals $annotated \\
        -I counts_files.list \\
        -imr OVERLAPPING_ONLY \\
        -O cohort.filtered.interval_list
    """
}

// Process to determine ploidy for the cohort
process DETERMINE_PLOIDY_COHORT {
    tag "Determine_Ploidy"
    label 'gatk'
    label 'large_mem'
    
    input: 
    path intervals
    path counts
    path priors
    
    output: 
    path "ploidy-calls", emit: ploidy_calls
    
    script:
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" DetermineGermlineContigPloidy \\
        -L $intervals \\
        -I counts_files.list \\
        --contig-ploidy-priors $priors \\
        --output . \\
        --output-prefix ploidy \\
        -imr OVERLAPPING_ONLY
    """
}

// Process to scatter intervals
process SCATTER_INTERVALS {
    tag "Scatter_Intervals"
    label 'gatk'
    
    input: 
    path intervals
    
    output: 
    path "temp_*/scattered.interval_list", emit: shards
    
    script:
    """
    gatk --java-options "-Xmx4g" IntervalListTools \\
        --INPUT $intervals \\
        --SUBDIVISION_MODE INTERVAL_COUNT \\
        --SCATTER_CONTENT ${params.scatter_count} \\
        --OUTPUT .
    """
}

// Process germline CNV calling for the cohort
process GERMLINE_CNV_CALLER_COHORT {
    tag "${interval_shard.parent.name}"
    label 'gatk'
    label 'gpu_or_high_cpu'
    
    input: 
    path interval_shard
    path annotated
    path counts
    path ploidy_calls
    
    output: 
    path "shard-calls", emit: calls
    path "shard-model", emit: model
    
    script:
    def sensitivity_params = "--class-coherence-length 1000.0 --cnv-coherence-length 1000.0 --enable-bias-factors false --interval-psi-scale 1.0E-6 --log-mean-bias-standard-deviation 0.01 --sample-psi-scale 1.0E-6"
    """
    echo "${counts.join('\n')}" > counts_files.list
    gatk --java-options "-Xmx16g" GermlineCNVCaller \\
        --run-mode COHORT \\
        -L $interval_shard \\
        -I counts_files.list \\
        --contig-ploidy-calls $ploidy_calls \\
        --annotated-intervals $annotated \\
        --output . \\
        --output-prefix shard \\
        -imr OVERLAPPING_ONLY \\
        $sensitivity_params
    """
}

// Process to post-process CNV calls for individual samples
process POSTPROCESS_CALLS {
    tag "$sample_id"
    label 'gatk'
    publishDir "${outdir}/out_GCNV/vcfs", mode: 'copy', overwrite: true
    
    input: 
    tuple val(index), val(sample_id)
    path model_shards
    path call_shards
    path ploidy_calls
    path dict
    
    output: 
    path "${sample_id}.genotyped_segments.vcf.gz", emit: final_vcf
    path "${sample_id}.genotyped_intervals.vcf.gz", emit: intervals_vcf
    
    script:
    def model_args = model_shards.collect { "--model-shard-path \$it" }.join(" ")
    def call_args = call_shards.collect { "--calls-shard-path \$it" }.join(" ")
    """
    gatk --java-options "-Xmx8g" PostprocessGermlineCNVCalls \\
        \$call_args \\
        \$model_args \\
        --contig-ploidy-calls $ploidy_calls \\
        --sample-index $index \\
        --allosomal-contig chrX --allosomal-contig chrY \\
        --sequence-dictionary $dict \\
        --output-genotyped-intervals ${sample_id}.genotyped_intervals.vcf.gz \\
        --output-genotyped-segments ${sample_id}.genotyped_segments.vcf.gz
    """
}
