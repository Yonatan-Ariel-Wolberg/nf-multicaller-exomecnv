#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// =====================================================================================
// GLOBAL PARAMETERS
// =====================================================================================
params {
    outdir = './output' // Change to your desired output path
    test_size = -1 // Specify the test size for performance tuning
}

// Define output directory
outdir = file(params.outdir, type: 'dir')

// =====================================================================================
// HELPER FUNCTION
// =====================================================================================
def is_large_run() {
    // Safely fetch test_size, defaulting to -1 if missing
    def size = params.get('test_size', -1).toInteger()
    return (size > 50 || size == -1)
}

// =====================================================================================
// PROCESSES FOR CNVKIT
// =====================================================================================

process GENERATE_ACCESS {
    tag "Access"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    
    input: 
    path fasta
    
    output: 
    path "access.hg38.bed", emit: access_bed
    
    script: 
    """
    cnvkit.py access $fasta -o access.hg38.bed
    """
}

process AUTOBIN {
    tag "Autobin"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    
    input: 
    path fasta
    path targets
    path access
    path refflat
    path bam 
    
    output: 
    path "*.target.bed", emit: target_bed
    path "*.antitarget.bed", emit: antitarget_bed
    
    script: 
    """
    cnvkit.py autobin -t $targets -g $access --fasta $fasta --annotate $refflat --short-names $bam
    """
}

process COVERAGE {
    tag "$sample_id"
    label 'cnvkit'
    memory { is_large_run() ? '16 GB' : '8 GB' }
    
    input: 
    tuple val(sample_id), path(bam), path(bai)
    path target_bed
    path antitarget_bed
    
    output: 
    tuple val(sample_id), path("${sample_id}.targetcoverage.cnn"), emit: target_cov
    tuple val(sample_id), path("${sample_id}.antitargetcoverage.cnn"), emit: antitarget_cov
    
    script:
    """
    cnvkit.py coverage $bam $target_bed -o ${sample_id}.targetcoverage.cnn
    cnvkit.py coverage $bam $antitarget_bed -o ${sample_id}.antitargetcoverage.cnn
    """
}

process CREATE_POOLED_REFERENCE {
    tag "Pooled_Ref"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/reference", mode: 'copy', overwrite: true
    memory { is_large_run() ? '64 GB' : '24 GB' }
    
    input: 
    path fasta
    path covs 
    
    output: 
    path "pooled_reference.cnn", emit: ref_cnn
    
    script: 
    """
    cnvkit.py reference ${covs} --fasta $fasta -o pooled_reference.cnn
    """
}

process CALL_CNV {
    tag "$sample_id"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/calls/${sample_id}", mode: 'copy', overwrite: true
    
    input: 
    tuple val(sample_id), path(t_cov), path(t_anticov)
    path reference
    
    output: 
    tuple val(sample_id), path("${sample_id}.cnr"), path("${sample_id}.cns"), emit: results
    path "*.pdf"
    
    script:
    """
    cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
    cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
    cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-scatter.pdf
    cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-diagram.pdf
    """
}

process EXPORT_RESULTS {
    tag "$sample_id"
    label 'cnvkit'
    publishDir "${outdir}/out_CNVKIT/vcfs", mode: 'copy', overwrite: true
    
    input: 
    tuple val(sample_id), path(cnr), path(cns)
    
    output: 
    path("${sample_id}_CNVKIT_output.vcf"), emit: vcf
    path("${sample_id}_calls.bed"), emit: bed
    
    script:
    """
    # Using specific naming convention so it routes cleanly into the consensus modules
    cnvkit.py export vcf $cns -i $sample_id -o ${sample_id}_CNVKIT_output.vcf
    cnvkit.py export bed $cns -i $sample_id -o ${sample_id}_calls.bed
    """
}
