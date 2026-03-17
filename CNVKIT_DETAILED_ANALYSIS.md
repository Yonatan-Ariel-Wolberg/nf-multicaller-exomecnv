# CNVKit Pipeline Analysis: Cohort-Only (Panel of Normals) Exome Sequencing

## 1. CORRECT CNVKIT WORKFLOW FOR EXOME DATA WITH COHORT ONLY

### Recommended Approach
For exome sequencing with **only tumor/proband samples (no normals)**, CNVKit provides two main strategies:

#### Option A: Build "Flat Reference" (Recommended for cohort-only)
When you have no normal samples available, create a "flat reference" - a uniform reference assuming equal coverage in all bins:

```bash
cnvkit.py reference -o FlatReference.cnn -f ucsc.hg38.fa -t targets.bed -a antitargets.bed
```

Key characteristics:
- Creates neutral copy number (log2 0.0) for each probe
- Still computes GC content and RepeatMasker information from FASTA
- Less accurate than pooled reference but useful without normals
- Can use large-scale CNVs should still be visible

#### Option B: Build Pooled Reference from Cohort Samples
Use all your cohort samples to build a pooled reference, then use this reference to analyze all samples:

```bash
# Step 1: Run coverage on all samples
for sample in *.bam; do
  cnvkit.py coverage $sample targets.target.bed -o ${sample%.bam}.targetcoverage.cnn
  cnvkit.py coverage $sample antitargets.antitarget.bed -o ${sample%.bam}.antitargetcoverage.cnn
done

# Step 2: Build pooled reference from ALL samples' coverage
cnvkit.py reference *targetcoverage.cnn *antitargetcoverage.cnn --fasta hg38.fa -o pooled_reference.cnn

# Step 3: Run fix, segment, call on each sample with the pooled reference
for sample in *.bam; do
  sample_id=${sample%.bam}
  cnvkit.py fix ${sample_id}.targetcoverage.cnn ${sample_id}.antitargetcoverage.cnn pooled_reference.cnn \
    -o ${sample_id}.cnr
  cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
  cnvkit.py call ${sample_id}.cns -o ${sample_id}.call.cns
done
```

**Why this works**: CNVKit can build a pooled reference from tumor/proband samples themselves. While not ideal, benchmarking shows a pooled reference performs better than a flat reference. The documentation explicitly states: "you can create a pooled reference even if matched tumor-normal pairs were sequenced -- our benchmarking showed that a pooled reference performed slightly better."

---

## 2. `cnvkit.py autobin` COMMAND

### Purpose
Quickly estimate read counts/depths in BAM files to estimate **reasonable on- and off-target bin sizes**.

### Output
- Target and antitarget BED files
- Table of estimated average read depths and recommended bin sizes

### Inputs and Requirements

**KEY POINT**: **If multiple BAMs are given, uses the BAM with median file size - NOT all BAMs**

```bash
# Correct usage
cnvkit.py autobin *.bam -t my_targets.bed -g access.hg38.bed
cnvkit.py autobin *.bam -m amplicon -t my_targets.bed
cnvkit.py autobin *.bam -m wgs -b 50000 -g access.hg38.bed --annotate refFlat.txt
```

### How It Works
1. Uses BAM index (.bai) to quickly determine total reads present
2. Random sampling of targeted regions to estimate average on-target read depth
3. Much faster than full `coverage` command
4. When multiple BAMs given: **selects the one with median file size** for bin size estimation

### Required Files
- BAM files (with .bai index) 
- For hybrid capture: `-t targets.bed` and `-g access.bed`
- For WGS: `-m wgs`, `-b bin_size`, `-g access.bed`
- For amplicon: `-m amplicon`, `-t targets.bed`

---

## 3. `cnvkit.py reference` COMMAND

### Purpose
Compile a copy-number reference from normal samples or create a flat reference.

### Inputs - CRITICAL

**The command takes BOTH target AND antitarget coverage files (.cnn):**

```bash
# Correct: all *.cnn files (both target and antitarget)
cnvkit.py reference *coverage.cnn -f ucsc.hg38.fa -o Reference.cnn

# Or explicitly:
cnvkit.py reference *targetcoverage.cnn *antitargetcoverage.cnn -f ucsc.hg38.fa -o Reference.cnn

# NOT just one type:
# WRONG: cnvkit.py reference *.targetcoverage.cnn
```

### Required Files
- **Both** `*.targetcoverage.cnn` and `*.antitargetcoverage.cnn` files from `coverage` command
- File name prefixes must match for pairs (e.g., `Sample1.targetcoverage.cnn` and `Sample1.antitargetcoverage.cnn`)
- Reference genome FASTA (optional but recommended for GC/RepeatMasker)

### What It Does
1. Median-centers each input sample
2. Performs bias corrections (GC, RepeatMasker) on each sample
3. Calculates weighted average (Tukey's biweight location) at each bin across all samples
4. Calculates spread (Tukey's biweight midvariance) - reliability estimate
5. Adds GC and repeat-masked proportion if FASTA provided
6. Outputs reference.cnn with:
   - "log2" column: expected read depth
   - "spread" column: reliability of estimate
   - "gc" column: GC fraction (if FASTA given)
   - "rmask" column: repeat-masked fraction (if FASTA given)

### For Cohort-Only
Include all samples' coverage files:
```bash
cnvkit.py reference sample1.targetcoverage.cnn sample1.antitargetcoverage.cnn \
                     sample2.targetcoverage.cnn sample2.antitargetcoverage.cnn \
                     ... -f hg38.fa -o pooled_reference.cnn
```

---

## 4. `cnvkit.py segment` vs `cnvkit.py call`

### `cnvkit.py segment` Command
- **INPUT**: .cnr file (copy number ratios from `fix` command)
- **OUTPUT**: .cns file (segmented copy numbers with breakpoints)
- **FUNCTION**: Infers discrete copy number segments using various segmentation algorithms
- **Algorithms**: CBS (default - Circular Binary Segmentation), Haar, HMM, etc.
- Runs independently on each chromosome arm
- Can be parallelized with `-p` option

Example:
```bash
cnvkit.py segment Sample.cnr -o Sample.cns
```

### `cnvkit.py call` Command
- **INPUT**: .cns file (segmented log2 ratios)
- **OUTPUT**: Same .cns file with additional "cn" column (absolute integer copy numbers)
- **FUNCTION**: Converts log2 ratios → absolute integer copy numbers
- **Methods**:
  - `threshold`: Apply fixed log2 ratio cutoffs for each copy number state
  - `clonal`: Simple rounding to nearest integer copy number
  - `none`: No copy number calling, just rescaling/re-centering

Examples:
```bash
cnvkit.py call Sample.cns -o Sample.call.cns
cnvkit.py call Sample.cns -m threshold -t=-1.1,-0.4,0.3,0.7 -o Sample.call.cns
cnvkit.py call Sample.cns -m clonal --purity 0.65 -o Sample.call.cns
```

### Typical Workflow After Segmentation
**YES, `call` is needed** before exporting to VCF:

```bash
# Step 1: Segment
cnvkit.py segment Sample.cnr -o Sample.cns

# Step 2: Call integer copy numbers
cnvkit.py call Sample.cns -o Sample.call.cns

# Step 3: Export to VCF
cnvkit.py export vcf Sample.call.cns -o Sample.cnv.vcf
```

The `call` command adds the "cn" (copy number) column that export vcf needs.

---

## 5. `cnvkit.py export vcf` COMMAND

### Purpose
Convert .cns files to VCF format for use in other tools.

### Key Parameters

#### `-i` / `--sample-id` FLAG
- **Does**: Sets the SAMPLE NAME in the VCF file
- **Required for per-sample VCFs**
- **Example**:
```bash
cnvkit.py export vcf Sample.cns -i "SampleID" -o Sample.cnv.vcf
```

#### Correct Full Command for Per-Sample VCFs
```bash
# Required: Use after 'call' command has added 'cn' column
cnvkit.py export vcf Sample.call.cns -i "SampleID" -o Sample.cnv.vcf

# With optional parameters:
cnvkit.py export vcf Sample.call.cns -i "SampleID" -x female -y -o Sample.cnv.vcf
```

### Other Important Parameters
- `-x` / `--sample-sex`: Specify sample sex (male/female) - affects X/Y chromosome handling
- `-y` / `--male-reference` / `--haploid-x-reference`: Use if reference was built with `-y`
- `--show all`: Include all segments (by default only shows regions with CN != ploidy)

### What Gets Exported
- Input must be .cns file with "cn" column (from `call` command)
- Creates VCF records only for segments where CN differs from expected ploidy
- For autosomes: CN != 2
- For sex chromosomes: depends on sample sex and reference sex

### Filter/Merge Before Export
The `call` command supports filtering and merging:
```bash
cnvkit.py call Sample.cns --filter ci --merge bic -o Sample.call.cns
cnvkit.py export vcf Sample.call.cns -i "SampleID" -o Sample.cnv.vcf
```

---

## 6. RECOMMENDED APPROACH FOR PER-SAMPLE VCFs

### Complete Workflow for Cohort Exome Sequencing

```bash
# ============ REFERENCE BUILDING (once for all samples) ============
# Step 1: Generate accessible regions
cnvkit.py access hg38.fa -o access.hg38.bed

# Step 2: Estimate bin sizes from median-sized sample
cnvkit.py autobin sample1.bam sample2.bam sample3.bam ... \
  -t baits.bed -g access.hg38.bed --annotate refFlat.txt \
  -o autobin_report.txt

# Step 3: Target/antitarget preparation
cnvkit.py target baits.bed --annotate refFlat.txt -o targets.bed
cnvkit.py antitarget targets.bed -o antitargets.bed

# ============ COVERAGE CALCULATION (per sample, parallelizable) ============
# Step 4: For each sample, calculate coverage
for sample in *.bam; do
  sample_id=${sample%.bam}
  cnvkit.py coverage $sample targets.bed -o ${sample_id}.targetcoverage.cnn
  cnvkit.py coverage $sample antitargets.bed -o ${sample_id}.antitargetcoverage.cnn
done

# ============ POOLED REFERENCE BUILDING (once) ============
# Step 5: Build pooled reference from ALL samples
cnvkit.py reference *targetcoverage.cnn *antitargetcoverage.cnn \
  --fasta hg38.fa -o pooled_reference.cnn

# ============ PER-SAMPLE CALLING (per sample, parallelizable) ============
# Step 6: For each sample, call CNVs
for sample in *.bam; do
  sample_id=${sample%.bam}
  
  # Normalize to reference
  cnvkit.py fix ${sample_id}.targetcoverage.cnn \
            ${sample_id}.antitargetcoverage.cnn \
            pooled_reference.cnn -o ${sample_id}.cnr
  
  # Segment
  cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
  
  # Call integer copy numbers
  cnvkit.py call ${sample_id}.cns -o ${sample_id}.call.cns
  
  # Export to VCF
  cnvkit.py export vcf ${sample_id}.call.cns \
    -i "${sample_id}" -o ${sample_id}.cnv.vcf
done
```

### Key Points for Cohort Approach
1. **Build reference once**: Include ALL samples' coverage files
2. **Use pooled reference for all**: Better than flat reference
3. **Parallelizable**: coverage and per-sample calling can run in parallel
4. **Consistent filtering**: Apply same filters to all samples for comparability

---

# ANALYSIS OF LOCAL NEXTFLOW MODULE: modules-cnvkit.nf

## ISSUE 1: AUTOBIN PROCESS (Line 45-65)

### Current Implementation
```nextflow
process AUTOBIN {
    input: 
    path fasta
    path targets
    path access
    path refflat
    path bam              // <-- PROBLEM: Takes variable number of BAMs
    
    script: 
    """
    cnvkit.py autobin -t $targets -g $access --fasta $fasta --annotate $refflat --short-names $bam
    """
}
```

### ISSUE IDENTIFIED
❌ **INCORRECT**: The process passes `$bam` (which expands to all BAM files) to autobin

✅ **CORRECT**: According to CNVKit documentation, `autobin` with multiple BAMs uses the **median-sized BAM only**

### FIX REQUIRED
```nextflow
process AUTOBIN {
    input: 
    path fasta
    path targets
    path access
    path refflat
    path bams              // Multiple BAMs from workflow
    
    script: 
    """
    # Select median-sized BAM file for bin size estimation
    MEDIAN_BAM=\$(python3 << 'PYTHON'
import os
bams = [f for f in "${bams.join(' ')}".split() if f.endswith('.bam')]
if len(bams) > 1:
    sizes = [(b, os.path.getsize(b)) for b in bams]
    sizes.sort(key=lambda x: x[1])
    print(sizes[len(sizes)//2][0])
else:
    print(bams[0])
PYTHON
)
    
    cnvkit.py autobin \$MEDIAN_BAM -t $targets -g $access --fasta $fasta --annotate $refflat --short-names
    """
}
```

**Alternatively, use CNVKit batch command directly** which handles this automatically:
```nextflow
process AUTOBIN_VIA_BATCH {
    input: 
    path fasta
    path targets
    path access
    path refflat
    path bams
    
    script: 
    """
    # Let batch command handle bin sizing from all BAMs
    cnvkit.py batch ${bams.join(' ')} -t $targets -g $access --fasta $fasta \
      --annotate $refflat --output-dir autobin_out/ -n
    """
}
```

---

## ISSUE 2: CREATE_POOLED_REFERENCE PROCESS (Line 88-105)

### Current Implementation
```nextflow
process CREATE_POOLED_REFERENCE {
    input: 
    path fasta
    path covs              // Mixed target and antitarget .cnn files
    
    script: 
    """
    cnvkit.py reference *.cnn --fasta $fasta -o pooled_reference.cnn
    """
}
```

### ISSUE IDENTIFIED
⚠️ **POTENTIALLY OKAY BUT UNCLEAR**: The command `*.cnn` will match ALL .cnn files

✅ **BEST PRACTICE**: Explicitly separate target and antitarget files

### WHY THIS WORKS (but could be better)
- CNVKit's `reference` command accepts both targetcoverage.cnn and antitargetcoverage.cnn
- File name prefixes automatically determine pairing (e.g., Sample1.targetcoverage.cnn pairs with Sample1.antitargetcoverage.cnn)
- Using `*.cnn` will include all files and CNVKit will handle pairing correctly

### RECOMMENDED FIX (Better Clarity)
```nextflow
process CREATE_POOLED_REFERENCE {
    input: 
    path fasta
    path target_covs       // Channel of *.targetcoverage.cnn
    path antitarget_covs   // Channel of *.antitargetcoverage.cnn
    
    script: 
    """
    cnvkit.py reference ${target_covs.join(' ')} ${antitarget_covs.join(' ')} \
      --fasta $fasta -o pooled_reference.cnn
    """
}
```

**Or, if files arrive mixed** (current approach):
```nextflow
process CREATE_POOLED_REFERENCE {
    input: 
    path fasta
    path covs              // Can be mixed target + antitarget
    
    script: 
    """
    # Current approach works - CNVKit will automatically pair by sample ID
    cnvkit.py reference *.cnn --fasta $fasta -o pooled_reference.cnn
    """
}
```

The pairing is automatic because file names like:
- Sample1.targetcoverage.cnn
- Sample1.antitargetcoverage.cnn
- Sample2.targetcoverage.cnn
- Sample2.antitargetcoverage.cnn

CNVKit matches them by prefix.

---

## ISSUE 3: CALL_CNV PROCESS - Workflow (Line 107-127)

### Current Implementation
```nextflow
process CALL_CNV {
    input: 
    tuple val(sample_id), path(t_cov), path(t_anticov)
    path reference
    
    script:
    """
    cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
    cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
    cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-scatter.pdf
    cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.cns -o ${sample_id}-diagram.pdf
    """
}
```

### ISSUE IDENTIFIED
❌ **MISSING**: No call to `cnvkit.py call` after segmentation!

### REQUIRED FIX
```nextflow
process CALL_CNV {
    input: 
    tuple val(sample_id), path(t_cov), path(t_anticov)
    path reference
    
    script:
    """
    # Step 1: Normalize to reference
    cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
    
    # Step 2: Segment
    cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
    
    # Step 3: CRITICAL - Call integer copy numbers for VCF export
    cnvkit.py call ${sample_id}.cns -o ${sample_id}.call.cns
    
    # Step 4: Visualizations
    cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.call.cns -o ${sample_id}-scatter.pdf
    cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.call.cns -o ${sample_id}-diagram.pdf
    """
}
```

### Why This Matters
- The `segment` command outputs log2 ratios in .cns files
- The `call` command converts these log2 ratios → integer copy numbers (adds "cn" column)
- The `export vcf` command REQUIRES the "cn" column from `call`
- Without `call`, VCF export will fail or produce incorrect output

**From CNVKit documentation workflow:**
```
cnvkit.py fix Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn my_reference.cnn -o Sample.cnr
cnvkit.py segment Sample.cnr -o Sample.cns
cnvkit.py call Sample.cns -o Sample.call.cns  ← REQUIRED before export
```

### Updated Output Declaration
```nextflow
output: 
    tuple val(sample_id), path("${sample_id}.cnr"), path("${sample_id}.call.cns"), emit: results
    path "*.pdf"
```

---

## ISSUE 4: EXPORT_RESULTS PROCESS (Line 129-147)

### Current Implementation
```nextflow
process EXPORT_RESULTS {
    input: 
    tuple val(sample_id), path(cnr), path(cns)
    
    script:
    """
    cnvkit.py export vcf $cns -i $sample_id -o ${sample_id}_CNVKIT_output.vcf
    cnvkit.py export bed $cns -i $sample_id -o ${sample_id}_calls.bed
    """
}
```

### ISSUE IDENTIFIED
✅ **CORRECT**: Uses `-i` flag for sample ID, which is correct

✅ **MOSTLY CORRECT**: The `-i` flag sets the SAMPLE column name in VCF

### POTENTIAL IMPROVEMENT
The current implementation will work IF the input `cns` file already has the "cn" column (from `call` command).

**To ensure correctness:**
```nextflow
process EXPORT_RESULTS {
    input: 
    tuple val(sample_id), path(cnr), path(call_cns)  // Note: call_cns from CALL_CNV output
    
    script:
    """
    # cnr file should be the raw copy number ratios
    # call_cns should be from 'cnvkit.py call' (has 'cn' column)
    
    cnvkit.py export vcf $call_cns -i "${sample_id}" -o ${sample_id}_CNVKIT_output.vcf
    cnvkit.py export bed $call_cns -i "${sample_id}" -o ${sample_id}_calls.bed
    """
}
```

### Optional: Add Filtering Before Export
For production pipelines, consider adding filters:
```nextflow
script:
"""
# Apply filters to remove uncertain segments
cnvkit.py call $cns --filter ci --merge bic -o ${sample_id}.filtered.cns

# Export filtered results
cnvkit.py export vcf ${sample_id}.filtered.cns -i "${sample_id}" \
  -o ${sample_id}_CNVKIT_output.vcf
cnvkit.py export bed ${sample_id}.filtered.cns -i "${sample_id}" \
  -o ${sample_id}_calls.bed
"""
```

---

## WORKFLOW INTEGRATION ISSUES

### Issue with COVERAGE→REFERENCE→CALL_CNV Chain

**Current workflow (Line 221-232):**
```nextflow
all_covs_ch = COVERAGE.out.target_cov
    .mix(COVERAGE.out.antitarget_cov)
    .map { sample_id, cov -> cov }
    .collect()

CREATE_POOLED_REFERENCE(fasta_ch, all_covs_ch)

cnv_input_ch = COVERAGE.out.target_cov
    .join(COVERAGE.out.antitarget_cov, by: 0)

CALL_CNV(cnv_input_ch, CREATE_POOLED_REFERENCE.out.ref_cnn)
```

✅ **This is CORRECT** - properly mixes and collects coverage files

---

## SUMMARY OF FIXES NEEDED

| Issue | Severity | Fix |
|-------|----------|-----|
| AUTOBIN takes all BAMs instead of median-sized one | 🔴 HIGH | Use only median-sized BAM or extract from filename/size |
| CREATE_POOLED_REFERENCE uses *.cnn (works but unclear) | 🟡 MEDIUM | Explicitly separate target/antitarget for clarity |
| CALL_CNV missing `cnvkit.py call` step | 🔴 CRITICAL | Add call command after segment |
| EXPORT_RESULTS assumes 'cn' column exists | 🟡 MEDIUM | Ensure input is from CALL_CNV output |
| No filtering before VCF export | 🟡 MEDIUM | Consider adding --filter ci --merge bic |

---

