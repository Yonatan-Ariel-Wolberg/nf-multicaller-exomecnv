# CNVKit Pipeline - Quick Reference Guide

## Question 1: Cohort-Only Workflow
**For exome sequencing with tumor/proband samples ONLY (no normals):**

✅ **BEST:** Pooled reference from cohort
```bash
# Calculate coverage for all samples
for sample in *.bam; do
  cnvkit.py coverage $sample targets.bed -o ${sample%.bam}.targetcoverage.cnn
  cnvkit.py coverage $sample antitargets.bed -o ${sample%.bam}.antitargetcoverage.cnn
done

# Build reference from all samples' coverage
cnvkit.py reference *.cnn --fasta hg38.fa -o pooled_ref.cnn
```

✅ **ALTERNATIVE:** Flat reference
```bash
cnvkit.py reference -o flat_ref.cnn -f hg38.fa -t targets.bed -a antitargets.bed
```

---

## Question 2: autobin Command

| Aspect | Details |
|--------|---------|
| **Purpose** | Estimate optimal bin sizes |
| **Input** | Multiple BAM files |
| **Key Behavior** | **Uses ONLY median-sized BAM** |
| **Outputs** | targets.bed, antitargets.bed |
| **Command** | `cnvkit.py autobin *.bam -t targets.bed -g access.bed` |

**⚠️ CRITICAL:** Does NOT process all BAMs - selects median-sized one

---

## Question 3: reference Command

| Aspect | Details |
|--------|---------|
| **Purpose** | Build reference from coverage files |
| **Input Files** | **BOTH** targetcoverage.cnn + antitargetcoverage.cnn |
| **Pairing** | Auto-paired by sample ID in filename |
| **Command** | `cnvkit.py reference *.cnn --fasta hg38.fa -o ref.cnn` |
| **Output Columns** | log2, spread, gc, rmask |

**✅ CORRECT:** All .cnn files (auto-pairs by name)

---

## Question 4: segment vs call

```
Input BAM
    ↓
[coverage] → *.targetcoverage.cnn + *.antitargetcoverage.cnn
    ↓
[fix] → *.cnr (log2 ratios per bin)
    ↓
[segment] → *.cns (breakpoints, log2 values)
    ↓
[call] ← CRITICAL STEP! → *.call.cns (adds "cn" column)
    ↓
[export vcf] → *.vcf (requires "cn" column!)
```

| Command | Input | Output | Purpose |
|---------|-------|--------|---------|
| **segment** | .cnr | .cns (log2) | Find breakpoints |
| **call** | .cns | .cns + "cn" | Convert to integers |

---

## Question 5: export vcf Flags

```bash
cnvkit.py export vcf Sample.call.cns -i "SampleID" -x female -o Sample.vcf
                    └─────────────────┘          └──────────┘
                    MUST have "cn" column        Sample name
```

| Flag | Purpose | Example |
|------|---------|---------|
| `-i` | Sample ID in VCF | `-i "Patient123"` |
| `-x` | Sample sex | `-x female` or `-x male` |
| `-y` | Haploid X reference | `-y` (if used in reference building) |

**⚠️ REQUIREMENT:** Input .cns MUST have "cn" column (from `call` command)

---

## Question 6: Complete Workflow

```bash
# PHASE 1: Setup (once)
cnvkit.py access hg38.fa -o access.hg38.bed
cnvkit.py autobin *.bam -t baits.bed -g access.hg38.bed
cnvkit.py target baits.bed -o targets.bed
cnvkit.py antitarget targets.bed -o antitargets.bed

# PHASE 2: Coverage (parallelizable)
for sample in *.bam; do
  cnvkit.py coverage $sample targets.bed -o ${sample%.bam}.targetcoverage.cnn
  cnvkit.py coverage $sample antitargets.bed -o ${sample%.bam}.antitargetcoverage.cnn
done

# PHASE 3: Reference (once)
cnvkit.py reference *.cnn --fasta hg38.fa -o pooled_ref.cnn

# PHASE 4: Per-sample calling (parallelizable)
for sample in *.bam; do
  id=${sample%.bam}
  cnvkit.py fix ${id}.targetcoverage.cnn ${id}.antitargetcoverage.cnn pooled_ref.cnn -o ${id}.cnr
  cnvkit.py segment ${id}.cnr -o ${id}.cns
  cnvkit.py call ${id}.cns -o ${id}.call.cns                          # ← CRITICAL!
  cnvkit.py export vcf ${id}.call.cns -i "$id" -o ${id}.cnv.vcf       # ← REQUIRES call
done
```

---

## Module Issues Checklist

### ✅ Current Issues Found

| Process | Issue | Severity | Status |
|---------|-------|----------|--------|
| AUTOBIN | Passes all BAMs (should use median-sized) | 🔴 HIGH | FIX NEEDED |
| CREATE_POOLED_REFERENCE | Uses *.cnn (works but unclear) | 🟡 MEDIUM | OK (could improve clarity) |
| CALL_CNV | **MISSING `cnvkit.py call`** | 🔴 CRITICAL | **FIX NEEDED** |
| EXPORT_RESULTS | Assumes "cn" column exists | 🟡 MEDIUM | Depends on CALL_CNV fix |

### 🔴 CRITICAL FIX REQUIRED

**In CALL_CNV process (Line 107-127):**

```nextflow
// CURRENT (WRONG)
script:
"""
cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.cns -o ...
cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.cns -o ...
"""

// CORRECT
script:
"""
cnvkit.py fix $t_cov $t_anticov $reference -o ${sample_id}.cnr
cnvkit.py segment ${sample_id}.cnr -o ${sample_id}.cns
cnvkit.py call ${sample_id}.cns -o ${sample_id}.call.cns    // ← ADD THIS!
cnvkit.py scatter ${sample_id}.cnr -s ${sample_id}.call.cns -o ...
cnvkit.py diagram ${sample_id}.cnr -s ${sample_id}.call.cns -o ...
"""
```

---

## Key Sources

- **Repository:** https://github.com/etal/cnvkit
- **Files Analyzed:**
  - `/cnvlib/batch.py` - Source code for batch workflow
  - `/doc/pipeline.rst` - Command documentation
  - `/doc/calling.rst` - Calling and filtering
  - `/doc/importexport.rst` - Export options

- **Key Documentation Quote:**
  > "you can create a pooled reference even if matched tumor-normal pairs were sequenced -- our benchmarking showed that a pooled reference performed slightly better"

---

## Quick Decision Tree

```
Do you have normal/control samples?
├─ YES → Use pooled reference from normals (best practice)
└─ NO  → Choose:
         ├─ Flat reference (simpler, less accurate)
         └─ Pooled reference from cohort samples (better, recommended)

Need VCF output?
├─ YES → Must run: fix → segment → CALL → export vcf
└─ NO  → Can skip call, but needed for integer copy numbers

Multiple BAM files?
├─ YES for autobin → CNVKit uses median-sized BAM
└─ Single BAM → Uses that BAM
```

