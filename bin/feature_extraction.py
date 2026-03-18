"""
Feature extraction for CNV machine-learning models.

Extracts a rich feature matrix from a SURVIVOR- or Truvari-merged SV VCF
together with per-caller normalised VCFs (produced by
``normalise_cnv_caller_quality_scores.py``) and optional genomic annotation
files (capture BED, BAM/CRAM, reference FASTA, mappability BED).

Features extracted
------------------
Structural
    chrom, start, end, cnv_size, cnv_type (1=DUP, 0=DEL/other)

Per-caller (one column per caller in ``tool_vcfs``)
    is_{caller}          -- 1 if caller supports this event, else 0
    qual_norm_{caller}   -- QUAL_norm score from the normalised VCF QUAL field
    raw_qual_{caller}    -- caller-native raw quality metric (non-redundant
                           with QUAL_norm): Q_SOME (CANOES/CLAMMS), SQ (XHMM),
                           QS (GATK-gCNV), CNQ (CNVkit), original QUAL (DRAGEN),
                           synthetic Phred score (INDELIBLE).

Concordance
    concordance -- number of callers supporting this event

Genomic annotations (NaN when the corresponding optional file is absent)
    n_probes     -- number of capture-target probes overlapping the CNV
    rd_ratio     -- mean read depth in CNV / mean read depth in flanking regions
    gc_content   -- GC fraction of the CNV interval
    mappability  -- weighted-mean mappability score over the CNV interval

Caller-specific secondary metrics (NaN when not present in the VCF)
    xhmm_rd        -- XHMM RD Z-score (INFO field)
    cnvkit_weight  -- CNVkit bin weight (INFO field)
    cnvkit_log2    -- CNVkit log2 ratio (INFO field)
    dragen_sm      -- DRAGEN sample median (INFO field)
    dragen_sd      -- DRAGEN sample SD (INFO field)

INDELIBLE split-read metrics (NaN when not matched in indelible_counts)
    total_sr, sr_entropy, mapq_avg, dual_split
"""

import pysam
import pandas as pd
import numpy as np
from collections import defaultdict


# ── BED / probe helpers ───────────────────────────────────────────────────────

def _load_bed(bed_file):
    """Return dict chrom -> list[(start, end)] loaded from a BED file."""
    intervals = defaultdict(list)
    with open(bed_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            intervals[chrom].append((start, end))
    return intervals


def _count_probes(chrom, start, end, bed_intervals):
    """Count BED target intervals that overlap the half-open interval [start, end)."""
    count = 0
    for iv_start, iv_end in bed_intervals.get(chrom, []):
        if iv_start < end and iv_end > start:
            count += 1
    return count


# ── Read-depth helpers ────────────────────────────────────────────────────────

def _mean_depth(bam, chrom, start, end):
    """Compute mean read depth over [start, end) using pysam pileup."""
    if start >= end:
        return 0.0
    depths = [
        col.nsegments
        for col in bam.pileup(chrom, start, end, truncate=True, min_base_quality=0)
    ]
    return float(np.mean(depths)) if depths else 0.0


def _rd_ratio(bam, chrom, start, end, flank=500):
    """Read-depth ratio: mean depth in CNV region / mean depth in flanking regions.

    Returns NaN when flanking depth is zero (avoids division by zero).
    """
    target_depth = _mean_depth(bam, chrom, start, end)
    left_start = max(0, start - flank)
    left_depth = _mean_depth(bam, chrom, left_start, start)
    right_depth = _mean_depth(bam, chrom, end, end + flank)
    flank_depth = (left_depth + right_depth) / 2.0
    if flank_depth == 0.0:
        return np.nan
    return target_depth / flank_depth


# ── GC content ────────────────────────────────────────────────────────────────

def _gc_content(fasta, chrom, start, end):
    """Fraction of G/C bases in the genomic interval [start, end)."""
    try:
        seq = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        return np.nan
    if not seq:
        return np.nan
    gc = sum(1 for c in seq if c in 'GC')
    return gc / len(seq)


# ── Mappability helpers ───────────────────────────────────────────────────────

def _load_mappability_bed(mappability_file):
    """Load a 4-column mappability BED (chrom, start, end, score)."""
    intervals = defaultdict(list)
    with open(mappability_file) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start, end, score = int(parts[1]), int(parts[2]), float(parts[3])
            intervals[chrom].append((start, end, score))
    return intervals


def _mean_mappability(chrom, start, end, mappability_intervals):
    """Weighted-mean mappability score over [start, end) from a mappability BED.

    Returns NaN when no intervals overlap or the mappability dict is empty.
    """
    if not mappability_intervals:
        return np.nan
    total_bases = 0
    weighted_sum = 0.0
    for iv_start, iv_end, score in mappability_intervals.get(chrom, []):
        ov_start = max(iv_start, start)
        ov_end = min(iv_end, end)
        if ov_end > ov_start:
            bases = ov_end - ov_start
            total_bases += bases
            weighted_sum += score * bases
    if total_bases == 0:
        return np.nan
    return weighted_sum / total_bases


# ── Per-caller raw quality extractor ─────────────────────────────────────────

def _get_raw_qual(record, caller):
    """Extract the caller-native quality metric from a (normalised) VCF record.

    Normalised VCFs produced by ``normalise_cnv_caller_quality_scores.py``
    retain all original FORMAT/INFO fields, so the raw scores are still
    accessible alongside the new QUAL_norm value in the QUAL column.

    Returns (metric_name: str, value: float) or (None, np.nan) if unavailable.

    Caller mappings
    ---------------
    canoes    -- Q_SOME FORMAT field
    clamms    -- Q_SOME FORMAT field
    xhmm      -- SQ FORMAT field
    gatk      -- QS FORMAT field
    cnvkit    -- CNQ FORMAT field
    dragen    -- original QUAL preserved in OQ FORMAT field by normalisation
    indelible -- synthetic Phred = SR_TOTAL * (AVG_MAPQ / 60.0) * 100  (INFO)
    """
    sample = next(iter(record.samples.values()), None)
    try:
        if caller == 'canoes':
            if sample is not None and 'Q_SOME' in sample:
                return 'Q_SOME', float(sample['Q_SOME'])
        elif caller == 'clamms':
            if sample is not None and 'Q_SOME' in sample:
                return 'Q_SOME', float(sample['Q_SOME'])
        elif caller == 'xhmm':
            if sample is not None and 'SQ' in sample:
                return 'SQ', float(sample['SQ'])
        elif caller == 'gatk':
            if sample is not None and 'QS' in sample:
                return 'QS', float(sample['QS'])
        elif caller == 'cnvkit':
            if sample is not None and 'CNQ' in sample:
                return 'CNQ', float(sample['CNQ'])
        elif caller == 'dragen':
            # After normalisation the original QUAL is stored in OQ FORMAT field.
            if sample is not None and 'OQ' in sample:
                return 'QUAL', float(sample['OQ'])
        elif caller == 'indelible':
            # Reconstruct the synthetic score used as the QUAL_norm input.
            if 'SR_TOTAL' in record.info and 'AVG_MAPQ' in record.info:
                sr_total = float(record.info['SR_TOTAL'])
                avg_mapq = float(record.info['AVG_MAPQ'])
                synthetic = sr_total * (avg_mapq / 60.0) * 100.0
                return 'SYNTHETIC_PHRED', synthetic
    except (ValueError, TypeError):
        pass
    return None, np.nan


# ── Main extraction function ──────────────────────────────────────────────────

def extract_normalized_features(
    merged_vcf,
    tool_vcfs,
    indelible_counts,
    merger_mode='survivor',
    sample_id=None,
    bed_file=None,
    bam_file=None,
    reference_fasta=None,
    mappability_file=None,
    rd_flank=500,
):
    """Extract a feature matrix from a merged SV VCF for ML model training/scoring.

    Parameters
    ----------
    merged_vcf : str
        Path to SURVIVOR- or Truvari-merged VCF.
    tool_vcfs : dict[str, str]
        Ordered mapping of caller name -> path to *normalised* per-caller VCF
        (produced by ``normalise_cnv_caller_quality_scores.py``).
        In SURVIVOR mode the dict order must match the SUPP_VEC bit positions.
    indelible_counts : pd.DataFrame
        INDELIBLE small-variant count table with columns:
        Start, Total_SR, Entropy, MAPQ_Avg, Dual_Split.
    merger_mode : {'survivor', 'truvari'}
        How the merged VCF was produced.
    sample_id : str, optional
        Sample identifier (not currently used for feature values but reserved
        for downstream labelling).
    bed_file : str, optional
        BED file of capture target regions used to count probes per CNV.
    bam_file : str, optional
        BAM or CRAM alignment file used to compute read-depth ratio.
    reference_fasta : str, optional
        Indexed FASTA reference used for GC-content calculation and (when
        bam_file is a CRAM) as the CRAM reference.
    mappability_file : str, optional
        4-column BED file (chrom, start, end, score) with mappability scores.
    rd_flank : int
        Flanking bases on each side used for RD-ratio calculation (default 500).

    Returns
    -------
    pd.DataFrame
        One row per merged SV record; columns described in the module docstring.
    """
    # ── load optional annotation sources ─────────────────────────────────
    bed_intervals = _load_bed(bed_file) if bed_file else {}
    mappability_intervals = (
        _load_mappability_bed(mappability_file) if mappability_file else {}
    )

    bam = None
    if bam_file:
        if bam_file.endswith('.cram'):
            bam = pysam.AlignmentFile(
                bam_file, 'rc', reference_filename=reference_fasta
            )
        else:
            bam = pysam.AlignmentFile(bam_file, 'rb')

    fasta = pysam.FastaFile(reference_fasta) if reference_fasta else None

    # ── open VCFs ─────────────────────────────────────────────────────────
    vcf_in = pysam.VariantFile(merged_vcf)
    tools = {k: pysam.VariantFile(v) for k, v in tool_vcfs.items()}
    caller_order = list(tool_vcfs.keys())

    all_records = []

    for record in vcf_in:
        v_data = {
            'chrom':    record.chrom,
            'start':    record.pos,
            'end':      record.stop,
            'cnv_size': record.stop - record.pos,
            'cnv_type': 1 if 'DUP' in str(record.info.get('SVTYPE', '')) else 0,
        }

        # ── GC content ────────────────────────────────────────────────────
        v_data['gc_content'] = (
            _gc_content(fasta, record.chrom, record.pos, record.stop)
            if fasta is not None
            else np.nan
        )

        # ── probe count ───────────────────────────────────────────────────
        v_data['n_probes'] = (
            _count_probes(record.chrom, record.pos, record.stop, bed_intervals)
            if bed_intervals
            else np.nan
        )

        # ── mappability ───────────────────────────────────────────────────
        v_data['mappability'] = _mean_mappability(
            record.chrom, record.pos, record.stop, mappability_intervals
        )

        # ── RD ratio ──────────────────────────────────────────────────────
        v_data['rd_ratio'] = (
            _rd_ratio(bam, record.chrom, record.pos, record.stop, rd_flank)
            if bam is not None
            else np.nan
        )

        # ── per-caller quality features ───────────────────────────────────
        if merger_mode == 'survivor':
            supp_vec = str(record.info.get('SUPP_VEC', '0' * len(caller_order)))
            v_data['concordance'] = sum(int(x) for x in supp_vec)

            for i, bit in enumerate(supp_vec):
                caller_name = caller_order[i] if i < len(caller_order) else f'tool_{i}'
                v_data[f'is_{caller_name}'] = int(bit)

                if bit == '1' and caller_name in tools:
                    matches = list(
                        tools[caller_name].fetch(
                            record.chrom, record.pos - 10, record.pos + 10
                        )
                    )
                    if matches:
                        orig = matches[0]
                        # QUAL_norm: normalised QUAL written by
                        # normalise_cnv_caller_quality_scores.py
                        v_data[f'qual_norm_{caller_name}'] = (
                            orig.qual if orig.qual is not None else np.nan
                        )
                        # Raw (non-redundant) caller-native quality metric
                        _, raw_val = _get_raw_qual(orig, caller_name)
                        v_data[f'raw_qual_{caller_name}'] = raw_val

                        # Caller-specific secondary INFO metrics
                        if caller_name == 'xhmm' and 'RD' in orig.info:
                            v_data['xhmm_rd'] = orig.info['RD']
                        if caller_name == 'cnvkit':
                            if 'weight' in orig.info:
                                v_data['cnvkit_weight'] = orig.info['weight']
                            if 'log2' in orig.info:
                                v_data['cnvkit_log2'] = orig.info['log2']
                        if caller_name == 'dragen':
                            if 'SM' in orig.info:
                                v_data['dragen_sm'] = orig.info['SM']
                            if 'SD' in orig.info:
                                v_data['dragen_sd'] = orig.info['SD']
                    else:
                        v_data[f'qual_norm_{caller_name}'] = np.nan
                        v_data[f'raw_qual_{caller_name}'] = np.nan
                else:
                    v_data[f'qual_norm_{caller_name}'] = np.nan
                    v_data[f'raw_qual_{caller_name}'] = np.nan

        elif merger_mode == 'truvari':
            var_id = str(record.id) if record.id else ''
            for caller_name in caller_order:
                is_supp = 1 if caller_name.upper() in var_id.upper() else 0
                v_data[f'is_{caller_name}'] = is_supp

                if is_supp and caller_name in tools:
                    matches = list(
                        tools[caller_name].fetch(
                            record.chrom, record.pos - 10, record.pos + 10
                        )
                    )
                    if matches:
                        orig = matches[0]
                        v_data[f'qual_norm_{caller_name}'] = (
                            orig.qual if orig.qual is not None else np.nan
                        )
                        _, raw_val = _get_raw_qual(orig, caller_name)
                        v_data[f'raw_qual_{caller_name}'] = raw_val

                        if caller_name == 'xhmm' and 'RD' in orig.info:
                            v_data['xhmm_rd'] = orig.info['RD']
                        if caller_name == 'cnvkit':
                            if 'weight' in orig.info:
                                v_data['cnvkit_weight'] = orig.info['weight']
                            if 'log2' in orig.info:
                                v_data['cnvkit_log2'] = orig.info['log2']
                        if caller_name == 'dragen':
                            if 'SM' in orig.info:
                                v_data['dragen_sm'] = orig.info['SM']
                            if 'SD' in orig.info:
                                v_data['dragen_sd'] = orig.info['SD']
                    else:
                        v_data[f'qual_norm_{caller_name}'] = np.nan
                        v_data[f'raw_qual_{caller_name}'] = np.nan
                else:
                    v_data[f'qual_norm_{caller_name}'] = np.nan
                    v_data[f'raw_qual_{caller_name}'] = np.nan

        # ── INDELIBLE split-read counts ───────────────────────────────────
        indelible_data = indelible_counts[indelible_counts['Start'] == record.pos]
        if not indelible_data.empty:
            v_data['total_sr']   = indelible_data['Total_SR'].values[0]
            v_data['sr_entropy'] = indelible_data['Entropy'].values[0]
            v_data['mapq_avg']   = indelible_data['MAPQ_Avg'].values[0]
            v_data['dual_split'] = indelible_data['Dual_Split'].values[0]
        else:
            v_data['total_sr']   = np.nan
            v_data['sr_entropy'] = np.nan
            v_data['mapq_avg']   = np.nan
            v_data['dual_split'] = np.nan

        all_records.append(v_data)

    # ── cleanup ───────────────────────────────────────────────────────────
    if bam is not None:
        bam.close()
    if fasta is not None:
        fasta.close()

    return pd.DataFrame(all_records)


# Example usage:
# tool_vcfs = {                               # order matches SURVIVOR SUPP_VEC bits
#     'canoes':    'sample_CANOES.normalised.vcf.gz',
#     'clamms':    'sample_CLAMMS.normalised.vcf.gz',
#     'xhmm':      'sample_XHMM.normalised.vcf.gz',
#     'gatk':      'sample_GCNV.normalised.vcf.gz',
#     'cnvkit':    'sample_CNVKIT.normalised.vcf.gz',
#     'dragen':    'sample_DRAGEN.normalised.vcf.gz',
#     'indelible': 'sample_INDELIBLE.normalised.vcf.gz',
# }
# indelible_counts = pd.read_csv('indelible_counts.tsv', sep='\t')
# df = extract_normalized_features(
#     'merged.vcf',
#     tool_vcfs,
#     indelible_counts,
#     merger_mode='survivor',
#     bed_file='capture_targets.bed',
#     bam_file='sample.cram',
#     reference_fasta='GRCh38.fa',
#     mappability_file='mappability_100mer.bed',
# )
