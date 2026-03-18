import pysam
import pandas as pd
import numpy as np
import subprocess

def normalize_q(score, tool):
    """Normalizes various tool quality scores to a 0-1 scale."""
    if score is None or np.isnan(score): 
        return 0.0
    if tool in ['xhmm', 'gatk', 'canoes', 'clamms']:
        return min(float(score) / 100.0, 1.0)
    if tool == 'dragen':
        return min(float(score) / 60.0, 1.0)
    if tool == 'cnvkit':
        return max(1.0 - float(score), 0.0)
    if tool == 'indelible':
        return float(score)  # Already 0-1 probability
    return 0.0

def extract_normalized_features(merged_vcf, tool_vcfs, indelible_counts, merger_mode='survivor', sample_id=None):
    vcf_in = pysam.VariantFile(merged_vcf)
    tools = {k: pysam.VariantFile(v) for k, v in tool_vcfs.items()}
    
    all_records = []  # Initialize an empty list to store results

    for record in vcf_in:
        v_data = {
            'chrom': record.chrom,
            'start': record.pos,
            'end': record.stop,
            'size': record.stop - record.pos,
            'sv_type': 1 if 'DUP' in str(record.info.get('SVTYPE', '')) else 0
        }
        
        caller_scores = []  # Initialize caller_scores for the current record

        # Extract different metrics based on the merger mode
        if merger_mode == 'survivor':
            supp_vec = str(record.info.get('SUPP_VEC', '0000000'))
            v_data['concordance'] = sum(int(x) for x in supp_vec)
            for i, bit in enumerate(supp_vec):
                v_data[f'is_tool_{i}'] = int(bit)
                if bit == '1' and i in tools:
                    matches = list(tools[i].fetch(record.chrom, record.pos - 10, record.pos + 10))
                    if matches:
                        orig = matches[0]  # Using the first match for simplification
                        # Mapping indices to tool names for normalization
                        tool_map = {0: 'xhmm', 1: 'gatk', 4: 'cnvkit', 5: 'dragen', 6: 'indelible'}
                        score = orig.info.get('SQ', orig.qual)  # Get the quality score
                        caller_scores.append(normalize_q(score, tool_map.get(i, 'other')))
                        
            # Extract additional features 
            if 'RD' in orig.info:  # XHMM RD (Z-score)
                v_data['xhmm_rd'] = orig.info['RD']

            if 'QS' in orig.info:  # GATK QS
                v_data['gatk_qs'] = orig.info['QS']

            if 'CNQ' in orig.info:  # GATK CNQ
                v_data['gatk_cnq'] = orig.info['CNQ']

            if 'weight' in orig.info:  # CNVkit Weight
                v_data['cnvkit_weight'] = orig.info['weight']

            if 'log2' in orig.info:  # CNVkit Log2
                v_data['cnvkit_log2'] = orig.info['log2']

            # DRAGEN related features
            if 'SM' in orig.info:  # DRAGEN SM
                v_data['dragen_sm'] = orig.info['SM']

            if 'SD' in orig.info:  # DRAGEN SD
                v_data['dragen_sd'] = orig.info['SD']

        elif merger_mode == 'truvari':
            var_id = str(record.id) if record.id else ""
            for name, t_vcf in tool_vcfs.items():
                is_supp = 1 if name.upper() in var_id.upper() else 0
                v_data[f'is_{name}'] = is_supp
                if is_supp:
                    matches = list(tools[name].fetch(record.chrom, record.pos - 10, record.pos + 10))
                    if matches:
                        caller_scores.append(normalize_q(matches[0].qual, name))

        # Handle indelible counts for small variants
        indelible_data = indelible_counts[indelible_counts['Start'] == record.start]
        if not indelible_data.empty:
            v_data['total_sr'] = indelible_data['Total_SR'].values[0]
            v_data['sr_entropy'] = indelible_data['Entropy'].values[0]
            v_data['mapq_avg'] = indelible_data['MAPQ_Avg'].values[0]
            v_data['dual_split'] = indelible_data['Dual_Split'].values[0]

        v_data['score_qual'] = max(caller_scores) if caller_scores else 0.0
        all_records.append(v_data)

    return pd.DataFrame(all_records)

# Example usage:
# tool_vcfs = {
#     'xhmm': 'path_to_xhmm.vcf',
#     'gatk': 'path_to_gatk.vcf',
#     'canoes': 'path_to_canoes.vcf',
#     'clamms': 'path_to_clamms.vcf',
#     'cnvkit': 'path_to_cnvkit.vcf',
#     'dragen': 'path_to_dragen.vcf',
#     'indelible': 'path_to_indelible.vcf'
# }
# indelible_counts = pd.read_csv('path_to_indelible_counts.csv')  # Load Indelible counts data
# df_features = extract_normalized_features('merged_vcf.vcf', tool_vcfs, indelible_counts, 'survivor')
