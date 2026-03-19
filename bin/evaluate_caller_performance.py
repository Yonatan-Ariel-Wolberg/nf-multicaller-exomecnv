import pandas as pd

def load_bed_file(bed_file):
    return pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end'])

def load_cnv_file(cnv_file):
    return pd.read_csv(cnv_file, sep='\t', header=None, names=['chr', 'start', 'end', 'sample'])

def find_overlaps(probes, cnvs):
    overlaps = []
    for index, probe in probes.iterrows():
        for _, cnv in cnvs.iterrows():
            if (probe['chr'] == cnv['chr'] and 
                probe['start'] < cnv['end'] and 
                probe['end'] > cnv['start']):
                overlaps.append((probe, cnv))
    return overlaps

def categorize_probes(probes, truth_cnv, callset_cnv):
    categorized = {
        'TP': [],
        'FN': [],
        'FP': [],
        'TN': []
    }
    
    # Find overlaps with truth set
    truth_overlaps = find_overlaps(probes, truth_cnv)
    for probe in probes.itertuples(index=False):
        is_in_truth = any((probe.chr == overlap[0]['chr'] and 
                           probe.start < overlap[0]['end'] and 
                           probe.end > overlap[0]['start']) for overlap in truth_overlaps)
        
        is_called = any((probe.chr == cnv['chr'] and 
                         probe.start < cnv['end'] and 
                         probe.end > cnv['start']) for cnv in callset_cnv.itertuples(index=False))
        
        if is_in_truth and is_called:
            categorized['TP'].append(probe)
        elif is_in_truth and not is_called:
            categorized['FN'].append(probe)
        elif not is_in_truth and is_called:
            categorized['FP'].append(probe)
        else:
            categorized['TN'].append(probe)
    
    return categorized

def compute_metrics(categorized):
    TP = len(categorized['TP'])
    FN = len(categorized['FN'])
    FP = len(categorized['FP'])
    TN = len(categorized['TN'])

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    
    return sensitivity, precision, TN

def main(truth_bed, callset_bed, probe_bed):
    probes = load_bed_file(probe_bed)
    truth_cnv = load_cnv_file(truth_bed)
    callset_cnv = load_cnv_file(callset_bed)

    categorized = categorize_probes(probes, truth_cnv, callset_cnv)
    sensitivity, precision, TN = compute_metrics(categorized)

    print(f"Sensitivity: {sensitivity}")
    print(f"Precision: {precision}")
    print(f"True Negatives: {TN}")

# Example usage
# main('truthset_cnv.bed', 'callset_cnv.bed', 'probes_sanger.bed')
