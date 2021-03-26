"""
Description:
    This script is run after `vep_cli.py`. It computes the disease impact
    score for each variant based on DeepSEA and SeqWeaver predicted
    difference scores, sorts the variants based on the disease impact scores,
    and outputs the results as TSV files.

Usage:
    sequence_class_projscores_computation.py <results-dir> <vcf>
    sequence_class_projscores_computation.py -h | --help

Options:
    <results-dir>    The results directory. Final outputs will be written to
                     a directory in <results-dir> called `final_outputs`.
                     The results directory was an input to `vep_cli.py` and
                     should contain the `sei`, `seqweaver_human`, and
                     `seqweaver_mouse` directories created from running
                     `vep_cli.py`.
    <vcf>            Name of VCF file. If you ran
                     `annotate_vcf_with_gene_info.py`, use the resulting VCF
                     from running that script.

"""
import os

from docopt import docopt
import h5py
import numpy as np
import pandas as pd


def get_features(filename):
    features = []
    with open(filename, 'r') as file_handle:
        for line in file_handle:
            features.append(line.strip())
    return features


def get_data(filename):
    fh = h5py.File(filename, 'r')
    data = fh["data"][()]
    fh.close()
    return np.abs(data)




def read_rowlabels_file(rowlabels, use_strand=False):
    labels = []
    with open(rowlabels, 'r') as file_handle:
        for line in file_handle:
            if 'contains_unk' in line:
                continue
            cols = line.strip().split('\t')
            chrom, pos, id, ref, alt, strand, ref_match, contains_unk = cols
            if not use_strand:
                strand = '.'
            labels.append((chrom, pos, id, ref, alt, strand, ref_match, contains_unk))
    return labels


def write_to_tsv(sorted_ixs, data, sequence_class_projscores, features, seqclass_features, rowlabels, output_file, output_sequence_class_file):
    assert len(sorted_ixs) == data.shape[0]
    rowlabels = np.array(rowlabels)[sorted_ixs]
    data = data[sorted_ixs, :]
    sequence_class_projscores = sequence_class_projscores[sorted_ixs,:]
    rows = []
    rows_data = []
    for ix, r in enumerate(rowlabels):
        ref_match = r[-2]
        contains_unk = r[-1]
        label = r[:-2]
        row = np.hstack([[ref_match, contains_unk], label, sequence_class_projscores[ix,:]]).tolist()
        row_data = np.hstack([[ref_match, contains_unk], label, data[ix, :]]).tolist()
        rows.append(row)
        rows_data.append(row_data)
    df = pd.DataFrame(rows, columns=  [
        'ref_match', 'contains_unk',
        'chrom', 'pos', 'id', 'ref', 'alt', 'strand',] + seqclass_features)
    check_columns = df.columns.tolist()
    check_columns.remove('strand')
    df.drop_duplicates(subset=check_columns, inplace=True)
    df.to_csv(output_sequence_class_file, sep='\t', index=False)
    
    df = pd.DataFrame(rows_data, columns=  [
        'ref_match', 'contains_unk',
        'chrom', 'pos', 'id', 'ref', 'alt', 'strand',] + features)
    check_columns = df.columns.tolist()
    check_columns.remove('strand')
    df.drop_duplicates(subset=check_columns, inplace=True)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    results_dir = arguments['<results-dir>']
    input_vcf = arguments['<vcf>']

    _, filename = os.path.split(input_vcf)
    filename_prefix = '.'.join(filename.split('.')[:-1])

    sei_dir = os.path.join("models", "sei")
    sei_features = get_features(os.path.join(
        sei_dir, "target.names"))

    seqclass_features = get_features(os.path.join(
        sei_dir, "seqclass.names"))
        
    rowlabels_filename = "{0}_row_labels.txt".format(filename_prefix)
    sei_rowlabels = os.path.join(
        results_dir, "sei", rowlabels_filename)
        
    sei_data_ref = get_data(os.path.join(results_dir, "sei", "{0}.ref_predictions.h5".format(filename_prefix)))
    sei_data_alt = get_data(os.path.join(results_dir, "sei", "{0}.alt_predictions.h5".format(filename_prefix)))

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_features.npy'))
    histone_inds = np.load(os.path.join(sei_dir, 'histone_inds.npy'))
    

    sei_data_ref_adjust = sei_data_ref.copy()
    sei_data_ref_adjust[:,histone_inds] = sei_data_ref_adjust[:,histone_inds] * ( (np.sum(sei_data_ref[:, histone_inds], axis=1)*0.5 + np.sum(sei_data_alt[:, histone_inds], axis=1)*0.5) / np.sum(sei_data_ref[:, histone_inds], axis=1))[:, None]

    sei_data_alt_adjust = sei_data_alt.copy()
    sei_data_alt_adjust[:,histone_inds] = sei_data_alt_adjust[:,histone_inds] * ((np.sum(sei_data_ref[:, histone_inds], axis=1)*0.5 + np.sum(sei_data_alt[:, histone_inds], axis=1)*0.5) / np.sum(sei_data_alt[:, histone_inds], axis=1))[:, None]


    refproj =  np.dot((sei_data_ref_adjust),clustervfeat.T)/np.linalg.norm(clustervfeat,axis=1)
    altproj =  np.dot((sei_data_alt_adjust),clustervfeat.T)/np.linalg.norm(clustervfeat,axis=1)
    diffproj = altproj[:,:40] - refproj[:,:40]

    max_abs_diff = np.abs(diffproj).max(axis=1)
    sorted_sei_ixs = np.argsort(max_abs_diff)[::-1]
    output_dir = os.path.join(results_dir, "final_outputs")
    os.makedirs(output_dir, exist_ok=True)

    write_to_tsv(sorted_sei_ixs,
                 sei_data_alt - sei_data_ref,
                 diffproj,
                 sei_features,
                 seqclass_features,
                 read_rowlabels_file(sei_rowlabels, use_strand=False),
                 os.path.join(output_dir, "sei_results.tsv"),
                 os.path.join(output_dir, "sequenceclass_results.tsv"))


