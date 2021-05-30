"""
Description:
    This script is run after `vep_cli.py`. It computes the sequence-class
    level variant effect scores based on Sei prediction, sorts the variants
    based on the maximum absolute scores across sequence classes and outputs
    the results as TSV files.

Usage:
    sequence_class.py <results-dir> <vcf>
    sequence_class.py -h | --help

Options:
    <results-dir>    The results directory.
    <vcf>            Name of VCF file.

"""
import os

from docopt import docopt
import h5py
import numpy as np
import pandas as pd


def get_targets(filename):
    targets = []
    with open(filename, 'r') as file_handle:
        for line in file_handle:
            targets.append(line.strip())
    return targets


def get_data(filename):
    fh = h5py.File(filename, 'r')
    data = fh["data"][()]
    fh.close()
    return data


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
            labels.append(
                (chrom, pos, id, ref, alt, strand,
                 ref_match, contains_unk))
    return labels


def write_to_tsv(max_abs_diff,
                 chromatin_profile_diffs,
                 sequence_class_projscores,
                 chromatin_profiles,
                 seqclass_names,
                 rowlabels,
                 output_chromatin_profile_file,
                 output_sequence_class_file):
    sorted_sc_abs_diff = np.sort(max_abs_diff)[::-1]
    sorted_ixs = np.argsort(max_abs_diff)[::-1]

    assert len(sorted_ixs) == chromatin_profile_diffs.shape[0]
    rowlabels = np.array(rowlabels)[sorted_ixs]
    chromatin_profile_diffs = chromatin_profile_diffs[sorted_ixs, :]
    sequence_class_projscores = sequence_class_projscores[sorted_ixs,:]
    chromatin_profile_rows = []
    sequence_class_rows = []
    for ix, r in enumerate(rowlabels):
        ref_match = r[-2]
        contains_unk = r[-1]
        label = r[:-2]
        chromatin_profile_row = np.hstack(
            [[sorted_sc_abs_diff[ix], ref_match, contains_unk],
             label,
             chromatin_profile_diffs[ix, :]]).tolist()
        sequence_class_row = np.hstack(
            [[sorted_sc_abs_diff[ix], ref_match, contains_unk],
             label,
             sequence_class_projscores[ix, :]]).tolist()
        chromatin_profile_rows.append(chromatin_profile_row)
        sequence_class_rows.append(sequence_class_row)

    del chromatin_profile_diffs
    colnames = [
        'seqclass_max_absdiff', 'ref_match', 'contains_unk',
        'chrom', 'pos', 'id', 'ref', 'alt', 'strand']

    sc_df = pd.DataFrame(
        sequence_class_rows, columns=colnames + seqclass_names)
    check_columns = sc_df.columns.tolist()
    check_columns.remove('strand')
    sc_df.drop_duplicates(subset=check_columns, inplace=True)
    sc_df.to_csv(output_sequence_class_file, sep='\t', index=False)

    cp_df = pd.DataFrame(
        chromatin_profile_rows, columns=colnames + chromatin_profiles)
    check_columns = cp_df.columns.tolist()
    check_columns.remove('strand')
    cp_df.drop_duplicates(subset=check_columns, inplace=True)
    if len(cp_df) > 10000:
        cp_df.to_csv(
            output_chromatin_profile_file, sep='\t', index=False, compression='gzip')
    else:
        cp_df.to_csv(
            output_chromatin_profile_file, sep='\t', index=False)



if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    results_dir = arguments['<results-dir>']
    input_vcf = arguments['<vcf>']

    _, filename = os.path.split(input_vcf)
    filename_prefix = '.'.join(filename.split('.')[:-1])

    sei_dir = "./model"
    chromatin_profiles = get_targets(os.path.join(sei_dir, "target.names"))
    seqclass_names = get_targets(os.path.join(sei_dir, "seqclass.names"))

    profile_pred_dir = os.path.join(results_dir, 'chromatin-profiles-hdf5')
    rowlabels_filename = "{0}_row_labels.txt".format(filename_prefix)
    chromatin_profile_rowlabels = os.path.join(profile_pred_dir, rowlabels_filename)

    chromatin_profile_ref = get_data(os.path.join(
        profile_pred_dir,
        "{0}.ref_predictions.h5".format(filename_prefix)))
    chromatin_profile_alt = get_data(os.path.join(
        profile_pred_dir,
        "{0}.alt_predictions.h5".format(filename_prefix)))

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))
    histone_inds = np.load(os.path.join(sei_dir, 'histone_inds.npy'))

    chromatin_profile_ref_adjust = chromatin_profile_ref.copy()
    chromatin_profile_ref_adjust[:,histone_inds] = \
        chromatin_profile_ref_adjust[:,histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_ref[:, histone_inds], axis=1))[:, None]

    chromatin_profile_alt_adjust = chromatin_profile_alt.copy()
    chromatin_profile_alt_adjust[:,histone_inds] = \
        chromatin_profile_alt_adjust[:,histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_alt[:, histone_inds], axis=1))[:, None]

    refproj = (np.dot(chromatin_profile_ref_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    altproj = (np.dot(chromatin_profile_alt_adjust,clustervfeat.T) /
               np.linalg.norm(clustervfeat,axis=1))
    diffproj = altproj[:,:40] - refproj[:,:40]

    max_abs_diff = np.abs(diffproj).max(axis=1)

    write_to_tsv(max_abs_diff,  # max sequence class score
                 chromatin_profile_alt - chromatin_profile_ref,  # chromatin profile diffs
                 diffproj,  # sequence class diffs
                 chromatin_profiles,  # chromatin profile targets
                 seqclass_names,  # sequence class names
                 read_rowlabels_file(chromatin_profile_rowlabels, use_strand=False),
                 os.path.join(results_dir, "chromatin_profile_diffs.tsv"),
                 os.path.join(results_dir, "sequence_class_scores.tsv"))


