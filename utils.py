import os

import h5py
import numpy as np
import pandas as pd


def get_targets(filename):
    """
    Load the names of all the chromatin profiles predicted by Sei
    """
    targets = []
    with open(filename, 'r') as file_handle:
        for line in file_handle:
            targets.append(line.strip())
    return targets


def get_data(filename):
    """
    Load HDF5 file of predictions into memory
    """
    fh = h5py.File(filename, 'r')
    data = fh["data"][()]
    fh.close()
    return data


def sc_projection(chromatin_profile_preds, clustervfeat):
    return (np.dot(chromatin_profile_pred, clustervfeat.T) /
            np.linalg.norm(clustervfeat, axis=1))


def sc_hnorm_varianteffect(chromatin_profile_ref, chromatin_profile_alt, clustervfeat, histone_inds):
    chromatin_profile_ref_adjust = chromatin_profile_ref.copy()
    chromatin_profile_ref_adjust[:, histone_inds] = \
        chromatin_profile_ref_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_ref[:, histone_inds], axis=1))[:, None]

    chromatin_profile_alt_adjust = chromatin_profile_alt.copy()
    chromatin_profile_alt_adjust[:, histone_inds] = \
        chromatin_profile_alt_adjust[:, histone_inds] * (
        (np.sum(chromatin_profile_ref[:, histone_inds], axis=1)*0.5 +
         np.sum(chromatin_profile_alt[:, histone_inds], axis=1)*0.5) /
        np.sum(chromatin_profile_alt[:, histone_inds], axis=1))[:, None]

    refproj = sc_projection(chromatin_profile_ref_adjust, clustervfeat)
    altproj = sc_projection(chromatin_profile_alt_adjust, clustervfeat)
    diffproj = altproj[:,:40] - refproj[:,:40]
    return diffproj


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



