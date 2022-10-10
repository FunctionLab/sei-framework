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
    return (np.dot(chromatin_profile_preds, clustervfeat.T) /
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


def get_filename_prefix(filename):
    """Filename must follow Selene output file conventions.
    """
    prefix = None
    if '.alt_predictions' in filename:
        prefix = filename.split('.alt_predictions')[0]
    elif '.ref_predictions' in filename:
        prefix = filename.split('.ref_predictions')[0]
    else:
        prefix = filename.split('_predictions')[0]
    return prefix


def write_to_tsv(max_abs_diff,
                 chromatin_profile_diffs,
                 sequence_class_projscores,
                 chromatin_profiles,
                 seqclass_names,
                 rowlabels,
                 output_chromatin_profile_file,
                 output_sequence_class_file):
    sorted_sc_abs_diff = np.sort(max_abs_diff)[::-1]
    sorted_maxsc_df = pd.DataFrame(sorted_sc_abs_diff,  # dataframe
                                   columns=['seqclass_max_absdiff'])

    sorted_ixs = np.argsort(max_abs_diff)[::-1]
    assert len(sorted_ixs) == chromatin_profile_diffs.shape[0]
    sorted_rowlabels = rowlabels.iloc[sorted_ixs]  # dataframe
    sorted_rowlabels.reset_index(inplace=True)
    rowlabel_columns = sorted_rowlabels.columns.tolist()

    # sorted now
    chromatin_profile_diffs = chromatin_profile_diffs[sorted_ixs, :]
    sequence_class_projscores = sequence_class_projscores[sorted_ixs,:]

    # dataframes
    sorted_profiles_df = pd.DataFrame(chromatin_profile_diffs, columns=chromatin_profiles)
    sorted_sc_df = pd.DataFrame(sequence_class_projscores, columns=seqclass_names)
    del chromatin_profile_diffs

    sei_df = pd.concat([sorted_maxsc_df, sorted_rowlabels, sorted_profiles_df],
                       axis=1)
    sc_df = pd.concat([sorted_maxsc_df, sorted_rowlabels, sorted_sc_df],
                      axis=1)
    sc_df[['seqclass_max_absdiff'] + rowlabel_columns + seqclass_names].to_csv(
        output_sequence_class_file, sep='\t', index=False)

    if len(sei_df) > 10000:
        sei_df[['seqclass_max_absdiff'] + rowlabel_columns + chromatin_profiles].to_csv(
            output_chromatin_profile_file, sep='\t', index=False, compression='gzip')
    else:
        sei_df[['seqclass_max_absdiff'] + rowlabel_columns + chromatin_profiles].to_csv(
            output_chromatin_profile_file, sep='\t', index=False)



