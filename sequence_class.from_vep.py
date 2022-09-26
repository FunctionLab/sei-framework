"""
Description:
    This script is run after `vep_cli.py`. It computes the sequence class-level
    variant effect scores based on Sei chromatin profile predictions,
    and by default, sorts the variants based on the maximum absolute scores
    across sequence classes and outputs the results as TSV files.

Usage:
    sequence_class.from_vep.py <results-dir> <vcf> [--no-tsv]
    sequence_class.from_vep.py -h | --help

Options:
    <results-dir>    The results directory.
    <vcf>            Name of VCF file.
    --no-tsv         The TSVs outputted sort the variants based on maximum
                     absolute scores across sequence classes and are intended
                     for easier perusal of the predictions for those more
                     familiar with this file type, but makes this script take
                     much longer to complete as a result. If you are comfortable
                     working with HDF5 and NPY files, you can suppress the TSV
                     output and use the files in `chromatin-profiles-hdf5`.

"""
import os

from docopt import docopt
import h5py
import numpy as np
import pandas as pd

from utils import get_targets, get_data
from utils import read_rowlabels_file
from utils import sc_hnorm_varianteffect
from utils import write_to_tsv


if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    results_dir = arguments['<results-dir>']
    input_vcf = arguments['<vcf>']
    no_tsv = arguments['--no-tsv']

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

    diffproj = sc_hnorm_varianteffect(
        chromatin_profile_ref,
        chromatin_profile_alt,
        clustervfeat,
        histone_inds)
    max_abs_diff = np.abs(diffproj).max(axis=1)

    np.save(os.path.join(profile_pred_dir, "sequence_class_scores.npy"), diffproj)
    if not no_tsv:
        write_to_tsv(max_abs_diff,  # max sequence class score
                     chromatin_profile_alt - chromatin_profile_ref,  # chromatin profile diffs
                     diffproj,  # sequence class diffs
                     chromatin_profiles,  # chromatin profile targets
                     seqclass_names,  # sequence class names
                     read_rowlabels_file(chromatin_profile_rowlabels, use_strand=False),
                     os.path.join(results_dir, "sorted.chromatin_profile_diffs.tsv"),
                     os.path.join(results_dir, "sorted.sequence_class_scores.tsv"))

