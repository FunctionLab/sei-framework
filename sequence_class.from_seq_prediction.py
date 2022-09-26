"""
Description:
    This script is run after `seq_prediction_cli.py`. It computes the
    sequence class-level scores based on Sei chromatin profile predictions,
    and by default, sorts the variants based on the maximum absolute scores
    across sequence classes and outputs the results as TSV files.

Usage:
    sequence_class.from_seq_prediction.py <results-dir>
                                          [--input=<input-preds>]
                                          [--ref=<ref-preds>] [--alt=<alt-preds>]
                                          [--no-tsv]
    sequence_class.from_seq_prediction.py -h | --help

Options:
    <results-dir>           The results directory.
    --input=<input-preds>   Path to sequence predictions to compute sequence
                            class projection scores (NO variant effect).
                            Use `--input` ONLY or `--ref` and `--alt` ONLY.
                            If `--input` is used `--ref` and `--alt` will be ignored.
    --ref=<ref-preds>       Path to reference sequence predictions to compute
                            sequence class variant effect scores.
                            `--alt` must also be used.
    --alt=<alt-preds>       Path to the alternate sequence predictions to
                            compute sequence class variant effect scores.
                            `--ref` must also be used.
    --no-tsv                The TSVs outputted sort the variants based on maximum
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
    input_pred = arguments['--input']
    ref_pred = arguments['--ref']
    alt_pred = arguments['--alt']
    no_tsv = arguments['--no-tsv']

    if input_pred is None and ref_pred is None:

    print(arguments)
    import sys
    sys.exit(0)

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

