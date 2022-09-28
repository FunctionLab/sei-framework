"""
Description:
    This script is run after `seq_prediction_cli.py`. It computes the
    sequence class-level scores based on Sei chromatin profile predictions,
    and by default, sorts the variants based on the maximum absolute scores
    across sequence classes and outputs the results as TSV files.

Usage:
    get_raw_sc_score.py <input-fp> <output-dir>
                        [--out-name=<out-name>]
                        [--no-tsv]
    get_raw_sc_score.py -h | --help

Options:
    <input-fp>             Path to Sei sequence predictions to compute raw
                           sequence class projection scores (NO variant effect).
    <output-dir>           The directory to output the sequence class scores
                           & TSVs (if `--no-tsv` is not used).
    --out-name=<out-name>  Specify an output filename prefix that all outputted
                           files will use. Otherwise, filenames will be based on
                           <input-fp>.
    --no-tsv               The TSVs outputted sort the variants based on maximum
                           absolute scores across sequence classes and are intended
                           for easier perusal of the predictions for those more
                           familiar with this file type, but makes this script take
                           much longer to complete as a result. If you are comfortable
                           working with HDF5 and NPY files, you can suppress the TSV
                           output and use the files in `chromatin-profiles-hdf5`.
    --no-chromprof-tsv     TODO
"""
import os

from docopt import docopt
import h5py
import numpy as np
import pandas as pd

from utils import get_filename_prefix, get_data, get_targets
from utils import sc_projection
from utils import write_to_tsv


if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    input_pred_file = arguments['<input-fp>']
    input_preds = get_data(input_pred_file)
    input_dir, input_fn = os.path.split(input_pred_file)

    input_prefix = get_filename_prefix(input_fn)
    rowlabels_file = os.path.join(input_dir, '{0}_row_labels.txt'.format(
        input_prefix))
    rowlabels = pd.read_csv(rowlabels_file, sep='\t')
    if len(rowlabels) != len(input_preds):
        raise ValueError(("Rowlabels file '{0}' does not have the same number "
                          "of rows as '{1}'").format(rowlabels_file,
                                                     input_pred_file))

    output_prefix = arguments['--out-name']
    if output_prefix is None:
        output_prefix = input_prefix
    print("Output files will start with {0}".format(output_prefix))

    no_tsv = arguments['--no-tsv']

    sei_dir = "./model"
    chromatin_profiles = get_targets(os.path.join(sei_dir, "target.names"))
    seqclass_names = get_targets(os.path.join(sei_dir, "seqclass.names"))

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))

    projscores = sc_projection(input_preds, clustervfeat)
    np.save(os.path.join(
        output, "{0}.raw_sequence_class_scores.npy".format(output_prefix)),
        projscores)

    if not no_tsv:
        write_to_tsv(max_abs_diff,  # max sequence class score
                     chromatin_profile_alt - chromatin_profile_ref,  # chromatin profile diffs
                     diffproj,  # sequence class diffs
                     chromatin_profiles,  # chromatin profile targets
                     seqclass_names,  # sequence class names
                     read_rowlabels_file(chromatin_profile_rowlabels, use_strand=False),
                     os.path.join(results_dir, "sorted.chromatin_profile_diffs.tsv"),
                     os.path.join(results_dir, "sorted.sequence_class_scores.tsv"))

