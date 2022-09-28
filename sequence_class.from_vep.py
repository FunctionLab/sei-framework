"""
Description:
    This script is run after `vep_cli.py`. It computes the sequence class-level
    variant effect scores based on Sei chromatin profile predictions,
    and by default, sorts the variants based on the maximum absolute scores
    across sequence classes and outputs the results as TSV files.

Usage:
    sequence_class.from_vep.py <ref-fp> <alt-fp> <output-dir>
                               [--out-name=<out-name>]
                               [--no-tsv]
    sequence_class.from_vep.py -h | --help

Options:
    <ref-fp>               Reference sequence Sei predictions file. Assumes
                           the row labels file is in the same directory. The resulting
                           variant effect sequence class scores will be computed as
                           alt - ref sequence class scores (adjusted for nucleosome
                           occupancy).
    <alt-fp>               Alternate sequence Sei predictions file.
    <output-dir>           The directory to output the sequence class scores
                           & TSVs (if `--no-tsv` is not used).
    --out-name=<out-name>  Specify an output filename prefix that all outputted
                           files will use. Otherwise, filenames will be based on
                           <alt-fp>.
    --no-tsv               The TSVs outputted sort the variants based on maximum
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
from utils import sc_hnorm_varianteffect
from utils import write_to_tsv


if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='1.0.0')
    ref_pred_file = arguments['<ref-fp>']
    alt_pred_file = arguments['<alt-fp>']

    # load predictions
    chromatin_profile_ref = get_data(ref_pred_file)
    chromatin_profile_alt = get_data(alt_pred_file)
    if len(chromatin_profile_ref) != len(chromatin_profile_alt):
        raise ValueError(("{0} and {1} have different number of rows: {2} vs {3}, "
                          "respectively.").format(ref_pred_file, alt_pred_file,
                                                  len(chromatin_profile_ref),
                                                  len(chromatin_profile_alt)))

    ref_dir, _ = os.path.split(ref_pred_file)
    alt_dir, _ = os.path.split(alt_pred_file)

    # checks if the ref/alt are from variant effect prediction (VCF)
    # or sequence prediction (BED or FASTA file inputs)
    seq_from = None
    alt_prefix = None
    if '.ref_predictions' in ref_pred_file and '.alt_predictions' in alt_pred_file:
        seq_from = 'VCF'
        alt_prefix = os.path.basename(alt_pred_file).split('.alt_predictions')[0]
    elif '.alt_predictions' in ref_pred_file and 'ref_predictions' in ref_pred_file:
        seq_from = 'VCF'
        alt_prefix = os.path.basename(alt_pred_file).split('.ref_predictions')[0]
    elif '_predictions' in ref_pred_file and '_predictions' in alt_pred_file:
        seq_from = 'BED/FASTA'
        alt_prefix = os.path.basename(alt_pred_file).split('_predictions')[0]
    rowlabels_file = os.path.join(alt_dir, '{0}_row_labels.txt'.format(alt_prefix))
    rowlabels = pd.read_csv(rowlabels_file, sep='\t')
    if len(rowlabels) != len(chromatin_profile_alt):
        raise ValueError(("Rowlabels file '{0}' does not have the same number "
                          "of rows as '{1}'").format(rowlabels_file, alt_pred_file))

    output_dir = arguments['<output-dir>']
    os.makedirs(output_dir, exist_ok=True)

    output_prefix = arguments['--out-name']
    if output_prefix is None:
        output_prefix = alt_prefix
    print("Output files will start with {0}".format(output_prefix))

    no_tsv = arguments['--no-tsv']

    sei_dir = "./model"
    chromatin_profiles = get_targets(os.path.join(sei_dir, "target.names"))
    seqclass_names = get_targets(os.path.join(sei_dir, "seqclass.names"))

    clustervfeat = np.load(os.path.join(sei_dir, 'projvec_targets.npy'))
    histone_inds = np.load(os.path.join(sei_dir, 'histone_inds.npy'))

    diffproj = sc_hnorm_varianteffect(
        chromatin_profile_ref,
        chromatin_profile_alt,
        clustervfeat,
        histone_inds)
    max_abs_diff = np.abs(diffproj).max(axis=1)

    np.save(os.path.join(
        output_dir, "{0}.sequence_class_scores.npy".format(output_prefix)), diffproj)

    if not no_tsv:
        write_to_tsv(max_abs_diff,  # max sequence class score
                     chromatin_profile_alt - chromatin_profile_ref,  # chromatin profile diffs
                     diffproj,  # sequence class diffs
                     chromatin_profiles,  # chromatin profile targets
                     seqclass_names,  # sequence class names
                     rowlabels,
                     os.path.join(output_dir,
                                  "sorted.{0}.chromatin_profile_diffs.tsv".format(
                                      output_prefix)),
                     os.path.join(output_dir,
                                  "sorted.{0}.sequence_class_scores.tsv".format(
                                      output_prefix)))

