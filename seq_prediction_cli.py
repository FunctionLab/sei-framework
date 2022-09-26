"""
Description:
    CLI for sequence prediction using the Sei deep learning model.
    Outputs Sei chromatin profile predictions.

Usage:
    seq_prediction_cli.py <seq-input> <output-dir> [--genome=<hg>] [--cuda]
    seq_prediction_cli.py -h | --help

Options:
    -h --help               Show this screen.
    <vcf>                   Input FASTA or BED file.
    <output-dir>            Output directory
    --genome=<hg>           If <seq-input> is a BED file, specify the reference
                              genome hg38 or hg19 [default: hg19]
    --cuda                  Run variant effect prediction on a CUDA-enabled
                              GPU

"""
import os

from docopt import docopt

from selene_sdk.sequences import Genome
from selene_sdk.utils import load_path
from selene_sdk.utils import parse_configs_and_run


def _finditem(obj, val):
    for k, v in obj.items():
        if hasattr(v, 'keywords'):
            _finditem(v.keywords, val)
        elif isinstance(v, dict):
            _finditem(v, val)
        elif isinstance(v, str) and '<PATH>' in v:
            obj[k] = v.replace('<PATH>', val)


if __name__ == "__main__":
    arguments = docopt(
        __doc__,
        version='0.0.0')

    seq_input = arguments["<seq-input>"]

    os.makedirs(arguments["<output-dir>"], exist_ok=True)

    # Assumes that the `models` directory is in the same directory as this
    # script. Please update this line if not.
    use_dir = os.path.dirname(os.path.abspath(__file__))
    use_cuda = arguments["--cuda"]

    sei_out = os.path.join(arguments["<output-dir>"], "chromatin-profiles-hdf5")
    os.makedirs(sei_out, exist_ok=True)

    configs = load_path("./model/sei_seq_prediction.yml", instantiate=False)
    _finditem(configs, use_dir)

    configs["prediction"].bind(input=seq_input, output_dir=sei_out)
    configs["analyze_sequences"].bind(use_cuda=use_cuda)

    if not seq_input.endswith('.fa') and not seq_input.endswith('.fasta'):
        hg_version = arguments["--genome"]
        genome = None
        if hg_version == 'hg38' or hg_version == 'hg19':
            genome = Genome(
                os.path.join('.', 'resources', '{0}_UCSC.fa'.format(hg_version)))
            configs["analyze_sequences"].bind(reference_sequence=genome)
        else:
            raise ValueError("--genome=<hg> must be 'hg19' or 'hg38'")

    parse_configs_and_run(configs)

