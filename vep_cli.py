"""
Description:
    CLI for variant effect prediction using the Sei deep learning model.
    Outputs both Sei chromatin profile predictions and sequence class scores.

Usage:
    vep_cli.py <vcf> <output-dir> [--genome=<hg>] [--cuda]
    vep_cli.py -h | --help

Options:
    -h --help               Show this screen.
    <vcf>                   Input VCF file
    <output-dir>            Output directory
    --genome=<hg>           hg38 or hg19 [default: hg19]
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

    os.makedirs(arguments["<output-dir>"], exist_ok=True)

    # Assumes that the `models` directory is in the same directory as this
    # script. Please update this line if not.
    use_dir = os.path.dirname(os.path.abspath(__file__))

    hg_version = arguments["--genome"]
    genome = None
    if hg_version == 'hg38' or hg_version == 'hg19':
        genome = Genome(
            os.path.join('.', 'resources', '{0}_UCSC.fa'.format(hg_version)))
    else:
        raise ValueError("--genome=<hg> must be 'hg19' or 'hg38'")

    use_cuda = arguments["--cuda"]

    def run_config(config_yml, output_dir):
        configs = load_path(config_yml, instantiate=False)
        _finditem(configs, use_dir)
        configs["analyze_sequences"].bind(
            reference_sequence=genome,
            use_cuda=use_cuda)
        configs["variant_effect_prediction"].update(
            vcf_files=[arguments["<vcf>"]],
            output_dir=output_dir)
        parse_configs_and_run(configs)

    sei_out = os.path.join(arguments["<output-dir>"], "chromatin-profiles-hdf5")
    os.makedirs(sei_out, exist_ok=True)
    run_config("./model/sei_prediction.yml", sei_out)

