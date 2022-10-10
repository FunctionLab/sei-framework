#!/usr/bin/env bash

#####################################################################
# Example script for running Sei deep learning model sequence
# prediction with Selene

# Usage:
# sh 1_sequence_prediction.sh <input-file> <genome> <output-dir> --cuda

# Please only specify hg38 or hg19 as input for <genome> if <input-file> is
# a BED file. If you are using a FASTA file of sequences, you can specify
# whatever genome version you wish for your own reference (will be printed
# as part of output for this script) but do not leave it empty.

# --cuda is optional, use if you are running on a CUDA-enabled
# GPU machine (see example_slurm_scripts/1_example_seqpred.slurm_gpu.sh)
#####################################################################

set -o errexit
set -o pipefail
set -o nounset

input_filepath="${1:-}"
hg_version="${2:-}"
out_dir="${3:-}"
cuda="${4:-}"

mkdir -p $out_dir

input_basename=$(basename $input_filepath)
cp $input_filepath $out_dir/

echo "Input argments: $input_filepath $out_dir $hg_version $cuda"

if [ "$cuda" = "--cuda" ]
then
    echo "use_cuda: True"
    python -u 1_sequence_prediction.py \
        "$out_dir/$input_basename" \
        $out_dir \
        --genome=${hg_version} \
        --cuda
else
    echo "use_cuda: False"
    python -u 1_sequence_prediction.py \
        "$out_dir/$input_basename" \
        $out_dir \
        --genome=${hg_version}
fi
