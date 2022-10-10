#!/usr/bin/env bash

#####################################################################
# Example script for computing the raw sequence class scores
# given Sei chromatin profile sequence predictions

# Usage:
# sh 2_raw_sc_score.sh <input-file> <output-dir>
#####################################################################

set -o errexit
set -o pipefail
set -o nounset

input_filepath="${1:-}"
out_dir="${2:-}"

mkdir -p $out_dir

echo "Input argments: $input_filepath $out_dir"

python -u 2_raw_sc_score.py $input_filepath $out_dir
