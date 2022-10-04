#!/usr/bin/env bash

#####################################################################
# Example script for computing the variant effect sequence class scores
# given Sei sequence predictions

# Usage:
# sh 2_varianteffect_sc_score.sh <ref-fp> <alt-fp> <output-dir>
#                                [--no-tsv]
#####################################################################

set -o errexit
set -o pipefail
set -o nounset

ref_fp="${1:-}"
alt_fp="${2:-}"
out_dir="${3:-}"
no_tsv="${4:-}"

mkdir -p $out_dir

echo "Input argments: $ref_fp $alt_fp $out_dir $no_tsv"

if [ "$no_tsv" = "--no-tsv" ]
then
    echo "--no-tsv flag is used"
    python -u 2_varianteffect_sc_score.py $ref_fp $alt_fp $out_dir --no-tsv
else
    echo "--no-tsv flag is not used"
    python -u 2_varianteffect_sc_score.py $ref_fp $alt_fp $out_dir
fi

