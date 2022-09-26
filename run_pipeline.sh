#!/usr/bin/env bash

###############################################################
# Example script for running the variant effect prediction
# pipeline for Sei and sequence classes.

# sh run_pipeline.sh <vcf> <hg> <output-dir> --cuda

# Please only specify hg38 or hg19 as input for <hg>.

# --cuda is optional, use if you are running on a CUDA-enabled
# GPU machine
###############################################################

set -o errexit
set -o pipefail
set -o nounset

vcf_filepath="${1:-}"
hg_version="${2:-}"
outdir="${3:-}"
cuda="${4:-}"

mkdir -p $outdir

vcf_basename=$(basename $vcf_filepath)
cp $vcf_filepath $outdir/

if [ "$cuda" = "--cuda" ]
then
    echo "use_cuda: True"
    python -u vep_cli.py "$outdir/$vcf_basename" \
                      $outdir \
                      --genome=${hg_version} \
                      --cuda
else
    echo "use_cuda: False"
    python -u vep_cli.py "$outdir/$vcf_basename" \
                      $outdir \
                      --genome=${hg_version}
fi

python sequence_class.py $outdir "$vcf_basename"
