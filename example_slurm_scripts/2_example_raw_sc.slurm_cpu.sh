#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --partition=ccb

# Example SLURM script for getting the raw sequence class scores
# for input sequence predictions (i.e. no variants)

input_preds="${1:-}"    # input preds path
outdir="${2:-}"         # output dir path

sh ./2_raw_sc_score.sh $input_preds $outdir
