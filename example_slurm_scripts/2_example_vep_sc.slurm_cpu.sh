#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --partition=ccb


# Example SLURM script for getting the variant effect sequence class scores

ref_preds="${1:-}"    # ref preds path
alt_preds="${2:-}"    # alt preds path
outdir="${3:-}"       # output dir path
no_tsv="${4:-}"       # --no-tsv flag

sh ../2_varianteffect_sc_score.sh $ref_preds $alt_preds $outdir $no_tsv
