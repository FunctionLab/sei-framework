#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --constraint=v100
#SBATCH --partition=gpu
#SBATCH -n 1
#SBATCH --mem 100G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<EMAIL>

# Example SLURM script for running the Sei framework code
# on variants in the input VCF file

vcf_filepath="${1:-}"  # path to vcf file
hg_version="${2:-}"    # hg19 or hg38
outdir="${3:-}"        # path to output dir

sh ./1_variant_effect_prediction.sh $vcf_filepath $hg_version $outdir --cuda
