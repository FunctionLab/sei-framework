#!/bin/bash
#SBATCH --time=2-00:00:00
#SBATCH --gres=gpu:1
#SBATCH --constraint=v100
#SBATCH --partition=gpu
#SBATCH -n 1
#SBATCH --mem 20000
#SBATCH -o o_run_vep_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<EMAIL>


vcf_filepath="${1:-}"
outdir="${2:-}"
reference_fa="${3:-}"


python -u vep_cli.py $vcf_filepath $outdir $reference_fa --cuda
