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
# on input sequences from a BED or FASTA file

input_filepath="${1:-}"  # path to input bed/fasta file
genome_version="${2:-}"  # genome version
outdir="${3:-}"          # path to output dir

sh ../1_sequence_prediction.sh $input_filepath $genome_version $outdir --cuda
