#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --gpus=v100-32gb:4
#SBATCH --partition=gpu
#SBATCH --constraint=v100
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --mem 100000




python3 -u ../selene_sdk/cli.py train.yml --lr=0.1
