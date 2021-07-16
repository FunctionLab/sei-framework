#!/bin/sh

# Trained Sei model
wget https://zenodo.org/record/4906997/files/sei_model.tar.gz

tar -xzvf sei_model.tar.gz

# Sei framework resources (FASTA files)
wget https://zenodo.org/record/4906962/files/sei_framework_resources.tar.gz

tar -xzvf sei_framework_resources.tar.gz
