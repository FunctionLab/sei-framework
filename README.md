# Sei framework

Welcome to the Sei framework repository! Sei is a framework for systematically predicting sequence regulatory activities and applying sequence information to human genetics data. Sei provides a global map from any sequence to regulatory activities, as represented by 40 sequence classes, and each sequence class integrates predictions for 21,907 chromatin profiles (transcription factor, histone marks, and chromatin accessibility profiles across a wide range of cell types).

This repository can be used to run the Sei model and get the Sei chromatin profile and sequence class predictions for an input VCF file.

We also provide information and instructions for [how to train the Sei deep learning sequence model](#training). 

### Requirements

Sei requires Python 3.6+ and Python packages PyTorch (>=1.0), Selene (>=0.5.0), and `docopt`. You can follow PyTorch installation steps [here](https://pytorch.org/get-started/locally/) and Selene installation steps [here](https://github.com/FunctionLab/selene). Install `docopt` with pip or conda (e.g. `conda install docopt`)

## Variant effect prediction

Sei predicts variant effects by comparing predictions from a pair of sequences carrying the reference allele and the alternative allele respectively. Our code provides both sequence-class level (40 classes) and chromatin profile-level (21,907 targets) predictions.

### Setup

Please download and extract the trained Sei model and `resources` (containing hg19 and hg38 FASTA files) `.tar.gz` files before proceeding:

```
sh ./download_data.sh
```

- [Sei model](https://doi.org/10.5281/zenodo.4906996)
- [Sei framework `resources` directory](https://doi.org/10.5281/zenodo.4906961)

### Usage

```
sh run_pipeline.sh <example-vcf> <hg-version> <output-dir> [--cuda]
```
Example command 
```
sh run_pipeline.sh test.vcf hg19 test-output --cuda
```

Arguments:
- `<example-vcf>`: Input VCF file
- `<hg-version>`: Reference FASTA. By default the framework accepts `hg19` or `hg38` coordinates.
- `<output-dir>`: Path to output directory (will be created if does not exist)
- `--cuda`: Optional, use this flag if running on a CUDA-enabled GPU.

We provide `test.vcf` (hg19 coordinates) so you can try running this command once you have installed all the requirements. Additionally, `run_pipeline.gpu_node.sh` is an example SLURM script with the same expected input arguments if you need to submit your job to a compute cluster. 

**Additional note**: we have added the capability of predicting variant effects from a pair of sequences in the `vep_cli_seq.py` script in the [`vep_seq`](https://github.com/FunctionLab/sei-framework/tree/vep_seq) development branch of this repo.

### Outputs

The following files and directories will be outputted:
-  `chromatin-profiles-hdf5`: directory
-  `chromatin_profiles_diffs.tsv`: chromatin profile prediction TSV file (**Note:** output file will be compressed if input has >10000 variants)
-  `sequence_classes_scores.tsv`: sequence class prediction TSV file 

The two `*.tsv` files are the final formatted outputs, while the `chromatin-profiles-hdf5` directory contains the intermediate HDF5 and row label files outputted from Selene from running the Sei deep learning model. 

You can use the HDF5 files directly if desired, but please keep in mind that the variants will not be ordered in the same way as the TSV files. (Please see the corresponding `*_row_labels.txt` file, for the variant labels.) 



### Sequence classes

Sequence classes are defined based on 30 million sequences tiling the genome and thus cover a wide range of sequence activities. To help interpretation, we grouped sequence classes into groups including P (Promoter), E (Enhancer), CTCF (CTCF-cohesin binding), TF (TF binding), PC (Polycomb-repressed), HET (Heterochromatin), TN (Transcription), and L (Low Signal) sequence classes. Please refer to our manuscript for a more detailed description of the sequence classes.


| Sequence class label |               Sequence class name | Rank by size | Group |
|---------------------:|----------------------------------:|-------------:|------:|
|                 PC1  |       Polycomb / Heterochromatin  |            0 |   PC  |
|                  L1  |                       Low signal  |            1 |    L  |
|                 TN1  |                    Transcription  |            2 |   TN  |
|                 TN2  |                    Transcription  |            3 |   TN  |
|                  L2  |                       Low signal  |            4 |    L  |
|                  E1  |                        Stem cell  |            5 |    E  |
|                  E2  |                     Multi-tissue  |            6 |    E  |
|                  E3  |               Brain / Melanocyte  |            7 |    E  |
|                  L3  |                       Low signal  |            8 |    L  |
|                  E4  |                     Multi-tissue  |            9 |    E  |
|                 TF1  |                    NANOG / FOXA1  |           10 |   TF  |
|                 HET1 |                  Heterochromatin  |           11 |  HET  |
|                  E5  |                      B-cell-like  |           12 |    E  |
|                  E6  |                  Weak epithelial  |           13 |    E  |
|                 TF2  |                            CEBPB  |           14 |   TF  |
|                 PC2  |                    Weak Polycomb  |           15 |   PC  |
|                  E7  |            Monocyte / Macrophage  |           16 |    E  |
|                  E8  |                Weak multi-tissue  |           17 |    E  |
|                  L4  |                       Low signal  |           18 |    L  |
|                 TF3  |                FOXA1 / AR / ESR1  |           19 |   TF  |
|                 PC3  |                         Polycomb  |           20 |   PC  |
|                 TN3  |                    Transcription  |           21 |   TN  |
|                  L5  |                       Low signal  |           22 |    L  |
|                 HET2 |                  Heterochromatin  |           23 |  HET  |
|                  L6  |                       Low signal  |           24 |    L  |
|                   P  |                         Promoter  |           25 |    P  |
|                  E9  |                Liver / Intestine  |           26 |    E  |
|                 CTCF |                     CTCF-Cohesin  |           27 |  CTCF |
|                 TN4  |                    Transcription  |           28 |   TN  |
|                 HET3 |                  Heterochromatin  |           29 |  HET  |
|                 E10  |                            Brain  |           30 |    E  |
|                 TF4  |                             OTX2  |           31 |   TF  |
|                 HET4 |                  Heterochromatin  |           32 |  HET  |
|                  L7  |                       Low signal  |           33 |    L  |
|                 PC4  | Polycomb / Bivalent stem cell Enh |           34 |   PC  |
|                 HET5 |                       Centromere  |           35 |  HET  |
|                 E11  |                           T-cell  |           36 |    E  |
|                 TF5  |                               AR  |           37 |   TF  |
|                 E12  |                Erythroblast-like  |           38 |    E  |
|                 HET6 |                       Centromere  |           39 |   HET |


## Training

The configuration file and script for running train is under the `train` directory. To run Sei deep learning sequence model training, you will need GPU computing capability (we run training on 4x Tesla V100 GPUs connected with NVLink). 

The training data is available [here](https://doi.org/10.5281/zenodo.4907037) should be downloaded and extracted into the `train` directory. **NOTE**: because the Sei training data contains processed files from the Cistrome Project, please first agree to the Cistrome Project [terms of usage](http://cistrome.org/db/#/bdown) before downloading the data:

```
cd ./train
sh ./download_data.sh  # in the train directory
```

The Sei training configuration YAML file is provided as the `train/train.yml` file. You can read more about the Selene command-line interface and configuration file formatting [here](https://selene.flatironinstitute.org/master/overview/cli.html#). You must use Selene version >0.5.0 to train this model ([release notes](https://github.com/FunctionLab/selene/blob/master/RELEASE_NOTES.md)). 

We also provide an example SLURM script `train.sh` for submitting a training job to a cluster.

## Help 
Please post in the Github issues or e-mail Kathy Chen (chen.kathleenm@gmail.com) with any questions about the repository, requests for more data, etc. 

## License

If you are interested in obtaining the software for commercial use, please contact Office of Technology Licensing, Princeton University (Laurie Tzodikov 609-258-7256, tzodikov@princeton.edu). 

```
Copyright (c) [2021] [The Trustees of Princeton University, The Simons Foundation, Inc. and The University of Texas Southwestern Medical Center]
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted for academic and research use only (subject to the limitations in the disclaimer below) provided that the following conditions are met:
     * Redistributions of source code must retain the above copyright notice,
     this list of conditions and the following disclaimer.
     * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer in the
     documentation and/or other materials provided with the distribution.
     * Neither the name of the copyright holders nor the names of its
     contributors may be used to endorse or promote products derived from this
     software without specific prior written permission.
NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
THIS LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```
