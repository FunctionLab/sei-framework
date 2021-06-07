# Sei framework

Welcome to the Sei framework repository! Sei is a framework for systematically predicting sequence regulatory activities and applying sequence information to human genetics data. Sei provides a global map from any sequence to regulatory activities, as represented by 40 sequence classes, and each sequence class integrates predictions for 21,907 chromatin profiles (transcription factor, histone marks, and chromatin accessibility profiles across a wide range of cell types).

This repository can be used to run the Sei model and get the Sei chromatin profile and sequence class predictions for an input VCF file.

We also provide information and instructions for [how to train the Sei deep learning sequence model](#train). 

### Requirements

Sei requires Python 3.6+ and Python packages PyTorch (>=1.0), Selene (>=0.5.0), and `docopt`. You can follow PyTorch installation steps [here](https://pytorch.org/get-started/locally/) and Selene installation steps [here](https://github.com/FunctionLab/selene). Install `docopt` with pip or conda (e.g. `conda install docopt`)

## Variant effect prediction

Sei predicts variant effects by comparing predictions from a pair of sequences carrying the reference allele and the alternative allele respectively. Our code provides both sequence-class level (40 classes) and chromatin profile-level (21,907 targets) predictions.

### Setup


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

The configuration file and script for running train is under the `train` directory. To run Sei deep learning sequence model training, you will need GPU computing capability (we run training on 4x Tesla V100 GPUs connected with NVLink). The training data can be downloaded by running the download script (`train/download_data.sh`) or manually downloading the training data tarball and extracting the files. **NOTE**: because the Sei training data contains processed files from the Cistrome Project, please first agree to the Cistrome Project [terms of usage](http://cistrome.org/db/#/bdown) before downloading the data here. 

The Sei training configuration YAML file is provided as the `train/train.yml` file. You can read more about the Selene command-line interface and configuration file formatting here (https://selene.flatironinstitute.org/master/overview/cli.html#). We also provide an example SLURM script `train.sh` for submitting a training job to a cluster.