# Sei framework

Welcome to the Sei framework repository! Sei is a framework for systematically predicting sequence regulatory activities and applying sequence information to human genetics data. Sei provides a global map from any sequence to regulatory activities, as represented by 40 sequence classes, and each sequence class integrates predictions for 21,907 chromatin profiles (transcription factor, histone marks, and chromatin accessibility profiles across a wide range of cell types).

This repository can be used to run the Sei model and get the Sei chromatin profile and sequence class predictions for an input VCF file.

We also provide information and instructions for [how to train the Sei deep learning sequence model](#training). 

### Requirements

Sei requires Python 3.6+ and Python packages PyTorch (>=1.0), Selene (>=0.5.0), and `docopt`. You can follow PyTorch installation steps [here](https://pytorch.org/get-started/locally/) and Selene installation steps [here](https://github.com/FunctionLab/selene). Install `docopt` with pip or conda (e.g. `conda install docopt`)

## Variant effect prediction

Sei predicts variant effects by comparing predictions from a pair of sequences carrying the reference allele and the alternative allele respectively. Our code provides both sequence-class level (40 classes) and chromatin profile-level (21,907 targets) predictions.

### Setup

Please download and extract the trained Sei model and `resources` (containing hg19 and hg38 FASTA files) `.tar.gz` files before proceeding.  

- [Sei model](https://www.dropbox.com/s/4q4kixk4roxw3v2/sei_model.tar.gz?dl=0)
- [Sei framework `resources` directory](https://www.dropbox.com/s/wam5tg6g3gpor5w/sei_framework_resources.tar.gz?dl=0)

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

The training data is available [here](https://www.dropbox.com/s/lonq7o8uwft7sbg/sei_training_data.tar.gz?dl=0) should be downloaded and extracted into the `train` directory. **NOTE**: because the Sei training data contains processed files from the Cistrome Project, please first agree to the Cistrome Project [terms of usage](http://cistrome.org/db/#/bdown) before downloading the data. 

The Sei training configuration YAML file is provided as the `train/train.yml` file. You can read more about the Selene command-line interface and configuration file formatting [here](https://selene.flatironinstitute.org/master/overview/cli.html#). We also provide an example SLURM script `train.sh` for submitting a training job to a cluster.

## Use Agreement

If you are interested in obtaining the software for commercial use, please contact Office of Technology Licensing, Princeton University (Laurie Tzodikov 609-258-7256, tzodikov@princeton.edu, or Linda Jan, 609-258-3653, ljan@princeton.edu). For academic use, downloading or using the software means you agree with the following Academic Use SOFTWARE Agreement.

```
PRINCETON Academic Use SOFTWARE Agreement

The Trustees of Princeton University, a non-profit educational corporation organized and existing under the
laws of the State of New Jersey with its Office of Technology Licensing at 87 Prospect Avenue, Princeton NJ
08544 (“PRINCETON”) is willing to make the Sei Framework Software (software for predicting chromatin profile and sequence class effects from sequence) (“SOFTWARE”) available to you
(“RECIPIENT”) under the following terms and conditions: 

1. The above SOFTWARE is the property of PRINCETON and is made available as a service to the research
community. The SOFTWARE is a research tool still in the development stage and is being provided “as is”
without any support, services or improvements. PRINCETON makes no representations and extends no
warranties of any kind, and expressly disclaims any representations or warranties of any kind
(including, but not limited to, any warranties of merchantability, fitness for a particular purpose, or
non-infringement).
2. The SOFTWARE will be used for teaching or not-for-profit research purposes only. The RECIPIENT agrees
not to use the Software for commercial purposes, or for diagnosis, treatment, cure, prevention or mitigation of
disease or any other conditions in man. The RECIPIENT agrees that the Software is not intended to substitute
for care by a licensed healthcare professional.
    3. The SOFTWARE will not be further distributed to others without PRINCETON’S prior written consent.
    The RECIPIENT shall promptly refer any request for the SOFTWARE to PRINCETON.
    4. RECIPIENT
    acknowledges that any programs or software created based on the SOFTWARE will be considered a
    derivative of SOFTWARE and owned by PRINCETON.
    5. The RECIPIENT agrees to acknowledge the source of the SOFTWARE in any publications.
    6. RECIPIENT will use the SOFTWARE in compliance with all applicable laws, policies and regulations
    including, but not limited to, any approvals, informed consent, patient confidentiality principles and US
    government and local export control regulations. RECIPIENT acknowledges that the SOFTWARE may
    not be exported to Cuba, Iran, North Korea, or Syria.
    7. In no event shall PRINCETON be responsible or liable to RECIPIENT or any third party for
    RECIPIENT’s activities under or related to this Agreement, including RECIPIENT’S use of the
    SOFTWARE. RECIPIENT agrees to indemnify, defend, and hold harmless PRINCETON (including its
    trustees, officers, faculty, employees, students and agents), from and against any and all claims,
    liabilities, or damages based on, arising out of, or relating to this Agreement, including RECIPIENT’s
    use of the SOFTWARE or RECIPIENT’S breach of this Agreement.
    8. All disputes regarding the construction, interpretation and the parties’ obligations under this Agreement
    shall be governed by the laws of the State of New Jersey, notwithstanding any of that state’s laws to the
    contrary. The venue and jurisdiction for the resolution of any such disputes shall be Mercer County,
    New Jersey
```

