<p align="center">
  <img height="200" src="images/logo.png">
</p>



Welcome to the Sei framework repository! Sei is a framework for systematically predicting sequence regulatory activities and applying sequence information to human genetics data. Sei provides a global map from any sequence to regulatory activities, as represented by 40 sequence classes, and each sequence class integrates predictions for 21,907 chromatin profiles (transcription factor, histone marks, and chromatin accessibility profiles across a wide range of cell types).

Sei is now published, you can read the manuscript [here](https://doi.org/10.1038/s41588-022-01102-2).

This repository can be used to run the Sei model and get the Sei chromatin profile and sequence class predictions for input sequences or variants.

We also provide information and instructions for [how to train the Sei deep learning sequence model](#training). 

### Requirements

Sei requires Python 3.6+ and Python packages PyTorch (>=1.0), Selene (>=0.5.0), and `docopt`. You can follow PyTorch installation steps [here](https://pytorch.org/get-started/locally/) and Selene installation steps [here](https://github.com/FunctionLab/selene). Install `docopt` with pip or conda (e.g. `conda install docopt`)

### Setup

Please download and extract the trained Sei model and `resources` (containing hg19 and hg38 FASTA files) `.tar.gz` files before proceeding:

```
sh ./download_data.sh
```

- [Sei model](https://doi.org/10.5281/zenodo.4906996)
- [Sei framework `resources` directory](https://doi.org/10.5281/zenodo.4906961)


### Chromatin profile prediction

1. The following scripts can be used to obtain Sei deep learning predictions for 21,907 chromatin profiles (please run on a GPU node):
(1) `1_sequence_prediction.py` (and corresponding bash script, `1_sequence_prediction.sh`): Accepts either a BED (`.bed`) or FASTA (`.fa`, `.fasta`) file as input and makes sequence predictions.

Example usage:
```
sh 1_sequence_prediction.sh <input-file> <genome> <output-dir> --cuda
```

Arguments:
- `<input-file>`: BED or FASTA input file
- `<genome>`: If you use a BED file as input, this must be either `hg19` or `hg38` as these are the FASTA reference genome files we provide by default. If you are using a FASTA file, you can specify whichever genome version you are using for logging purposes. 
- `<output-dir>`: Path to output directory (will be created if does not exist)
- `--cuda`: Optional, use this flag if running on a CUDA-enabled GPU.

You can run `python 1_sequence_prediction.py -h` for the full documentation of inputs.

2. `1_variant_effect_prediction.py` (and corresponding bash script, `1_variant_effect_prediction.sh`): Accepts a VCF file as input and makes variant effect predictions.

Example usage:
```
sh 1_variant_effect_prediction.sh <vcf> <hg> <output-dir> [--cuda]
```

Arguments:
- `<vcf>`: VCF file
- `<hg>`: Either hg19 or hg38
- `<output-dir>`: Path to output directory (will be created if does not exist)
- `--cuda`: Optional, use this flag if running on a CUDA-enabled GPU.

You can run `python 1_variant_effect_prediction.py -h` for the full documentation of inputs.

These scripts will output the chromatin profile predictions as HDF5 files to a subdirectory `chromatin-profiles-hdf5` in your specified output directory. 

See `example_slurm_scripts/1_example_seqpred.slurm_gpu.sh` and `example_slurm_scripts/1_example_vep.slurm_gpu.sh` for sample scripts for running chromatin profile prediction on SLURM.


### Sequence class prediction

Sequence class scores can be obtained from Sei chromatin profile predictions. There are 2 types of scores that can be computed:

- Raw sequence class scores: For sequences only. Raw sequence class scores are projection scores of chromatin profile predictions projected on the unit-length vectors representing each sequence class. This is an intermediate score originally developed for variant score prediction and is made available for use for developing downstream analyses or applications, such as using them as a sequence representation. **Note** our manuscript uses the Louvain community clustering, whole-genome sequence class annotation of the human genome whenever we apply sequence classes to reference genome sequences, and we encourage the use of these annotations over the raw sequence class scores when possible. Sequence class annotations for hg38 and hg19 (lifted over from hg38) are available for download from [this Zenodo record](10.5281/zenodo.7113989). 
- Sequence class variant effect score (nucleosome-occupancy-adjusted): For variants only. Computed as alt - ref of the raw sequence class scores **adjusted for nucleosome occupancy, i.e. histone normalized**. To better represent predicted variant effects on histone marks, it is necessary to normalize for nucleosome occupancy (for example, a LoF mutation near the TSS can decrease H3K4me3 modification level while increasing nucleosome occupancy, resulting in an overall increase in observed H3K4me3 quantity). Therefore, for variant effect computation, we used the sum of all histone profile predictions as an approximation to nucleosome occupancy and adjusted all histone mark predictions to remove the impact of nucleosome occupancy change (nonhistone mark predictions are unchanged). See manuscript methods for more detail.

#### Sequence prediction

Example usage:
```
sh 2_raw_sc_score.sh <input-file> <output-dir>
```

Arguments:
- `<input-file>`: Path to the Sei `_predictions.h5` file.
- `<output-dir>`: Path to output directory (will be created if does not exist)

You can run `python 2_raw_sc_score.py -h` for the full documentation of inputs.

#### Variant effect prediction

Example usage:
```
sh 2_varianteffect_sc_score.sh <ref-fp> <alt-fp> <output-dir> [--no-tsv]
```

Arguments:
- `<ref-fp>`: Path to the Sei `.ref_predictions.h5` file.
- `<alt-fp>`: Path to the Sei `.alt_predictions.h5` file.
- `<output-dir>`: Path to output directory (will be created if does not exist)
- `--no-tsv`: Optional flag if you'd like to suppress the outputted TSV files (see the next section 'Example variant effect prediction run' for more information).

You can run `python 2_varianteffect_sc_score.py -h` for the full documentation of inputs.

### Example variant effect prediction run:

We provide `test.vcf` (hg19 coordinates) so you can try running this command once you have installed all the requirements. Additionally, `example_slurm_scripts` contains example scripts with the same expected input arguments if you need to submit your job to a compute cluster. 

Example command run on GPU:
```
sh 1_variant_effect_prediction.sh test.vcf hg19 ./test_outputs --cuda
```

Example command run on CPU:
```
sh 2_varianteffect_sc_score.sh ./test_outputs/chromatin-profiles-hdf5/test.ref_predictions.h5 \
                               ./test_outputs/chromatin-profiles-hdf5/test.alt_predictions.h5 \
                               ./test_outputs
```
You can add `--no-tsv` to this command to suppress the TSV file outputs if you are comfortable working with HDF5 and NPY files. Note you will need to match the rows to the `test_row_labels.txt` file in `./test_outputs/chromatin-profiles.hdf5` and the columns to `./model/target.names` (chromatin profile HDF5 files) and `./model/seqclass.names` (sequence class NPY file). 


Expected outputs:
-  `chromatin-profiles-hdf5`: directory containing HDF5 Sei predictions files and the corresponding `test_row_labels.txt` file. 
- `sorted.test.chromatin_profile_diffs.tsv`: chromatin profile prediction TSV file (**Note:** output file will be compressed if input has >10000 variants), sorted by max absolute sequence class score. 
- `sorted.test.sequence_class_scores.tsv`: sequence class prediction TSV file, sorted by max absolute sequence class scores.
- `test.sequence_class_scores.npy`: sequence class scores NPY file, note this is NOT sorted and will be ordered in the same way as `chromatin-profiles-hdf5/test_row_labels.txt` file.

## Sequence classes

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

The training data is available [here](https://doi.org/10.5281/zenodo.4907037) should be downloaded and extracted into the `train` directory. 

**NOTE**: because the Sei training data contains processed files from the Cistrome Project, please first agree to the Cistrome Project [terms of usage](http://cistrome.org/db/#/bdown) before downloading the data:

```
cd ./train
sh ./download_data.sh  # in the train directory
```

The Sei training configuration YAML file is provided as the `train/train.yml` file. You can read more about the Selene command-line interface and configuration file formatting [here](https://selene.flatironinstitute.org/master/overview/cli.html#). 

You must use Selene version >0.5.0 to train this model ([release notes](https://github.com/FunctionLab/selene/blob/master/RELEASE_NOTES.md)). 

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
