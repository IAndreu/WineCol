# WineCol

Method to estimate the colour of ancient grape berries from genotype data using a support vector machine. 

This tool takes as input an ancient gene VvmybA1, a transcriptional regulator of anthocyanin biosynthesis (principal pigment of grapes). Having trained a model with more than 400 sequences of modern VvmybA1 genes, the tool returns the estimated phenotypic colour of the ancient sample. In order to deal with the low coverage of the ancient sequences, the tool offers the option of imputing the input sequence by inferring the genotypes using an haplotype reference panel. 

#### Contents
1. Prerequisites
2. Installation
3. Usage    
4. Example


### 1. Prerequisites

The requiered dependencies necessary for running the tool without imputation are:
- **Python**: Download the newest version available from https://www.python.org/downloads/
- **Biopython**: Download from https://biopython.org/wiki/Download
- **Scikit-learn**: Download from: https://scikit-learn.org/stable/install.html

The requiered dependencies necessary for performing the imputation are:
- **Java 8**
- **GATK3.6**: Download from https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2
- **IMPUTE2**: Download from https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#download

### 2. Installation

Download this git repository.

### 3. Usage

```
usage: WineCol.py [-h] [-train TRAIN] [-gret1 {0,1}] [-impute IMPUTE] ancient

Estimates color of grapes from genotype data using a support vector machine

positional arguments:
  ancient         Fasta file with the aligned sequence of the myAb1 gene whose color will be predicted.

optional arguments:
  -h, --help      show this help message and exit
  -train TRAIN    Fasta file consisting of a multiple sequence alignment of the myAb1 genes used to train the model. Default is the one in the 'data/' folder.
  -gret1 {0,1}    Presence of Gret1 retrotransposon (0:No, 1:Yes).
  -impute IMPUTE  Perform imputation of the input data. PATH to BAM file of the ancient sequence.
  ```
If imputation is requiered, download all the files starting with "grape12Xv2" from https://sid.erda.dk/sharelink/hSM6HJzoha and place the files in the "Data/" folder. Also, specify the paths for GenomeAnalysisTK and impute2 in the file WineCol_imputation.py.

Run example:

```
python3 WineCol.py Data/MDV14_US13525_P7.fa
``` 
