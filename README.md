# WineCol

Method to estimate the colour of ancient grape berries from genotype data using a support vector machine.


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


Run example:

```
python3 WineCol.py Data/MDV14_US13525_P7.fa
``` 
