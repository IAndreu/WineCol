import argparse
import WineCol_model

parser = argparse.ArgumentParser(description='Estimates color of grapes from genotype data using a support vector machine')
parser.add_argument('ancient', type = str,
                     help = 'Fasta file with the aligned sequence of the myAb1 gene whose color will be predicted.')

parser.add_argument('-train', type = str,
                     help = "Fasta file consisting of a multiple sequence alignment of the myAb1 genes used to train the model. Default is the one in the 'data/' folder.")

parser.add_argument("-gret1", type=int, choices=[0, 1],
                    help="Presence of Gret1 retrotransposon (0:No, 1:Yes).")
parser.add_argument("-impute", type=str,
                    help="Perform imputation of the input data. PATH to BAM file of the ancient sequence.")

args = parser.parse_args()
ancient = args.ancient


if args.train:
    train = str(args.train)
else:
    train = "Data/MSA_grapes_ordered.fa"
if args.gret1==0:
    gret1 = 0
elif args.gret1==1:
    gret1=1
else:
    gret1 = 3
if args.impute:
    impute = str(args.impute)
else:
    impute = "No"
    
    
WineCol_model.predict_grape(train, ancient, gret1, impute)
