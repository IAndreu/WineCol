import argparse
import WineCol_model

parser = argparse.ArgumentParser(description='Estimates color of grapes from genotype data using a support vector machine')
parser.add_argument('train', type = str,
                     help = 'Fasta file consisting of a multiple sequence alignment of the myAb1 genes used to train the model')
parser.add_argument('test', type = str,
                     help = 'Fasta file with the sequence/s of the myAb1 gene whose color will be estimated')

parser.add_argument("-snp", type=int, choices=[0, 1],
                    help="Use only the polymorphisms of the input data. Default 1 (0:No, 1:Yes)")
parser.add_argument("-e", type=int, choices=[0, 1, 2],
                    help="encoder to use. Default 0")
parser.add_argument("-gret1", type=int, choices=[0, 1],
                    help="Take into account Gret1 retrotransposon presence. Default 0 (0:No, 1:Yes)")
parser.add_argument("-impute2", type=int, choices=[0, 1],
                    help="Perform imputation of the input data (not implemented yet). Default 0 (0:No, 1:Yes)")

args = parser.parse_args()
train = args.train
test = args.test
if args.snp:
    snp = str(args.snp)
else:
    snp = 1
if args.gret1:
    gret1 = str(args.gret)
else:
    gret1 = 0
if args.e:
    encoder = str(args.e)
else:
    encoder = 0
if args.impute2:
    impute2 = str(args.impute2)
else:
    impute2 = 0
    
    
WineCol_model.predict_grape(train, test, encoder)
