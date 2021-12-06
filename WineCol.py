import argparse
import WineCol_model

parser = argparse.ArgumentParser(description='Estimates color of grapes from genotype data using a support vector machine')
parser.add_argument('train', type = str,
                     help = 'Fasta file consisting of a multiple sequence alignment of the myAb1 genes used to train the model')
parser.add_argument('test', type = str,
                     help = 'Fasta file with the sequence/s of the myAb1 gene whose color will be estimated')

parser.add_argument("-e", type=int, choices=[0, 1, 2],
                    help="encoder to use. Default 0")
args = parser.parse_args()
train = args.train
test = args.test
if args.e:
    encoder = str(args.e)
else:
    encoder = 0

WineCol_model.predict_grape(train, test, encoder)
