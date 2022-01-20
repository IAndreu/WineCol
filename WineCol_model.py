from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn import svm, datasets
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from matplotlib.text import Annotation
import numpy as np
import pandas as pd
from Bio import SeqIO
import os
import WineCol_imputation
from imblearn.over_sampling import SMOTE


# polymorphic positions of the gene
positions = [129, 144, 147, 160, 168, 170, 190, 224, 231, 244, 258, 358, 360, 372, 376, 380, 394, 395, 418, 425, 447, 453, 468, 485, 490, 518, 551, 557, 562, 563, 570, 573, 575, 577, 580, 594, 603, 622, 631, 649, 678, 702, 715, 724, 760, 761, 771, 776, 805, 821, 822, 825, 826, 868, 872, 920, 1014, 1024, 1051, 1073, 1084, 1087, 1111, 1133, 1145, 1154, 1155, 1158, 1160, 1171, 1175, 1180, 1213, 1219, 1221, 1228, 1230, 1256, 1257, 1263, 1275, 1293, 1326, 1340, 1341, 1343, 1345, 1348, 1349, 1352, 1535, 1537, 1551, 1552, 1553, 1560, 1563, 1574, 1585, 1595, 1596, 1597, 1610, 1636, 1637, 1650, 1659, 1667, 1689, 1691, 1698, 1718, 1719, 1727, 1729, 1731, 1733, 1739, 1755, 1763, 1773, 1812, 1827, 1838, 1882, 1896, 1902, 1907, 1911, 1917, 1926, 1940, 1941, 1970, 1972, 1973, 1974]

# positions with SNPs
imputed_snp_positions = [129, 144, 147, 160, 168, 190, 231, 244, 358, 360, 372, 376, 380, 394, 395, 418, 425, 447, 453, 468, 485, 490, 518, 551, 557, 562, 563, 570, 573, 575, 577, 580, 594, 603, 622, 631, 649, 678, 702, 715, 724, 760, 761, 821, 822, 825, 868, 872, 920, 1014, 1024, 1051, 1073, 1084, 1087, 1111, 1133, 1145, 1154, 1155, 1158, 1160, 1171, 1175, 1180, 1213, 1219, 1221, 1228, 1230, 1256, 1257, 1263, 1275, 1293, 1326, 1340, 1341, 1343, 1345, 1348, 1349, 1535, 1537, 1551, 1552, 1553, 1560, 1563, 1574, 1585, 1595, 1596, 1597, 1610, 1636, 1637, 1650, 1659, 1667, 1689, 1691, 1718, 1719, 1727, 1729, 1731, 1733, 1739, 1755, 1763, 1773, 1827, 1882, 1896, 1902, 1907, 1917, 1940, 1941]

def read_seqs(MSA):
    fasta_sequences = SeqIO.parse(open(MSA),'fasta')
    seqs={'B':[], 'N':[], 'NR':[], 'Rg':[], 'Rs':[]}
    for fasta in fasta_sequences:
        for color in seqs.keys():
            if fasta.id.endswith("_"+color):
                seqs[color].append((fasta.id, str(fasta.seq)))
    return seqs
                
def CountFrequency(my_list):
    # Creating an empty dictionary
    freq = {}
    length = len(my_list)
    for item in my_list:
        if (item in freq):
            freq[item] += 1
        else:
            freq[item] = 1
    for key, value in freq.items():
        freq[key]/=length    
    return(freq)

def get_frequencies(seqs):
    msa_len=len(seqs['B'][0][1])
    positions={i:[] for i in range(msa_len)}
    for color in seqs:
        for sequence in seqs[color]:
            for i in range(msa_len):
                positions[i].append(sequence[1][i])

    frequencies={i:[] for i in range(msa_len)}
    for i in positions:
        frequencies[i]=CountFrequency(positions[i])
    return frequencies

def impute2_to_seq(file, freq):
    imputed_file = [i.split() for i in open(file,'r').readlines()]
    imputed_file.reverse()
    seq = []
    strand={'A':'T','C':'G','T':'A','G':'C','N':'N','-':'-'}
    for idx ,i in enumerate(imputed_snp_positions):
        dominant_allele= strand[list(freq[i])[0]]
        if imputed_file[idx][5]==imputed_file[idx][6]:
            if imputed_file[idx][3+int(imputed_file[idx][5])]==dominant_allele:
                seq.append(0)
            else:
                seq.append(2)
        else:
            seq.append(1)
    return seq 
    
def encoder0(freq, pos):
    freq=list(freq)
    for poly in freq:
        if poly not in ['L','-','N']:
            higher_freq=poly
            break
    if pos=='-' or pos=='L': # if is a gap -10
        return -1
    elif pos=='N': # if is an N
        return 0
    elif pos==higher_freq:
        return 0
    else:
        if pos in ['A','C','T','G']:
            return 2
        else:
            return 1
                
# Encode the training sequences and obtain the unique representations per color
def encode_train(seqs, frequencies, encoder, gret1):
    colors_code={'B':1, 'N':2, 'NR':3, 'Rg':4, 'Rs':5}
    X=[]
    y=[]
    labels=[]
    for color in seqs:
        for sequence in seqs[color]:
            seq=[]
            for i, nucl in enumerate(str(sequence[1])):
                if i in positions:
                    seq.append(encoder(frequencies[i],nucl))
            X.append(seq)
            y.append(colors_code[color])
            labels.append(sequence[0])
    if gret1!=3: # use gret information
        for i in range(len(y)):
            if int(y[i])==1:
                X[i].append(100)
            else:
                X[i].append(0)
    colors={'white':set(), 'black':set(), 'darkred':set(),  'red':set(), 'pink':set()}
    col_sym={'_B':'white', '_N':'black', 'NR':'darkred',  'Rg':'red', 'Rs':'pink'}
    for i in range(len(labels)):
        for j in range(len(labels)):
            if list(X[i])==list(X[j]):
                    colors[col_sym[labels[i][-2:]]].add(labels[j])
    id_seq={}
    for i in range(len(labels)):
        id_seq[labels[i]]=X[i]

    new_color={}
    for col in colors:
        new_color[col]=[]
        for seq in colors[col]:
            new_color[col].append(id_seq[seq])

    white = list(np.unique(np.array(new_color['white']),axis=0))
    black = list(np.unique(np.array(new_color['black']),axis=0))
    darkred = list(np.unique(np.array(new_color['darkred']),axis=0))
    red = list(np.unique(np.array(new_color['red']),axis=0))
    pink = list(np.unique(np.array(new_color['pink']),axis=0))
    X = np.array(white+black+darkred+red+pink)
    y = np.array([1 for i in range(len(white))]+[2 for i in range(len(black))]+[3 for i in range(len(darkred))]+[4 for i in range(len(red))]+[5 for i in range(len(pink))])
    
    return X, y

def encode_ancient(sequence, frequencies, encoder, gret1, imputed_snps):
    seq=[]
    for idx, pos in enumerate(sequence):
        seq.append(encoder(frequencies[idx], pos))
    
    if imputed_snps:
        for idx, snp in enumerate(imputed_snp_positions):
            seq[snp]=imputed_snps[idx]

    seq= [seq[i] for i in positions]        
    if int(gret1)==0:
        seq.append(0)
    if int(gret1)==1:
        seq.append(100)

    seq = np.array(seq)
    
    return seq



def predict_grape(train, ancient, gret1, impute):
    id_to_color={1:'White', 2:'Black', 3:'Black reddish', 4:'Red', 5:'Rose'}  
    seqs_train = read_seqs(train)
    frequencies = get_frequencies(seqs_train)
    X_train,y_train = encode_train(seqs_train, frequencies, encoder0, int(gret1))
    oversample = SMOTE()
    X_train,y_train = oversample.fit_resample(X_train,y_train)
    ancient_seq = str(list(SeqIO.parse(open(ancient),'fasta'))[0].seq)
    if impute!='No':
        WineCol_imputation.impute_ancient(impute, ancient)
        imputed_snps = impute2_to_seq("output/"+ancient.split('.')[0]+"_haps", frequencies)
    else:
        imputed_snps=[]
    
    X_test = encode_ancient(ancient_seq,frequencies, encoder0, gret1, imputed_snps) 
    
    rbf = svm.SVC(kernel='rbf', gamma=1, C=10, decision_function_shape='ovo', class_weight='balanced').fit(X_train, y_train)
    rbf_pred =rbf.predict(X_test.reshape(1, -1))
    
    print('The estimated color of your berry is... %s \n' % ( id_to_color[int(rbf_pred)]))
