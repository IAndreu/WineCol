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



def read_seqs(MSA):
    fasta_sequences = SeqIO.parse(open(MSA),'fasta')
    seqs={'B':[], 'G':[], 'N':[], 'NR':[], 'Nt':[], 'Rg':[], 'Rs':[], 'nd':[]}
    for fasta in fasta_sequences:
        c=0
        for color in seqs.keys():
            if fasta.id.endswith("_"+color):
                seqs[color].append((fasta.id, str(fasta.seq)))
                c=100
            c+=1
            if c==8:
                seqs['nd'].append((fasta.id, str(fasta.seq)))
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

        
def encoder0(freq, pos):
    code={"M":['A','C'],"Y":['C','T'],"R":['A','G'],"K":['G','T'],"W":['A','T'],"S":['C','G']}
    distance_code={1:10,2:17,3:20}
    freq=list(freq)
    for poly in freq:
        if poly not in ['L','-','N']:
            higher_freq=poly
            break
    if pos=='-' or pos=='L': # if is a gap -10
        return -10
    elif pos=='N': # if is an N 1
        return 1
    elif pos==higher_freq: # if is the nucleotide with higher frequency (implement in the future the restriction that the most frequent can be A/C)
        return 1
    else:
        distance=0
        for poly in freq:
            if poly==pos and pos not in ['A','C','T','G']:
                if higher_freq in code[poly]:
                    return 7
                else:
                    return 13
            if poly==pos and pos in ['A','C','T','G']:
                return distance_code[distance]
            if poly in ['A','C','T','G']:
                distance+=1
                
def encode_train(seqs, frequencies, encoder):
    colors={'B':1, 'G':2, 'N':3, 'NR':4, 'Nt':5, 'Rg':6, 'Rs':7, 'nd':8}        
    X=[]
    y=[]
    for color in seqs:
        for sequence in seqs[color]:
            seq=[]
            for i, nucl in enumerate(str(sequence[1])):
                if encoder==0:
                    seq.append(encoder0(frequencies[i],nucl))
            X.append(seq)
            y.append(colors[color])
    X = np.array(X)
    y = np.array(y)
    
    return X, y

def encode_test(test, frequencies, encoder):
    fasta_sequences = SeqIO.parse(open(test),'fasta')
    X=[]
    ids=[]
    for fasta in fasta_sequences:
        seq=[]
        for i,nucl in enumerate(str(fasta.seq)):
            if encoder==0:
                seq.append(encoder0(frequencies[i],nucl))
        X.append(seq)
        ids.append(fasta.id)
    return X, ids



def predict_grape(train, test, encoder):
    id_to_color={1:'White', 2:'Gray', 3:'Black', 4:'Black reddish', 5:'unknown', 6:'Red', 7:'Rose', 8:'unknown'}  
    seqs_train = read_seqs(train)
    frequencies = get_frequencies(seqs_train)
    X_train,y_train = encode_train(seqs_train, frequencies, encoder)
    X_test, ids = encode_test(test, frequencies, encoder)
    linear = svm.SVC(kernel='linear', C=0.1, decision_function_shape='ovo').fit(X_train, y_train)
    linear_pred = linear.predict(X_test)
    for i, j in enumerate(linear_pred):
        print('The estimated color of your grape %s is... %s \n' % (ids[i], id_to_color[j]))
