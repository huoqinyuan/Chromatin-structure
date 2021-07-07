from Bio import Seq
import numpy as np

def single_seq_analyzer(dinucleotide,seq,window):
    length =len(seq)
    num = length//window
    one_density = np.zeros(num)
    for i in range(num):
        seq_tmp = seq[i*window:(i+1)*window]
        Ncount = seq_tmp.count("N")
        count_tmp = seq_tmp.count(dinucleotide)
        one_density[i] = count_tmp/(float(window)-Ncount+0.00000000001)
    return one_density


