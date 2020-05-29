import sys
from Bio import SeqIO,Seq,SeqRecord,SeqFeature
import numpy as np
import matplotlib.pyplot as plt
from single_seq_analyzer import single_seq_analyzer
fln_gbff = sys.argv[1]   #gbff
fln_fna = sys.argv[2]    #fna

chrms = list(SeqIO.parse(fln_gbff,"genbank"))    #Seq_record of each chromosome
seqs = list(SeqIO.parse(fln_fna,"fasta"))      #Seq
name = chrms[0].description.split()[0]+' '+chrms[0].description.split()[1]
length = 8000
bodylen = 8000
dinu_density = []
window = 40
num = (length+bodylen)//window
dinucleotide = "AA"

for i in range(len(chrms)):
    gene_features = list(chrms[i].features)
    chrmlen = len(seqs[i].seq)
    gene_record_list = []
    target_seq_list = []
    for j in range(0,len(gene_features)):
        if gene_features[j].type == "gene":
            gene_record_list.append(gene_features[j])
            start_pos = int(gene_features[j].location.start)
            end_pos = int(gene_features[j].location.end)
            if start_pos > length and (start_pos+bodylen) < (chrmlen-100):
                target_seq = seqs[i].seq[start_pos - length:start_pos + bodylen]
                if gene_features[j].location.strand == -1:
                    target_seq.reverse_complement()
                target_seq.upper()
                target_seq_list.append(target_seq)


    for seq in target_seq_list:
        one_density = single_seq_analyzer(dinucleotide,seq,window)
        dinu_density.append(one_density)


dinu_density = np.array(dinu_density)

dinu_avg = np.zeros(num)
for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density]
    tmp = np.array(tmp)
    dinu_avg[i] = tmp.mean()

Xline = np.linspace(-(length//window),(bodylen//window),num)
fig, ax = plt.subplots()
ax.plot(Xline,dinu_avg)
ax.set_title(name)
ax.set_xlabel("distance/40bp")
ax.set_ylabel("%s density" % (dinucleotide))
fig.tight_layout()
plt.savefig("%s_%s_%s" % (dinucleotide,name,window))
plt.clf()



