import sys
from Bio import SeqIO,Seq,SeqRecord,SeqFeature
import numpy as np
import matplotlib.pyplot as plt
from single_seq_analyzer import single_seq_analyzer
import gzip

fln_gbff = sys.argv[1]   #gbff
fln_fna = sys.argv[2]    #fna

chrms = list(SeqIO.parse(fln_gbff,"genbank"))    #Seq_record of each chromosome
seqs = list(SeqIO.parse(fln_fna,"fasta"))      #Seq
name = chrms[0].description.split()[0]+' '+chrms[0].description.split()[1]
length = 8000
bodylen = 8000
dinu_density_AA = []
dinu_density_CG = []
dinu_density_AT = []
dinu_density_CC = []
window = 40
num = (length+bodylen)//window

print("a")

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
        one_density_AA = single_seq_analyzer("AA",seq,window)
        dinu_density_AA.append(one_density_AA)
        one_density_CG = single_seq_analyzer("CG",seq,window)
        dinu_density_CG.append(one_density_CG)
        one_density_AT = single_seq_analyzer("AT", seq, window)
        dinu_density_AT.append(one_density_AT)
        one_density_CC = single_seq_analyzer("CC", seq, window)
        dinu_density_CC.append(one_density_CC)


dinu_density_AA = np.array(dinu_density_AA)
dinu_density_CG = np.array(dinu_density_CG)
dinu_density_AT = np.array(dinu_density_AT)
dinu_density_CC = np.array(dinu_density_CC)

print("b")

dinu_avg_AA = np.zeros(num)
for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_AA]
    tmp = np.array(tmp)
    dinu_avg_AA[i] = tmp.mean()

dinu_avg_CG = np.zeros(num)
for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_CG]
    tmp = np.array(tmp)
    dinu_avg_CG[i] = tmp.mean()

dinu_avg_AT = np.zeros(num)
for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_AT]
    tmp = np.array(tmp)
    dinu_avg_AT[i] = tmp.mean()

dinu_avg_CC = np.zeros(num)
for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_CC]
    tmp = np.array(tmp)
    dinu_avg_CC[i] = tmp.mean()

Xline = np.linspace(-(length//window),(bodylen//window),num)
fig, ax = plt.subplots()
ax.plot(Xline,dinu_avg_AA,label = "ApA density")
ax.plot(Xline,dinu_avg_CG, label = "CpG density")
ax.plot(Xline,dinu_avg_AT,label = "ApT density")
ax.plot(Xline,dinu_avg_CC, label = "CpC density")
ax.legend(loc = 7,bbox_to_anchor=(0, -0.005))
ax.set_title(name)
ax.set_xlabel("distance/40bp")
ax.set_ylabel("dinucleotide density")
ax.grid()
fig.tight_layout()
plt.savefig("4dinucleotide_density_%s_%s" % (name,window))
plt.clf()














