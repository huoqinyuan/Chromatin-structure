import sys
from Bio import SeqIO,Seq,SeqRecord,SeqFeature
import numpy as np
import matplotlib.pyplot as plt
from single_seq_analyzer import single_seq_analyzer

growth_list = []
stress_list = []

for line in open("GS_list.txt"):
    line = line.split()
    if line[0] == "TFIID-dominated":
        growth_list.append(line[1])
    elif line[0] == "SAGA-dominated":
        stress_list.append(line[1])



chrms = list(SeqIO.parse("GCF_000146045.2_R64_genomic.gbff","genbank"))    #Seq_record of each chromosome
seqs = list(SeqIO.parse("GCF_000146045.2_R64_genomic.fna","fasta"))      #Seq
name = chrms[0].description.split()[0]+' '+chrms[0].description.split()[1]
length = 8000
bodylen = 8000
dinu_density_growth = []
dinu_density_stress = []
window = 40
num = (length+bodylen)//window
dinucleotide = "CA"
types = ['source', 'rRNA', 'telomere', 'ncRNA', 'mRNA', 'centromere', 'CDS', 'repeat_region', 'tRNA', 'misc_RNA', 'rep_origin', 'misc_feature', 'mobile_element']

target_seq_list_growth = []
target_seq_list_stress = []

for i in range(len(chrms)):
    gene_features = list(chrms[i].features)
    chrmlen = len(seqs[i].seq)
    SeqIO.write(chrms[i],"genes_tmp.gbff","genbank")
    genes_tmp = list(open("genes_tmp.gbff"))
    count = 0
    record_list_growth = []
    record_list_stress = []
    for j in range(50,len(genes_tmp)):
        if (genes_tmp[j].split()[0] in types) and (".." in genes_tmp[j].split()[1]):
            count += 1
        elif (genes_tmp[j].split()[0] == "gene") and (".." in genes_tmp[j].split()[1]):
            count += 1
            for k in range(1,2):
                if "/locus_tag" in genes_tmp[j+k]:
                    gen_name = genes_tmp[j+k].strip().strip("/locus_tag=").strip('"')
                    if gen_name in growth_list:
                        record_list_growth.append(count)
                    elif gen_name in stress_list:
                        record_list_stress.append(count)

    for index in record_list_growth:
        start_pos = int(gene_features[index].location.start)
        end_pos = int(gene_features[index].location.end)
        if start_pos > length and (start_pos + bodylen) < (chrmlen - 100):
            target_seq = seqs[i].seq[start_pos - length:start_pos + bodylen]
            if gene_features[index].location.strand == -1:
                target_seq.reverse_complement()
            target_seq.upper()
            target_seq_list_growth.append(target_seq)
    for index in record_list_stress:
        start_pos = int(gene_features[index].location.start)
        end_pos = int(gene_features[index].location.end)
        if start_pos > length and (start_pos + bodylen) < (chrmlen - 100):
            target_seq = seqs[i].seq[start_pos - length:start_pos + bodylen]
            if gene_features[index].location.strand == -1:
                target_seq.reverse_complement()
            target_seq.upper()
            target_seq_list_stress.append(target_seq)


for seq in target_seq_list_growth:
    one_density_growth = single_seq_analyzer(dinucleotide, seq, window)
    dinu_density_growth.append(one_density_growth)

for seq in target_seq_list_stress:
    one_density_stress = single_seq_analyzer(dinucleotide, seq, window)
    dinu_density_stress.append(one_density_stress)

dinu_avg_growth = np.zeros(num)
dinu_avg_stress = np.zeros(num)

for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_growth]
    tmp = np.array(tmp)
    dinu_avg_growth[i] = tmp.mean()

for i in range(num):
    tmp = [dinu[i] for dinu in dinu_density_stress]
    tmp = np.array(tmp)
    dinu_avg_stress[i] = tmp.mean()

Xline = np.linspace(-(length//window),(bodylen//window),num)
fig, ax = plt.subplots()
ax.plot(Xline,dinu_avg_growth,label = 'growth genes')
ax.plot(Xline,dinu_avg_stress,label = "stress genes")
ax.set_title("Saccharomyces cerevisiae growth&stress genes %s density" % (dinucleotide))
ax.set_xlabel("distance/40bp")
ax.set_ylabel(dinucleotide+" density")
ax.legend()
ax.grid()
fig.tight_layout()
plt.savefig("%s_density_growth&stress_%s" % (dinucleotide,window))
plt.clf()

