from gimmemotifs.fasta import Fasta
f = Fasta("XENTR_9.1_Xenbase.fa")

m = {}
with open("rna2name.txt") as fin:
    for line in fin:
        rna,name = line.strip().split("\t")
        m[rna] = name

for name,seq in f.items():
    print(">{}::{}".format(m[name.split(":")[0]], name.split(":")[0]))
    print(seq)
