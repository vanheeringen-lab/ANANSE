from gimmemotifs.motif import read_motifs

xenbase_id_file = "/data/xenopus/xenbase/XenbaseGenepageToGeneIdMapping.txt"
anno = "XENTR_9.1_Xenbase.xt9_0.bed" 
with open(xenbase_id_file) as f:
    gene_names = [line.split("\t")[1] for line in f]

with open(anno) as f:
    gene_names += [line.split("\t")[3] for line in f]

with open("gimme.vertebrate.v3.2.xt.pwm") as f:
    motifs = read_motifs(f)
    for motif in motifs:
        motif.factors = [f.lower() for f in motif.factors if 
                f.lower() in gene_names]

with open("gimme.vertebrate.v3.3.xt.pwm", "w") as fout:
    for motif in motifs:
        fout.write("{}\n".format(motif.to_pwm()))

with open("gimme.vertebrate.v3.3.xt.motif2factors.txt", "w") as fout:
    for motif in motifs:
        fout.write("{}\t{}\n".format(
            motif.id,
            ",".join(motif.factors)
            ))

