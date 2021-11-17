setup
```bash
outdir=tests/data/GRCz11_chr9
cd $outdir
```

subset genome & gene annotation to chromosome 9
```bash
genomepy install GRCz11 -a -p ensembl -g . -r '>9'
```

subset bam to (some of) chromosome 9
```bash
bam=~/GRCz11-GSM1859511.samtools-coordinate.bam
samtools view -h $bam 9 -s 0.1 > chr9_reads.sam
samtools view -H chr9_reads.sam | grep @HD > chr9.sam
samtools view -H chr9_reads.sam | grep SN:9 >> chr9.sam
samtools view chr9_reads.sam >> chr9.sam
samtools view -b chr9.sam > chr9.bam
samtools index chr9.bam
rm chr9_reads.sam
rm chr9.sam
```

subset regions to chromosome 9
```bash
cat $regionsfile | grep -G ^9 > GRCz11_chr9_regions.bed
```

subset pfmfile & factor2motifs.txt 
(to pou factors as they are found on chr9 and have relatively few motifs)
```python
pfm = "/home/siebrenf/GRCz11.gimme.vertebrate.v5.0.pfm"
chr9_pfm = "tests/data/GRCz11_chr9/GRCz11_chr9.pfm"

# list of TFs on chromosome 9
chr9_tfs = ['hoxd3a', 'neurod1', 'creb1b', 'hoxd4a', 'ccnt2a', 'tfdp1a', 'sp5a', 'olig1', 'nr1i2', 'zic2a', 'pknox1.1', 'shox', 'nab1b', 'en1a', 'pou2f1b', 'sox1a', 'pou1f1', 'nfe2l2a', 'stat4', 'hoxd11a', 'pou3f3a', 'dlx1a', 'ikzf2', 'sp9', 'olig2', 'dlx2a', 'atf2', 'tbx15', 'nhlh2', 'hoxd9a', 'twist2', 'pou4f3', 'gli2a', 'fev', 'nr4a2a', 'etv5a', 'sp3a', 'itgb2', 'klf5b', 'tbr1b', 'zeb2a', 'zic5', 'evx2', 'klf12b', 'znf148', 'stat1b', 'gabpa', 'hoxd10a', 'ube2al', 'tfcp2l1', 'sox21b', 'hoxd12a', 'klf7b']

# pou subset
pou_tfs = ['pou2f1b', 'pou1f1', 'pou3f3a', 'pou4f3']
pou_motifs = ['GM.5.0.Homeodomain_POU.0009', 'GM.5.0.Homeodomain_POU.0008', 'GM.5.0.Homeodomain_POU.0002', 'GM.5.0.Homeodomain_POU.0025', 'GM.5.0.Homeodomain_POU.0012', 'GM.5.0.Homeodomain_POU.0003', 'GM.5.0.Mixed.0032', 'GM.5.0.Mixed.0066', 'GM.5.0.Mixed.0064', 'GM.5.0.Mixed.0101', 'GM.5.0.Homeodomain_POU.0015', 'GM.5.0.Homeodomain_POU.0023', 'GM.5.0.Mixed.0057', 'GM.5.0.Mixed.0054', 'GM.5.0.Homeodomain_POU.0001', 'GM.5.0.Homeodomain_POU.0014', 'GM.5.0.Homeodomain_POU.0026', 'GM.5.0.Mixed.0036', 'GM.5.0.Homeodomain_POU.0024', 'GM.5.0.Homeodomain_POU.0004', 'GM.5.0.Homeodomain_POU.0016', 'GM.5.0.Mixed.0062', 'GM.5.0.Homeodomain_POU.0006', 'GM.5.0.Mixed.0095', 'GM.5.0.Homeodomain_POU.0011', 'GM.5.0.Homeodomain_POU.0021', 'GM.5.0.Homeodomain_POU.0010', 'GM.5.0.Homeodomain_POU.0018', 'GM.5.0.Homeodomain_POU.0005']

keep = []
rec = False
with open(pfm) as f:
    for line in f:
        if line.startswith(">"):
            line = line.strip()
            if line[1:] in pou_motifs:
                rec = True
                keep.append(line + "\n")
            else:
                rec = False
        elif rec:
            keep.append(line)

with open(chr9_pfm, "w") as f:
    for line in keep:
        f.write(line)

```

subset pfmscorefile to chromosome 9 & pou motifs
```python
pfmscorefile = "/home/siebrenf/GRCz11.scan.tsv"
chr9_pfmscorefile = "tests/data/GRCz11_chr9/GRCz11_chr9.pfm"

import pandas as pd

pou_motifs = ['GM.5.0.Homeodomain_POU.0009', 'GM.5.0.Homeodomain_POU.0008', 'GM.5.0.Homeodomain_POU.0002', 'GM.5.0.Homeodomain_POU.0025', 'GM.5.0.Homeodomain_POU.0012', 'GM.5.0.Homeodomain_POU.0003', 'GM.5.0.Mixed.0032', 'GM.5.0.Mixed.0066', 'GM.5.0.Mixed.0064', 'GM.5.0.Mixed.0101', 'GM.5.0.Homeodomain_POU.0015', 'GM.5.0.Homeodomain_POU.0023', 'GM.5.0.Mixed.0057', 'GM.5.0.Mixed.0054', 'GM.5.0.Homeodomain_POU.0001', 'GM.5.0.Homeodomain_POU.0014', 'GM.5.0.Homeodomain_POU.0026', 'GM.5.0.Mixed.0036', 'GM.5.0.Homeodomain_POU.0024', 'GM.5.0.Homeodomain_POU.0004', 'GM.5.0.Homeodomain_POU.0016', 'GM.5.0.Mixed.0062', 'GM.5.0.Homeodomain_POU.0006', 'GM.5.0.Mixed.0095', 'GM.5.0.Homeodomain_POU.0011', 'GM.5.0.Homeodomain_POU.0021', 'GM.5.0.Homeodomain_POU.0010', 'GM.5.0.Homeodomain_POU.0018', 'GM.5.0.Homeodomain_POU.0005']
df = pd.read_csv(pfmscorefile, sep="\t", comment="#", index_col=0)
df = df[pou_motifs]
df = df[df.index.str.startswith("9:")]
df.to_csv(chr9_pfmscorefile, sep="\t")

```

subset expression file to chromosome 9
```python
qfile = "~/quant.sf"
chr9_qfile = "tests/data/GRCz11_chr9/expression_chr9_quant.sf"

import pandas as pd
import genomepy


# expression file contains transcript_ids
df = pd.read_csv(qfile, sep="\t")

a = genomepy.Annotation("tests/data/GRCz11_chr9/GRCz11")
ids = a.gtf.attribute.str.split('transcript_id "', expand=True)
ids = ids[1].dropna()
chr9_tids = list(set(ids.str.split('";', expand=True)[0]))

df = df.set_index("Name")
df2 = df[df.index.isin(chr9_tids)]
df2.to_csv(chr9_qfile, sep="\t")

```

create differential expression file for chromosome 9
```python
chr9_qfile = "tests/data/GRCz11_chr9/expression_chr9_quant.sf"
chr9_defile = "tests/data/GRCz11_chr9/chr9_diffexp.tsv"

import pandas as pd
import random


df = pd.read_csv(chr9_qfile, sep="\t", index_col=0)
df.index.name = None

df["log2FoldChange"] = [10*random.random() for _ in df.index]
df["padj"] = [random.random() for _ in df.index]

df = df[["log2FoldChange", "padj"]]
df.sort_values("padj", inplace=True)
df.to_csv(chr9_defile, sep="\t")

```

testrun on a subset of pou2f1b pou1f1 pou3f3a pou4f3
```bash
ananse binding \
-A chr9.bam \
-r GRCz11_chr9_regions.bed \
-f pou2f1b pou1f1 pou3f3a \
-g GRCz11/GRCz11.fa \
-p GRCz11_chr9.pfm \
--pfmscorefile GRCz11_chr9_scan.bed \
-n 1 \
-o GRCz11_binding
```

```bash
ananse network \
-b GRCz11_binding/binding.h5 \
-e expression_chr9_quant.sf \
-g GRCz11/GRCz11.fa \
--full-output \
-n 1 \
-o GRCz11_network
```

```bash
ananse influence \
-t GRCz11_network \
-d chr9_diffexp.tsv \
-a GRCz11/GRCz11.fa \
--full-output \
-n 1 \
-o GRCz11_influence
```
