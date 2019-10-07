## Working with binding
The `Binding` class include all functions used to infer TF-binding site for all enhancer peaks.

```
import sys 
import importlib
sys.path.append('../')

import pandas as pd
import dask.dataframe as dd

import grns.binding
```

```
gene_bed = "/home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed"
peak_bed = "data/krt_enhancer.bed"
pwmfile = "../data/gimme.vertebrate.v5.1.pfm"


a=grns.binding.Binding(genome="hg38", 
                        gene_bed= gene_bed, pwmfile=pwmfile)


filter_bed = a.clear_peaks(peak_bed)

pwm_weight = a.get_PWMScore(filter_bed.name)
pwm = dd.read_csv(pwm_weight.name, sep="\t")

peak_weight = a.get_peakRPKM(filter_bed.name)
peak = dd.read_csv(peak_weight.name, sep="\t")

table=a.get_binding_score(pwm, peak)

```
