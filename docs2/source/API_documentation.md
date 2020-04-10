## API document

### Working with `Binding` Class
The `Binding` class include all functions used to infer TF-binding site for all enhancer peaks.

```python
import pandas as pd
import dask.dataframe as dd

from ananse import binding
```

```python
gene_bed = "/data/hg38_genes.bed"
peak_bed = "data/krt_enhancer.bed"
pfmfile = "/data/gimme.vertebrate.v5.1.pfm"
```

A new `Binding` object `bid`
```python
bid=binding.Binding(genome="hg38", gene_bed= gene_bed, pfmfile=pfmfile)
```

Using `clear_peak()` function, we can filter the peaks in promoter ranges. 
```python
filter_bed = bid.clear_peak(peak_bed)
```

The `get_PWMScore()` function, calculate the motif z-score of all TFs in all peaks.
```python
pfm_weight = bid.get_PWMScore(filter_bed)
pfm = dd.read_csv(pfm_weight, sep="\t")
```

Load enhancer peak intensity.
```python
peak_weight = bid.get_peakRPKM(filter_bed)
peak = dd.read_csv(peak_weight, sep="\t")
```

Infer the TF-binding score.
```python
table=bid.get_binding_score(pfm, peak)
```

The `peak_bed` file is the enhancer peak file, and `"./"` is the output dir. With `run_binding()` function, we can infer the binding sites of all TFs in `peak_bed` enhancer ranges.

```python
bid.run_binding(peak_bed,"./")
```

### Working with `Network` Class
The `Network` class include all functions used to infer cell type specific gene regulatory network.

```python
import pandas as pd
import dask.dataframe as dd

from ananse import network
```

```python
fin_expression = ["data/KRT_rep1_TPM.txt", "data/KRT_rep2_TPM.txt"]
binding = "results/binding.txt"
gene_bed = "/data/hg38_genes.bed"
peak_bed = "data/krt_enhancer.bed"
pfmfile = "/data/gimme.vertebrate.v5.1.pfm"
corrfiles = "expressioncorrelation.txt"
```

A new `Network` object `net`
```python
net=network.Network(genome="hg38", gene_bed= gene_bed, pfmfile=pfmfile)
```

Using `clear_peak()` function, we can filter the peaks in promoter ranges. 
```python
filter_bed = net.clear_peak(peak_bed)
```

Using `get_promoter_dataframe()` function to overlap the enhancer with promoter, and using `get_gene_dataframe()` function to overlap the enhancer with gene.
```python
prom = net.get_promoter_dataframe(filter_bed)
p = net.get_gene_dataframe(filter_bed)
```

Using `weight()` function to establish weight function 100 kb up and downstream to the TSS.
```python
weight = net.distance_weight()
```

Using `aggregate_binding()` function to merge binding data with gene.
```python
features = net.aggregate_binding(binding, prom, p, weight)
```

Using `aggregate_binding()` function to establish and merge expression data with features, and using `get_correlation()` function to establsih and merge expression data with features.
```python
expression_file = net.get_expression(fin_expression, features)
corr_file = net.get_correlation(corrfiles, features)

other = [
    expression_file,
    corr_file,
]

featurefile = net.join_features(features, other)
```

Using `create_network()` function to establish gene regulatory network.
```python
net.create_network(featurefile, outfile)
```

### Working with `influence` Class

