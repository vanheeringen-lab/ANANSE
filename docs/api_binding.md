## Working with `Binding` class
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
