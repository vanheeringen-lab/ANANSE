## Working with `Binding` class
The `Binding` class include all functions used to infer TF-binding site for all enhancer peaks.

```python
import sys 
import importlib
sys.path.append('../')

import pandas as pd
import dask.dataframe as dd

from grns import binding
```

```python
gene_bed = "/home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed"
peak_bed = "data/krt_enhancer.bed"
pwmfile = "../data/gimme.vertebrate.v5.1.pfm"
```

A new `Binding` object `a`
```python
a=binding.Binding(genome="hg38", gene_bed= gene_bed, pwmfile=pwmfile)
```

Using `clear_peaks()` function, we can filter the peaks in promoter ranges. 
```python
filter_bed = a.clear_peaks(peak_bed)
```

The `get_PWMScore()` function, calculate the motif z-score of all TFs in all peaks.
```python
pwm_weight = a.get_PWMScore(filter_bed)
pwm = dd.read_csv(pwm_weight, sep="\t")
```

Load enhancer peak intensity.
```python
peak_weight = a.get_peakRPKM(filter_bed)
peak = dd.read_csv(peak_weight, sep="\t")
```

Infer the TF-binding score.
```python
table=a.get_binding_score(pwm, peak)
```

The `peak_bed` file is the enhancer peak file, and `"./"` is the output dir. With `run_binding()` function, we can infer the binding sites of all TFs in `peak_bed` enhancer ranges.

```python
xxx=a.run_binding(peak_bed,"./")
```
