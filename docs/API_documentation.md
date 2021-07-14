## API document

### Prediction of transcription factor binding

This uses the example data described [here](examples.html).

```
from ananse.peakpredictor import predict_peaks
import glob

outdir = "heart.binding"
histone_files = glob.glob("ANANSE_example_data/atac/heart*bam")
atac_files = glob.glob("ANANSE_example_data/H3K27ac/heart*bam")
reference = "ANANSE.REMAP.model.v1.0/"

def binding(args):
    predict_peaks(
        outdir,
        atac_bams=atac_files,
        histone_bams=histone_files,
        reference=reference,
        genome="hg38",
    )
```

### Working with `Network` Class
The `Network` class include all functions used to infer cell type specific gene regulatory network.

```python
from ananse.network import Network
from dask.distributed import Client

b = ananse.network.Network(
    genome="hg38", 
    include_promoter=True,
    include_enhancer=True
)

# Simple default dask

with Client() as client:
    b.run_network(
        binding="heart.binding/binding.h5",
        fin_expression=glob.glob("heart.binding/heart*TPM.txt"),
        outfile="heart.network.txt",

```


### Working with `influence` Class

