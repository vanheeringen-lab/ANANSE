# Examples

## Prerequisites

1. Please install ANANSE if you haven't done so already, using the [installation documentation](installation.md).

2. If you installed ANANSE via conda, which is recommended, make sure to activate the environment before you run it.

```
conda activate ananse
```

3. Install the `hg38` genome and annotation via genomepy.

```
genomepy install hg38 --annotation
```

4. Download the ANANSE model that incorporates REMAP data.

```
mkdir -p ANANSE.REMAP.model.v1.0
wget https://zenodo.org/record/4768075/files/ANANSE.REMAP.model.v1.0.tgz
tar xvzf ANANSE.REMAP.model.v1.0.tgz -C ANANSE.REMAP.model.v1.0
rm ANANSE.REMAP.model.v1.0.tgz
```

5. Download and unpack the example data:

```
wget https://zenodo.org/record/4769814/files/ANANSE_example_data.tgz
tar xvzf ANANSE_example_data.tgz
rm ANANSE_example_data.tgz
```


##  Data

To run a full ANANSE analysis you will need:

* ATAC-seq and/or H3K27ac ChIP-seq BAM files.
* RNA-seq expression quantification file (which contains TPM).
* RNA-seq differential expresson file (DEseq2).

These files are present in the example data for fibroblasts and for primary heart tissue:

```
$ tree ANANSE_example_data/

ANANSE_example_data/
├── ATAC
│   ├── fibroblast_ATAC_rep1.bam
│   ├── fibroblast_ATAC_rep1.bam.bai
│   ├── fibroblast_ATAC_rep2.bam
│   ├── fibroblast_ATAC_rep2.bam.bai
│   ├── heart_ATAC_rep1.bam
│   ├── heart_ATAC_rep1.bam.bai
│   ├── heart_ATAC_rep2.bam
│   └── heart_ATAC_rep2.bam.bai
├── H3K27ac
│   ├── fibroblast_H3K27ac_rep1.bam
│   ├── heart_H3K27ac_rep1.bam
│   └── heart_H3K27ac_rep1.bam.bai
├── README.txt
└── RNAseq
    ├── fibroblast2heart_degenes.csv
    ├── fibroblast_rep1_TPM.txt
    ├── fibroblast_rep2_TPM.txt
    ├── heart_rep1_TPM.txt
    └── heart_rep2_TPM.txt

3 directories, 17 files
```


Details for the different steps are described below.

### Prediction of TF binding

To predict TF binding, the main data that you'll need is ATAC-seq and/or H3K27ac ChIP-seq BAM files for all samples, conditions or cell types that you want to analyze. The TF binding prediction works best with both ATAC-seq and H3K27ac, but either of the two will still work.
You can use any number of replicates. The data of different replicates will be combined.
You can supply your own set of (putative) cis-regulatory regions (CREs), however, for this tutorial we will use a set of cis-regulatory regions that was created by combining all REMAP ChIP-seq peaks (see [here](ananse-cre.md)). Using this set of regions, the binding prediction will be more accurate.

Let's create predict some binding!

```
ananse binding -H ANANSE_example_data/H3K27ac/fibroblast*bam -A ANANSE_example_data/ATAC/fibroblast*bam -R ANANSE.REMAP.model.v1.0/ -o fibroblast.binding
ananse binding -H ANANSE_example_data/H3K27ac/heart*bam -A ANANSE_example_data/ATAC/heart*bam -R ANANSE.REMAP.model.v1.0/ -o heart.binding
```

### Gene regulatory network inference

Note: `ananse network` needs ~12-15GB of memory for a typical analysis of a human network.

```
ananse network -b  fibroblast.binding/binding.h5 -e ANANSE_example_data/RNAseq/fibroblast*TPM.txt -n 8 > fibroblast.network.txt
ananse network -b  heart.binding/binding.h5 -e ANANSE_example_data/RNAseq/heart*TPM.txt -n 8 > heart.network.txt
```

### Influence score calculation

```
ananse influence -s fibroblast.network.txt -t heart.network.txt -d ANANSE_example_data/RNAseq/fibroblast2heart_degenes.csv -p -o fibroblast2heart.influence.txt 
```



