# ANANSE: ANalysis Algorithm for Networks Specified by Enhancers

### Prediction of key transcription factors in cell fate determination using enhancer networks

![](/pic/Fig2.jpg)
> (A), Illustration of all the data used to predict key TFs in cell conversion. Those data include the enhancer database from ATAC-seq, DNase-seq or p300 ChIP-seq, the motif score of all TFs and the gene expression data of each cell type from RNA-seq. (B), The predicted cell-type specific TF binding profiles from enhancer database and TF's motif score in each cell type. (C), The predicted cell-type specific GRN based on TF/Gene binding, TF/Gene expression and its' distance. (D), The GRN difference between two interested cell types. (E), The ranked influence score of all TFs calculated from GRN.


## Quick start

### Installation

The most straightforward way to install ANANSE is by using [bioconda](https://bioconda.github.io/).

If you have not used bioconda before, first install [conda](https://docs.continuum.io/anaconda/) and then set up the necessary channels (in this order!). You only have to do this once.

```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

Now you can install ANANSE:

```
# Install all dependencies
$ conda create -n ananse python=3 gimmemotifs networkx chest dask pytables
adjusttext

# Activate the environment
$ conda activate ananse 

# Upgrade gimmemotifs to development version
$ pip install git+https://github.com/vanheeringen-lab/gimmemotifs.git@develop
```

For most of the analyses it is beneficial to use as many threads as possible for the motif analysis. This is configured by the GimmeMotifs config file. If you haven't done so, run `gimme`, which will create a new GimmeMotifs config.

```
$ gimme
```

Now edit the file `~/.config/gimmemotifs/gimmemotifs.cfg`, and change the `ncpus` parameter.

### Data files.

You need to download the genome of interest.

```
$ genomepy install hg38 UCSC --annotation
```

### Easy installation

* Activate the environment

```
$ conda activate ananse 
```

* Install `ananse` development version package from github
```bash
$ pip install git+https://github.com/vanheeringen-lab/ANANSE.git@develop
```

Python 3 is the required. Don't forget to activate the environment with `conda activate ananse` whenever you want to use `grns`.

### API documentation

* The ***python API documentation*** of this package can be found [***here***](/docs/api.md).


### Build binding network

In this command, the `-r/--fin_rpkm` and `-o/--outfile` is the required arguments. The `-r` input is a enhancer peak bed file. This is the example of this file:
```bash
chr2	148881617	148881817	7
chr7	145997204	145997404	4
chr13	109424160	109424360	20
chr14	32484901	32485101	2
chr3	93470526	93470726	315
chrY	56728139	56728339	50
chr15	82513939	82514139	36
chr2	49229595	49229795	3
chr12	41363656	41363856	2
chr8	69689994	69690194	13
```
The first column is chromosome name, the second and third column is the start and end point of peak. We recommend all peaks have 200bp. The fourth column is intensity of the peak, it could be RPKM or equivalent value.

* Example:
```
$ ananse binding  -r data/krt_enhancer.bed \
                    -o results/binding.txt \
                    -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                    -g hg38 \
                    -p ../data/gimme.vertebrate.v5.1.pfm
```

* All the optional arguments:
```
$ ananse binding -h

usage: ananse [-h] <subcommand> [options] binding [-h] [-p FILE] -r FILE
                                                    [-a BED] [-g GENOME] -o
                                                    FILE [-f NAME] [-d NAME]
optional arguments:
  -h, --help            show this help message and exit
  -p FILE, --pwmfile FILE
                        PWM
  -r FILE, --fin_rpkm FILE
                        BED file with RPKM (or equivalent) value on the 4th column
  -a BED, --annotation BED
                        Gene annotation in BED12 format
  -g GENOME, --genome GENOME
                        Genome
  -o FILE, --outfile FILE
                        Output file
  -f NAME, --filter_promoter NAME
                        Filter promoters, True or False, input should be
                        either 'True' or 'False'. (Default setting: True; if 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.)
  -d NAME, --keep_detail NAME
                        Keep detail files, True or False, input should be
                        either 'True' or 'False'. (Default setting: True)
```

### Built interaction network

* Example:
```
$ ananse interaction  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                        -r data/krt_enhancer.bed \
                        -o results/full_features.txt \
                        -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                        -g hg38 \
                        -b results/binding.txt \
                        -c /home/qxu/projects/regulatoryNetwork/history/cell_trans/human_gene_correlation/expressioncorrelation.txt \
                        -p ../data/gimme.vertebrate.v5.1.pfm
```
* Input
```
-e One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM; 
-o The folder to save results;
-a 12 columns BED file with gene annotation;
-g Genome;
-b The binding network from binding.py;
-c All gene correlation file;
-p motifs file (optional; if provided there should also be a motif2factors.txt).
```

### Built GRN

* Example:
```
$ ananse network -f results/full_features.txt -o results/full_network.txt
```
* Input
```
-f The interaction network from interaction.py;
-o The folder to save results.
```

### Infer influence score

* Example:
```
$ ananse influence    -a results/full_network.txt \
                        -e data/FB_rep1_TPM.txt \
                        -d data/FB2KRT_degenes.csv \
                        -o results/FB2KRT.txt
                        -p False
```
* Input
```
-b The network in first cell (optional);
-a The network in second cell;
-e The gene expression in first cell (optional);
-d The differential expression table between two cells; 
-o The result file.
-p Plot influence.

```

## Help

* The preferred way to get support is through the Github issues page

