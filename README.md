# XXX Network

Prediction of key transcription factors in cell fate determination using enhancer networks

# Quick start

## Dependencies

### Install Anaconda and activate bioconda channels.

```
# Install all dependencies
$ conda create -n regnetwork python=3 gimmemotifs genomepy networkx chest dask pytables

# Activate the environment
$ conda activate regnetwork

# install dependency via pip
$ pip install adjustText
```

### Run `gimme` to create a new GimmeMotifs config.

```
$ gimme
```

### Edit the file `~/.config/gimmemotifs/gimmemotifs.cfg`, and change the `ncpus` parameter.

### Download the genome of interest.

```
$ genomepy install Xenopus_tropicalis_v9.1 NCBI
```

## Built binding network

### Example:
```
$ python ../grns/binding.py -r data/krt_enhancer.bed \
                            -o results \
                            -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                            -g hg38 \
                            -p ../data/gimme.vertebrate.v5.1.pfm
```
### input

```
-r The enhancer BED file with enhancers (200bp) with RPKM (or equivalent) value in 4th column;
-o The folder to save results;
-a 12 columns BED file with gene annotation;
-g Genome;
-p Motifs file (optional; if provided there should also be a motif2factors.txt).
-f Filter promoters, input should be either 'True' or 'False'. (Default setting: True; if 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.)
-d Keep temporary files, input should be either 'True' or 'False'. (Default setting: True)
```

## Built interaction network

### Example:
```
$ python ../grns/interaction.py -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                                -o results \
                                -a /home/qxu/.local/share/genomes/hg38/hg38_gffbed_piroteinCoding.bed \
                                -g hg38 \
                                -b results/binding.predicted.h5 \
                                -c /home/qxu/projects/regulatoryNetwork/history/cell_trans/human_gene_correlation/expressioncorrelation.txt \
                                -p ../data/gimme.vertebrate.v5.1.pfm
```
### input
```
-e One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM; 
-o The folder to save results;
-a 12 columns BED file with gene annotation;
-g Genome;
-b The binding network from binding.py;
-c All gene correlation file;
-p motifs file (optional; if provided there should also be a motif2factors.txt).
```

## Built GRN

### Example:
```
$ python ../grns/builtnetwork.py -f results/full_features.h5 -o results
```
### input
```
-f The interaction network from interaction.py;
-o The folder to save results.
```

# Help

