## Command-line reference

In general, a full analysis with ANANSE will consist of the following steps:

1. Generate binding networks for a *source* and a *target* cell type using `ananse binding`.
2. Generate a gene regulatory network (GRN) for the *source* and *target* cell type using `ananse network`, based on the binding network from step 1.
3. Run `ananse influence` based on the GRN of the *source* cell type and the GRN of the *target* cell type (step 2)

### ananse binding

The `ananse binding` command predicts binding for all transcription factors with an associated motif.

Example command:

``` bash
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding
```

This command will use `hg38` by default and quantify the ATAC and H3K27ac signal in [putative enhancers](https://doi.org/10.5281/zenodo.4066424). This signal will be used together with the motif scores to predict transcription factor binding. If you use `hg38` and the default enhancer regions you can get significantly better predictions by using a more detailed model that also includes the average ChIP-seq signal. You will first have to download this model (you only need to do this once):

```
DATA_DIR=/path/to/data  # set to your preference
mkdir -p $DATA_DIR/ANANSE.REMAP.model.v1.0
cd $DATA_DIR/ANANSE.REMAP.model.v1.0
wget https://zenodo.org/record/4768075/files/ANANSE.REMAP.model.v1.0.tgz
tar xvzf ANANSE.REMAP.model.v1.0.tgz
rm ANANSE.REMAP.model.v1.0.tgz
```

You can now specify this model when running `ananse binding`:

``` bash
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding \
               -R $DATA_DIR/ANANSE.REMAP.model.v1.0
```

Note that the BAM file(s) should be indexed, for instance with `samtools index`.

Alternatively, you can use your own set of enhancer regions or putative cis-regulatory elements. For instance, to use the candidate cis-Regulatory Elements by ENCODE (from [SCREEN](https://screen.encodeproject.org/)):

``` bash
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding \
               -r https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed
```

Or to use the union of peaks from several ATAC-seq experiments:

``` bash
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding \
               -r ATAC.cell_type1.narrowPeak ATAC.cell_type2.narrowPeak
```

#### Other species

Please note that human and mouse are supported out-of-the-box, but for other species you will need a custom motif to transcription factor mapping. Take care to use the same gene symbols for your transcription factors as are present in your gene expression file. One way to create the TF to motif mapping is to use the `motif2factors` command from [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs). See the documentation [here](https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors). This command uses a relatively simple approach based on orthogroups. 


#### Output files

The output directory of `ananse binding` will contain one or two files, depending on the input data used. The file called `binding.h5` contains:

* The predicted binding probability for every TF in every enhancer. For every transcription factor there is a key with the transcription factor name and a value that contains the predicted binding probability.
* The predicted factor activity (based on the TF motif activity) for all transcriptionfactors. The motif activity is based on the contribution of motif score to the enhancer activity (the coefficient from a linear regression;  see [The FANTOM Consortium & Riken Omics Science Center 2009](https://doi.org/10.1038/ng.375) and [Balwierz et al. 2014](https://doi.org/10.1101/gr.169508.113)). In the implementation of ANANSE the motif activity is calculated based on a regression using the ATAC-seq signal and/or a regression using the H3K27ac signal.
* Depending on the input data, the file will also contain `atac` and/or `h3k27ac`, which provide the quantified and normalized signal of ATAC-seq and/or H3K27ac ChIP-seq data in the enhancer regions. If the `ANANSE.REMAP.model.v1.0` model is used, these files will also contain the relative signal. 

Mostly, you would want to use the `binding.h5` as input for `ananse network`. If you want to access the information in the file for other purposes you can
use the `ananse view` command.

If you provided multiple peaks or BED files as input, the output directory will also contain a BED files with the merged regions, which are used as potential enhancer.

#### Full options

Usage:

```
ananse [-h] <command> [options] binding [-A BAM [BAM ...]] [-H BAM [BAM ...]] [-o] [-R PATH] [-r  [...]] [-g GENOME]
                                               [-d] [-p] [-f [TF [TF ...]]] [-n NCPUS] [-h]
```

Required arguments:

```
  -A BAM [BAM ...], --atac-bams BAM [BAM ...]
                        ATAC-seq input BAM file(s), can be used alone or in combination with the -H option
  -H BAM [BAM ...], --histone-bams BAM [BAM ...]
                        H3K27ac ChIP-seq input BAM file(s), can be used alone or in combibation with the -A option
```

Optional arguments:

```
  -o , --outdir         Directory where you wish to store the output (default: ./ANANSE_binding)
  -R PATH, --reference PATH
                        Path to reference data directory
  -r  [ ...], --regionfiles  [ ...]
                        One or more BED format files with putative enhancer regions (e.g. BED, narrowPeak, broadPeak)
  -g GENOME, --genome GENOME
                        Genome name, as managed by genomepy, or path to the genome FASTA used to align the bams and peaks
                        to
  -d , --dist-func      bam reads are normalized to the selected distribution (default: an empirical distribution)
  -p , --pfmfile        PFM file of the transcription factors to search for (default gimme.vertebrate.v5.0)
  -f [TF [TF ...]], --factors [TF [TF ...]]
                        Transcription factors to use. Either a space-separated list or a file with one TF per line.
  -n NCPUS, --ncpus NCPUS
                        Number of processes to use for motif scanning
  --pfmscorefile        use precomputed gimmemotifs scores (gimme scan -T -g
                        GENOME INPUTFILE)
  -h, --help            Show this help message and exit
```


### ananse network

The `ananse network` command infers a cell type-specific GRN based on the predicted TF binding sites using `ananse binding` and the expression levels of both TFs as well as their target genes. TF-gene interaction scores, the edge weights in the network, are calculated based on the predicted TF binding probability in enhancers and the distance between the enhancers and the target gene, the predicted TF activity and the expression level of both TF and the target gene.

Note: `ananse network` needs ~12-15GB of memory for a typical analysis of a human network.

#### Full options


Usage: 

```
ananse [-h] <command> [options] network -b FILE -e FILE -o FILE [-g NAME] [-a BED] [-n NCORE] [--include-promoter]
                                        [--include-enhancer] [-h]
```

Required arguments:

```
  -b FILE, --binding FILE
                        TF binding prediction file (binding.h5 from ananse binding).
  -e FILE, --expression FILE
                        Expression scores. Should have gene names in the first column and should contain a column named
                        tpm. Both the quant.sf from salmon or the abundances.tsv from kallisto will work fine.
  -o FILE, --outfile	Name of output file.
```

Optional arguments:

```
  -g NAME, --genome NAME
                        Genome (genomepy name or FASTA file).
  -a BED, --annotation BED
                        Gene annotation in BED12 format. Not necessary if you use a genome that was installed using
                        genomepy with the --annotation flag.
  -n NCORE, --ncore NCORE
                        Number of core used.
  --include-promoter, --exclude-promoter
                        Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference. By default promoter peaks
                        are included.
  --include-enhancer, --exclude-enhancer
                        Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference. By default enhancer peaks
                        are included.
  -h, --help            show this help message and exit
```


### ananse influence

To calculate the influence score for the transition from a *source* cell type (`-s` or `--source`) to a *target* cell type (`t` or `--target`), `ananse influence` uses the GRNs for both cell types, predicted by `ananse network`. For each network, the top 100k interactions are selected, based on the rank of the interaction scores (edge weights). Using the differential GRN, capturing the difference between the two networks, a local network is built for each TF, up to a maximal number of three edges. Using this network, the influence score is calculated based on 1) the edge distance from the TF of interest to the target gene, 2) the predicted interaction score and 3) the change in expression between the source cell type and the target cell type.

The number of edges used is 100,000 by default but in some cases you'll get better performance with more edges. You can, for instance, try to increase to 500,000 edges by using `-i 500_000`. 

Example command: 

``` bash
$ ananse influence  -s results/FB_network.txt \
                    -t results/KRT_network.txt \
                    -d data/FB2KRT_degenes.tsv \
                    -o results/influence.txt \
                    -n 8
```

#### Full options

Usage: 

```
ananse [-h] <command> [options] influence -t FILE -d FILE -o FILE [-s FILE] [-i EDGES] [-p] [-n NCORE] [-h]
```

Required arguments:

```
  -t FILE, --target FILE
                        Network of target cell type.
  -d FILE, --degenes FILE
                        File with differential gene expression (DEseq2 output file).
  -o FILE, --outfile FILE
                        Output file.     
```                        

Optional arguments:

```
  -s FILE, --source FILE
                        Network of source cell type.
  -i EDGES, --edges EDGES
                        Number of top edges used (default is 100,000).
  -p, --plot            Create influence score plot.
  -n NCORE, --ncore NCORE
                        Number of cores to use.
  -h, --help            show this help message and exit
```


### ananse view

Convert the binding probabilities from  the `binding.h5` output of `ananse binding` to tab-separated text format.
This is only necessary if you want to use the output for other purposes, as `ananse network` can only use the `.h5` file.
Converting all factors may take some time and memory!

Example command to extract the binding probabilities for all factors in *wide* format (one column per TF):

``` bash
$ ananse view binding.h5 -o binding.tsv
```

Example command to extract the binding probabilities for TP53 and TP63 in *long* format, print to stdout:


``` bash
$ ananse view binding.h5 -f TP53 TP63 -F long
```

```
loc     factor  prob
chr1:181357-181557      TP53    0.2432
chr1:267938-268138      TP53    0.2568
chr1:586086-586286      TP53    0.2452
chr1:605299-605499      TP53    0.2445
chr1:629832-630032      TP53    0.994
chr1:631289-631489      TP53    0.965
chr1:631431-631631      TP53    0.689
chr1:633919-634119      TP53    0.9966
chr1:778533-778733      TP53    0.777
[...many more rows...]
```

#### Full options

```
usage: ananse [-h] <command> [options] view [-o FILE] [-f [TF [TF ...]]] [-F FORMAT] FILE

positional arguments:
  FILE                  input binding.h5 file

optional arguments:
  -o FILE, --outfile FILE
                        outputfile (tab-separated text, default: stdout)
  -f [TF [TF ...]], --factors [TF [TF ...]]
                        name(s) of transcription factors (default: all)
  -F FORMAT, --format FORMAT
                        format: wide or long (default: wide)
```