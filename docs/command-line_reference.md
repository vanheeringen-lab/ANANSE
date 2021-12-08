## Command-line reference

In general, a full analysis with ANANSE will consist of the following steps:

1. Generate binding networks for a *source*, and a *target* cell type using `ananse binding`.
2. Generate a gene regulatory network (GRN) for the *source* and *target* cell type using `ananse network`, based on the binding network from step 1.
3. Run `ananse influence` based on the GRN of the *source* cell type and the GRN of the *target* cell type (step 2)

### ananse binding

The `ananse binding` command predicts binding for all transcription factors with an associated motif.
It quantifies the ATAC and/or H3K27ac signal in [putative enhancers](https://doi.org/10.5281/zenodo.4066424). 
This signal will be used together with motif scores to predict transcription factor binding. 

Example command (using `hg38`):

```shell
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding
```

Example command (using `GRCz11`):

```shell
ananse binding -A ATAC.rep1.bam ATAC.rep2.bam \
               -H H3K27ac.rep1.bam H3K27ac.rep2.bam \
               -g GRCz11.fa \
               -r ATAC_macs2.narrowPeak \
               -p GRCz11.pfm \
               -o GRCz11_binding
```

Notes:

* The BAM file(s) should be indexed, for instance with `samtools index`.
* If you are using `hg38` or `hg19`, a genome fasta is required.
* For other genomes, a genome fasta with matching annotation is required.
    * [please see here](../input_data/#genome) for more info.
    * for other genomes, **additional steps are required** (see below)

#### Human (hg38) only: REMAP data

??? Expand
    <p></p>

    If you use `hg38`, you can get a set of default enhancer regions, and significantly better predictions, 
    by using a more detailed model (that also includes the average ChIP-seq signal).

    You will first have to download this model (you only need to do this once):
    ```shell
    DATA_DIR=/path/to/data  # set to your preference
    mkdir -p $DATA_DIR/ANANSE.REMAP.model.v1.0
    cd $DATA_DIR/ANANSE.REMAP.model.v1.0
    wget https://zenodo.org/record/4768075/files/ANANSE.REMAP.model.v1.0.tgz
    tar xvzf ANANSE.REMAP.model.v1.0.tgz
    rm ANANSE.REMAP.model.v1.0.tgz
    ```
    
    You can then specify this model when running `ananse binding`:
    ```shell
    ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
                   -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
                   -o IPS.binding \
                   -R $DATA_DIR/ANANSE.REMAP.model.v1.0
    ```

<p></p> [comment]: <- (empty line)

#### What, where?

Unless you use the `hg38` REMAP data, you will need 
a **genome** with matching **gene annotation**, 
a set of **regions** (where to look), 
and a **motif database** (what to look for) with motif to transcription factor mapping.
Take care to select a gene annotation with the same gene symbols for your transcription factors 
as are present in your (differential) gene expression files.

##### --regions
You can use your own set of enhancer regions or putative cis-regulatory elements. 
This will limit the analysis to ATAC/ChIP peaks in these regions.

For instance, to use the candidate cis-Regulatory Elements by ENCODE (from [SCREEN](https://screen.encodeproject.org/)):
```shell
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding \
               -r https://api.wenglab.org/screen_v13/fdownloads/GRCh38-ccREs.bed
```

Or to use the union of peaks from several ATAC-seq experiments:
```shell
ananse binding -A IPS.ATAC.rep1.bam IPS.ATAC.rep2.bam \
               -H IPS.H3K27ac.rep1.bam IPS.H3K27ac.rep2.bam \
               -o IPS.binding \
               -r ATAC.cell_type1.narrowPeak ATAC.cell_type2.narrowPeak
```

If you are using `seq2science` for your ATAC-/ChIP-seq analysis, you can use a counts table directly.

##### --genome
The filepath to the genome FASTA used to align your BAMs to (e.g. `~/genomes/GRCz11/GRCz11.fa`).

If you installed the genome with genomepy in a specific directory, the directory will work too (e.g. `~/genomes/GRCz11`).
Or, if you installed the genome with genomepy in the default directory, the genome name will do (`GRCz11`).

##### --pfmfile
The default motif database supports human and mouse genes out-of-the-box.

For other species you need a custom TF to motif mapping.
One way to create this mapping is to use the `gimme motif2factors` command from [GimmeMotifs](https://github.com/vanheeringen-lab/gimmemotifs). 
See the documentation [here](https://gimmemotifs.readthedocs.io/en/stable/reference.html#command-gimme-motif2factors). 
This command uses a relatively simple approach based on orthogroups.

If you specify a file, e.g. `mymotifs.pfm`, the command expects `mymotifs.motif2factors.txt` in the same folder.

<p></p> [comment]: <- (empty line)

#### Additional options

??? Expand
    <p></p>

    ##### --pfmscorefile
    If you expect to be running `ananse binding` multiple times on the same regions, you can precompute the motif scores once.
    This can be done with command `gimme scan -Tz --gc -g GENOME REGIONS > SCAN.tsv`.

    * GENOME is the genome fasta file you used to align the BAM files to.
    * REGIONS can be one FASTA, BED (including narrowPeak) or TSV with `chr:start-end` in the first column
    (similar to `--regions` above).
    * SCAN.tsv is your desired output filepath.
    
    ##### --jaccard-cutoff
    Not every transcription factor is present in the ANANSE database.
    For factors lacking a model, you can use another factor that shares one or more binding motifs.
    The jaccard-cutoff sets the minimum similarity score (between 0 for no similarity, 1 for perfect similarity).
    Use a higher value e.g. 0.2 to improve robustness of the selected model.

<p></p> [comment]: <- (empty line)

#### Output files

??? Expand

    The output directory of `ananse binding` will contain one or two files, depending on the input data used. 
    The file called `binding.h5` contains:

    * The predicted binding probability for every TF in every enhancer. 
    For every transcription factor there is a key with the transcription factor name and a value that contains the predicted binding probability.
    
    * The predicted factor activity (based on the TF motif activity) for all transcription factors. 
    The motif activity is based on the contribution of motif score to the enhancer activity (the coefficient from a linear regression; 
    see [The FANTOM Consortium & Riken Omics Science Center 2009](https://doi.org/10.1038/ng.375) and [Balwierz et al. 2014](https://doi.org/10.1101/gr.169508.113)). 
    In the implementation of ANANSE the motif activity is calculated based on a regression using the ATAC-seq signal and/or a regression using the H3K27ac signal.
    
    * Depending on the input data, the file will also contain `atac` and/or `h3k27ac`, which provide the quantified and normalized signal of ATAC-seq and/or H3K27ac ChIP-seq data in the enhancer regions. 
    If the `ANANSE.REMAP.model.v1.0` model is used, these files will also contain the relative signal. 
    
    Mostly, you would want to use the `binding.h5` as input for `ananse network`. 
    If you want to access the information in the file for other purposes you can use the `ananse view` command.
    
    If you provided multiple peaks or BED files as input, the output directory will also contain a BED files with the merged regions, which are used as potential enhancer.

<p></p> [comment]: <- (empty line)

#### Full options

```shell
usage: ananse [-h] <command> [options] binding [-A BAM [BAM ...]] [-H BAM [BAM ...]] [-g NAME] [-r FILE [FILE ...]] [-p FILE] [-o DIR] [-R DIR]
                                               [--pfmscorefile FILE] [-t [TF ...]] [--jaccard-cutoff FLOAT] [-n INT] [-h]

Required arguments:
  -A BAM [BAM ...], --atac-bams BAM [BAM ...]
                        ATAC-seq input BAM file(s), can be used alone or in combination with the -H option
  -H BAM [BAM ...], --histone-bams BAM [BAM ...]
                        H3K27ac ChIP-seq input BAM file(s), can be used alone or in combination with the -A option
  -g NAME, --genome NAME
                        Genome (genomepy name or FASTA file) used to align the BAMs and regions to (default: hg38)

Required arguments (optional for hg38):
  -r FILE [FILE ...], --regions FILE [FILE ...]
                        Regions to analyse. Can be one or more BED format files (e.g. BED, narrowPeak, broadPeak) or one file with one region per line
                        (e.g. 'chr1:100-200') or a space-separated list. Optional if a pfmscorefile is provided (used to filter those regions instead)
  -p FILE, --pfmfile FILE
                        PFM file of the transcription factors to search for (default: gimme.vertebrate.v5.0)

Optional arguments:
  -o DIR, --outdir DIR  Directory where you wish to store the output (default: ./ANANSE_binding)
  -R DIR, --reference DIR
                        Path to reference data directory
  --pfmscorefile FILE   Use precomputed gimmemotifs scores (gimme scan -Tz --gc -g GENOME REGIONS > SCAN.tsv)
  -t [TF ...], --tfs [TF ...]
                        Filter Transcription Factors to use (default: all in motif2factors.txt). Either a space-separated list or one or more files with
                        one TF per line
  --jaccard-cutoff FLOAT
                        TFs with a jaccard motif similarity >= the cutoff can be used as backup model. 0: any similarity, 1: perfect similarity (default
                        is 0.1)
  -n INT, --ncore INT   Number of cores to use.
  -h, --help            show this help message and exit
```

### ananse network

The `ananse network` command infers a cell type-specific GRN based on the predicted TF binding sites,
using `ananse binding` and the expression levels of both TFs and their target genes. 
TF-gene interaction scores and network edge weights are calculated based on four metrics:

1. the predicted TF binding probability in enhancers,
2. the distance between the enhancers and the target gene, 
3. the predicted TF activity, 
4. and the expression level of both TF and the target gene.

Note: `ananse network` needs ~12-15GB of memory for a typical analysis of a human network.

Example command:
```shell
$ ananse network -b ANANSE_binding/binding.h5 \
                 -e quant.sf \
                 -g GRCz11/GRCz11.fa \
                 -n 1 \
                 -o ANANSE_network.tsv
```

##### Gene expression

Gene expression in TPM (abundances). 
These can be extracted from one or more files, such as Salmon quant.sf files, or seq2science count-TPM tables.
The name(s) of the column(s) to extract can be specified (case-insensitive) with the `--column` argument.

##### Gene symbols

The gene symbols used in the expression file must match with those in the pfmfile.
For Human and mouse data, these are HGNC names. 
For a pfmfile generated with `gimme motif2factors`, these are usually gene names (or gene ids, depending on your gene annotation GTF file).

If you provide a genomepy name (with annotation) to the `--genome` and/or `--annotation` argument, 
the symbols in the expression file will be converted to gene names automatically (if the original names don't match).
If gene symbols are converted from transcripts/gene ids to gene names, expression levels are summed for duplicates.

#### Full options

```shell
usage: ananse [-h] <command> [options] network [-e FILE [FILE ...]] [-g NAME] [-a BED] [-o FILE] [-t [TF ...]] [-r FILE] [-c COL] [-f]
                                               [--include-promoter] [--include-enhancer] [-n INT] [-h]
                                               FILE

required arguments:
  FILE                  TF binding prediction file (ANANSE binding output)
  -e FILE [FILE ...], --expression FILE [FILE ...]
                        Gene expression file(s) with genes as first column and expression column(s) in TPM values (column name(s) specified below). Genes
                        must be in HGNC symbols, unless a genomepy gene annotation is provided. Files can include transcript-level TPMs (the quant.sf
                        from salmon or the abundances.tsv from kallisto), or gene-level TPMs tables (summarized with e.g. tximeta).

Required arguments (optional for hg38):
  -g NAME, --genome NAME
                        Genome (genomepy name or FASTA file) used to align the BAMs and regions to (default: hg38)
  -a BED, --annotation BED
                        Gene annotation (genomepy name or BED12 file) used to quantify expression levels. Optional when a genomepy genome (with
                        annotation) is providing.

optional arguments:
  -o FILE, --outfile FILE
                        Name of the output network file (default: ./ANANSE_network.tsv)
  -t [TF ...], --tfs [TF ...]
                        Filter Transcription Factors to use (default: all in motif2factors.txt). Either a space-separated list or one or more files with
                        one TF per line
  -r FILE, --regions FILE
                        Filter regions to use (default: all in binding.h5). Either one region/BED format file or a space-separated list.
  -c COL, --column COL  One or more (case insensitive) column names to extract from the expression file(s) (default: tpm)
  -f, --full-output     Export the full GRN output to the output file
  --include-promoter, --exclude-promoter
                        Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference. By default promoter peaks are included.
  --include-enhancer, --exclude-enhancer
                        Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference. By default enhancer peaks are included.
  -n INT, --ncore INT   Number of cores to use.
  -h, --help            show this help message and exit
```

### ananse influence

To calculate the influence score for the transition from a *source* cell type (`-s` or `--source`) to a *target* cell type (`t` or `--target`), 
`ananse influence` uses the GRNs for both cell types, predicted by `ananse network`. 
For each network, the top 100k interactions are selected, based on the rank of the interaction scores (edge weights). 
Using the differential GRN, capturing the difference between the two networks, a local network is built for each TF, up to a maximal number of three edges. 
Using this network, the influence score is calculated based on:

1. the edge distance from the TF of interest to the target gene, 
2. the predicted interaction score and 
3. the change in expression between the source cell type and the target cell type.

The number of edges used is 100,000 by default but in some cases you'll get better performance with more edges. 
You can, for instance, try to increase to 500,000 edges by using `-i 500_000`. 

Example command: 

```shell
$ ananse influence  -s results/FB_network.txt \
                    -t results/KRT_network.txt \
                    -d data/FB2KRT_degenes.tsv \
                    -o results/influence.txt \
                    -n 8
```

##### Gene symbols

The gene symbols used in the differential expression file must match with those in the pfmfile.
For Human and mouse data, these are HGNC names. 
For a pfmfile generated with `gimme motif2factors`, these are usually gene names (or gene ids, depending on your gene annotation GTF file).

If you provide a genomepy name (with annotation) to the `--annotation` argument, 
the symbols in the differential expression file will be converted to gene names automatically (if the original names don't match).
If gene symbols are converted from transcripts/gene ids to gene names, the lowest adjusted p-value is selected for duplicates.

#### Full options

```shell
usage: ananse [-h] <command> [options] influence -t FILE -d FILE [-s FILE] [-o FILE] [-f] [-a GTF] [-i INT] [-j FLOAT] [-c STR] [-n INT] [-h]

required arguments:
  -t FILE, --target FILE
                        Network of target cell type.
  -d FILE, --degenes FILE
                        File with differential gene expression (DEseq2 output file). Genes must be in HGNC symbols, unless a genomepy gene annotation is
                        provided.

optional arguments:
  -s FILE, --source FILE
                        Network of source cell type.
  -o FILE, --outfile FILE
                        Name of the output influence file (default: ./ANANSE_influence.tsv)
  -f, --full-output     Export the full GRN output to the output file
  -a GTF, --annotation GTF
                        Gene annotation (genomepy name or GTF file) used to quantify expression levels
  -i INT, --edges INT   Number of top edges used (default: 100.000).
  -j FLOAT, --padj FLOAT
                        Adjusted p-value below which genes classify as differential (default: 0.05).
  -c STR, --GRNsort-column STR
                        Column of GRN file sorted to select top interactions, 'prob' by default
  -n INT, --ncore INT   Number of cores to use.
  -h, --help            show this help message and exit
```

### ananse plot

Plot the output of `ananse influence` to a dotplot.
When providing an `influence_diffnetwork.tsv` also plot a GRN image of the top TFs their interactions.

Example command: 

```shell
$ ananse plot -i results/influence.tsv \
              --diff-network results/influence_diffnetwork.tsv \
              -o results
```

#### Full options

```shell
usage: ananse [-h] <command> [options] plot [-d FILE] [-o DIR] [--edge-info EDGE_INFO] [--edge-min EDGE_MIN] [--node-placement NETWORK_ALGORITHM]
                                            [--n-tfs N_TFS] [-c CMAP] [-f] [-t FTYPE] [-h]
                                            FILE

required arguments:
  FILE                  TF influence file (ANANSE influence output)
  -d FILE, --diff-network FILE
                        TF influence diffnetwork file (also ANANSE influence output)

optional arguments:
  -o DIR, --outdir DIR  Directory where you wish to store the output (default: ./ANANSE_plot)
  --edge-info EDGE_INFO
                        Column to use for edges of GRN, default: 'weight'. When full_output is specified, options are 'wb_diff' ,'tf_act_diff',
                        'tf_expr_diff', 'tg_expr_diff'
  --edge-min EDGE_MIN   Minimum value for an edge to be included in the GRN image
  --node-placement NETWORK_ALGORITHM
                        pyviz cluster algorithm used for node placement, options include: neato, dot, fdp, twopi, sfdp, circo
  --n-tfs N_TFS         Amount of TFs to plot in the GRN, default is top 20 differential TFs
  -c CMAP, --cmap CMAP  matlotlib colour library
  -f, --full-output     Select if the diffnetwork is a full output file
  -t FTYPE, --type FTYPE
                        Specify the output filetype (default: pdf)
  -h, --help            show this help message and exit
```

### ananse view

Convert the binding probabilities from  the `binding.h5` output of `ananse binding` to tab-separated text format.
This is only necessary if you want to visualize the output, or use it for other purposes, as `ananse network` can only use the `.h5` file.
Converting all factors may take some time and memory!

Visualize the first 10 regions and motifs:

```shell
$ ananse view binding.h5 -n 10
```

List the regions or transcription factors in the binding file:

```shell
$ ananse view binding.h5 -lr
$ ananse view binding.h5 -lt
```

Example command to extract the binding probabilities for all factors in *wide* format (one column per TF):

```shell
$ ananse view binding.h5 -o binding.tsv
```

Example command to extract the binding probabilities for TP53 and TP63 in *long* format, print to stdout:

```shell
$ ananse view binding.h5 -f TP53 TP63 -F long

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

```shell
usage: ananse [-h] <command> [options] view [-o FILE] [-t [TF ...]] [-r [REGION ...]] [-F FORMAT] [-n INT] [-lr] [-lt] [-h] FILE

Explore the contents of an ANANSE binding file.

required arguments:
  FILE                  TF binding prediction file (ANANSE binding output)

optional arguments:
  -o FILE, --outfile FILE
                        Output file (tab-separated text, default: stdout)
  -t [TF ...], --tfs [TF ...]
                        Transcription factor(s) to display (default: all)
  -r [REGION ...], --regions [REGION ...]
                        Region(s) to display (default: all)
  -F FORMAT, --format FORMAT
                        Display format: wide (n columns) or long (3 columns) (default: wide)
  -n INT                Number of regions and tfs to display (default: all)
  -lr, --list-regions   Return a list of regions
  -lt, --list-tfs       Return a list of transcription factors
  -h, --help            show this help message and exit
```