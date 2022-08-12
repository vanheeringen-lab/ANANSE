## Input data

To run ANANSE you need the following data:

* A genome with a matching gene annotation
* For `ananse binding`: enhancer regions (optional for `hg38`)
* For `ananse binding`: enhancer activity from either/both:
    * ATAC-seq data
        * one or more indexed BAM file(s) or
        * one counts table with reads for each peak per sample
    * H3K27ac ChIP-seq data
        * one or more indexed BAM files(s) or
        * one counts table with reads per peak
* For `ananse network`: gene expression quantification from either:
    * one or more quantification files with one TPM column
    * one counts/TPM table with expression for each gene per sample
* For `ananse influence`: gene differential expression (DESeq2 output)
* A Motif database

### Genome

To run `ANANSE` on your sample(s), a matching genome and gene annotation is necessary.

#### genomepy

We recommend that you download these with [genomepy](https://github.com/vanheeringen-lab/genomepy).
If you installed ANANSE with conda, genomepy is present in the ANANSE environment.

To install a genome with `genomepy` you can use this command, 
which will download both the genome `FASTA` and the gene annotation `BED` and `GTF` files.

```shell
# activate ananse environment
conda activate ananse

# install the human genome from UCSC
genomepy install hg38 --provider UCSC --annotation

# alternatively, install the human genome from Ensembl
genomepy install GRCh38.p13 --provider Ensembl --annotation
```

You can then use the `-g` argument to specify your genome by name when running the ANANSE tools: `-g hg38`.

#### I have a genome

Alternatively, you can specify the genome manually. In this case you'll need:

* A genome FASTA file.
* A gene annotation in BED12 format.
  * the chromosome names must match with the genome
  * the gene names must be HGNC symbols

You can then specify the genome as follows: `-g /path/to/genome.fa`.
For `ananse network`, you need to add the annotation BED file with `-a /path/to/annotation.bed`.

You can use genomepy to copy these files into a genomepy-install with
```shell
genomepy install -p local /path/to/genome.fa --Local-path-to-annotation /path/to/annotation.gtf
```

### Enhancer regions

The `ananse binding` command requires a set of putative enhancer regions.

For human data (`hg38` only), ANANSE provides a [database of cis-regulatory regions](https://doi.org/10.5281/zenodo.4066424) (putative enhancers).
Which is used by default. You can still specify your own set of putative enhancers.

To define enhancer regions, you can use any type of genome-wide data that is associated with enhancers and *gives narrow peaks*.
H3K27ac signal, for instance, would not work well, as peaks from a H3K27ac ChIP-seq experiment are too broad to provide a precise region for the motif analysis.
Examples of data that is suitable would be DNaseI, ATAC-seq or EP300 ChIP-seq.
The enhancer regions should be specified as a BED file, centered at the summit.
You can also provide one or more `narrowPeak` files, for instance from [MACS2](https://github.com/taoliu/MACS).
In this case all the peaks will be merged and overlapping peaks will be centered at the summit of the highest peak.
If you analyze multiple samples with ANANSE, you should include all the peaks from every sample, thereby creating a union of all the peaks in your experiment.

You can use [seq2science](https://github.com/vanheeringen-lab/seq2science) to map your data.
The (raw) ATAC-/ChIP-seq counts table can be used as input directly.

### Enhancer activity

ANANSE can use ATAC-seq and/or H3K27ac ChIP-seq data to predict binding.
Using both will give the best performance, however, either one will also work well (see Fig. 3A in the ANANSE paper).
These data should be mapped to the relevant genome, with duplicates reads marked (or removed).
For both types of data you can supply one or more files (replicates), which will be averaged.

For each data type, enhancer activity can be given in two forms: one or more BAMs or one counts table.

###### BAMs

The BAM file(s) should be sorted indexed, for instance with `samtools sort` and `samtools index`.

###### Counts

A counts table must be a tab-separated file with peaks in the first column, and reads per peak for each sample.

You can use [seq2science](https://github.com/vanheeringen-lab/seq2science) to map your data.
The (raw) ATAC-seq counts table can be used as input directly.
The (raw) H3K27Ac ChIP-seq counts table can also be used directly, but it is suggested that peaks are merged and centered on a 2000 bp window.
This can be automated by adding the following options to the seq2science ChIP-seq config.yaml:
```yaml
# suggested H3K27Ac ChIP-seq settings:
slop: 1000
peak_windowsize: 1000
```

Example of a counts table:

|               | sample1 | sample2 | sample3 |
| -----------   | ------- | ------- | ------- |
| 9:2802-3002   | 8.0     | 6.0     | 6.0     |
| 9:3612-3812   | 5.0     | 12.0    | 5.0     |
| 9:16114-16314 | 2.0     | 12.0    | 4.0     |

You can specify which samples to use from this file with `--columns` (case-insensitive. By default, all columns are used). 
For example: `--columns SAMPLE1 sample3` will use samples 1 and 3, but ignores sample 2.

### Expression data

Expression data normally comes from an RNA-seq experiment.
The file(s) must be tab-separated, with genes in the first column.
We suggest using `gene_name`s for genes, and expression scores in `TPM`.

Notes: 

* If you are using a genomepy genome, genes are automatically converted if required.
  * So you could use genes described as `transcript_id`, `transcript_name` or `gene_id` as well.

Expression data can be given in two forms: one or more quantification files or one (TPM) counts table.

###### Counts

A counts table must be a tab-separated file with genes in the first column, and expression for each gene per sample.

Example of a counts table:

|        | sample1 | sample2 | sample3 |
| -----  | ------- | ------- | ------- |
| A1BG   | 8.0     | 6.0     | 6.0     |
| A1CF   | 5.0     | 12.0    | 5.0     |
| A2ML1  | 2.0     | 12.0    | 4.0     |

You **must** specify which samples to use from this file with `--columns` (case-insensitive. By default, column `tpm` is used). 
For example: `--columns SAMPLE1 sample3` will use samples 1 and 3, but ignores sample 2.

###### quantification file(s)

Quantification files are generated by pseudo-aligners such as salmon or kallisto, 
followed by summing transcript-level TPMs to gene-level TPMS using [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html).

Example of a quantification file:

|           | tpm      |
| --------- | -------- |
| A1BG      | 9.579771 |
| A1CF      | 0.0223   |
| A2ML1     | 664.452  |

If your files have a different column name for gene expression levels, you can specify this with `--columns` (case-insensitive. By default, column `tpm` is used). 
For example: `--columns TpM` will use the `tpm` column. This may be useful if your files have a different column name (e.g. `counts`)

### Differential expression data

The differential expression data normally comes from an RNA-seq experiment.
This file can be created using, for instance, [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

The file must be tab-separated, with genes (we suggest using `gene_name`s) in the first column, 
a column named `padj` and a column named `log2FoldChange`.
Column `log2FoldChange` is the **log2 of the fold change** of each gene between conditions.
Column `padj` is the **adjusted p-value** of the differential gene expression test.

Notes: 

* If you are using a genomepy genome, genes are automatically converted if required.
  * So you could use genes described as `transcript_id`, `transcript_name` or `gene_id` as well.
* The `log2FoldChange` should be a **positive number** if this gene is upregulated from the source to the target cell type, and **negative number** if this gene is downregulated from the source to the target cell type.

Example of a differential expression file:

|           | log2FoldChange    | padj   |
| --------- | ----------------- | ------ |
| ANPEP     | 7.44242618323665  | 0.001  |
| CD24      | -8.44520139575174 | 0      |
| COL1A2    | 8.265875689393    | 0.0123 |

### Motif database

By default ANANSE uses a non-redundant, clustered database of known vertebrate motifs: `gimme.vertebrate.v5.0`.
These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources.
[Large-scale benchmarks](https://www.biorxiv.org/content/10.1101/474403v1.full) using ChIP-seq peaks show that this database shows good performance and should be a good default choice.

Alternatively, you can use any of the other motif databases included with [GimmeMotifs](https://gimmemotifs.readthedocs.io/en/master/overview.html#motif-databases) (such as `JASPAR2020_vertebrates` or `HOMER`).
If you would like to use your own motif database, please make sure you create the following two files:
1) a motif file and 2) a file containing mapping of motifs to factors.
The motif file should contain positional frequency matrices and should end with the `.pfm` extension.
The motif mapping file should have the same base name as  the `motif` file and end with `.motif2factors.txt` instead of `.pfm`.
Examples are shown below.

**Please note:** if you want to use ANANSE for a species other than human or mouse, you will have to generate this database yourself, because the gene names with the database won't overlap.
Don't worry, it is easy when doing it with gimmemotifs' in-built command `gimme motif2factors`.
You can already generate a motif database for e.g. zebrafish (danRer11 assembly) with just this command:

```
gimme motif2factors --new-reference danRer11 --outdir mynewdatabase
```

Take a look at the [gimme motif2factors](https://gimmemotifs.readthedocs.io/en/master/reference.html#command-gimme-motif2factors) docs for more options

#### Motif file

```    
# Comments are allowd
>GM.5.0.Sox.0001
0.7213	0.0793	0.1103	0.0891
0.9259	0.0072	0.0062	0.0607
0.0048	0.9203	0.0077	0.0672
0.9859	0.0030	0.0030	0.0081
0.9778	0.0043	0.0128	0.0051
0.1484	0.0050	0.0168	0.8299
```

#### Motif2factors file  

```
Motif	Factor	Evidence	Curated
GM.5.0.Sox.0001	SRY	JASPAR	Y
GM.5.0.Sox.0001	SOX9	Transfac	Y
GM.5.0.Sox.0001	Sox9	Transfac	N
GM.5.0.Sox.0001	SOX9	SELEX	Y
GM.5.0.Sox.0001	Sox9	SELEX	N
GM.5.0.Sox.0001	SOX13	ChIP-seq	Y
GM.5.0.Sox.0001	SOX9	ChIP-seq	Y
GM.5.0.Sox.0001	Sox9	ChIP-seq	N
GM.5.0.Sox.0001	SRY	SELEX	Y
```
