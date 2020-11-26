## Command-line reference

* **All the example input data and result files can be found here:**
  [http://mbdata.science.ru.nl/qxu/ananse/ananse.html](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).

In general, an analysis with ANANSE will consist of the following steps:

1. Generate ANANSE input enhancer file for *target* cell type using `ananse enhancer`.
2. Generate binding network for *target* cell type using `ananse binding`, based on the enhancer file from step 1.
3. Generate gene regulatory network (GRN) for *target* cell type using `ananse network`, based on the binding network from step 2.
4. Generate ANANSE input enhancer file for *source* cell type using `ananse enhancer`.
5. Generate a binding network for a *source* cell type using `ananse binding`, based on the enhancer file from step 4.
6. Generate a GRN for a *source* cell type using `ananse network`, based on the binding network from step 5.
7. Run `ananse influence` based on the GRN of the *source* cell type (step 4) and the GRN of the *target* cell type (step 6)

### Establish ANANSE input enhancer file: ananse enhancer

The following command generate KRT enhancer file from KRT H3K27ac ChIP-seq BAM file and KRT H3K27ac ChIP-seq BoardPeak file (***Only for hg38***). It is also possible generate this enhancer file from 1) p300 ChIP-seq BAM file and p3000 ChIP-seq narrowPeak file; 2) H3K27ac BAM file and ATAC-seq narrowPeak file (***for other genome***).

Example command:

``` bash
$ ananse enhancer -g hg38 -t hg38H3K27ac \
                  -b data/KRT_H3K27ac.sorted.bam \
                  -p data/KRT_H3K27ac.broadPeak \
                  -o data/KRT_enhancer.bed
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-t, --etype`  
    Enhancer type, hg38H3K27ac, p300, or ATAC. If you would like to run ANANSE in human data, we recommend you using hg38 genome and H3k27ac data as enhancer type. And this **hg38H3K27ac*** option ***only provide for hg38!** For other genome or human data does not have H3K27ac, you can set `-t` to `p300` or `ATAC`. 

!!! note 
    There is 3 type of enhancer data: `hg38H3K27ac`, `p300`, or `ATAC`.  
    
    * For human with H3K27ac ChIP-seq data, using `hg38H3K27ac`: 1, hg38 genome; 2, H3K27ac ChIP-seq BAM file; 3, H3K27ac ChIP-seq BoardPeak file.  
    * For p300 ChIP-seq data, using `p300`: 1, p300 ChIP-seq BAM file; 2, p3000 ChIP-seq narrowPeak file.  
    * For ATAC-seq data, using `ATAC`: 1, H3K27ac BAM file; 2, ATAC-seq narrowPeak file.  

* `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. You can find the method to generate genome files in the section [Input data](input_data/#genome). The default genome is `hg38`.   * `-b, --bam_input`. The H3K27ac or p300 ChIP-seq bam file.
* `-p, --epeak`  
    The H3K27ac ChIP-seq broadPeak, or the p300 ChIP-seq / ATAC-seq narrowPeak.
* `-o, --bed_output`  
    The output enhancer file.

**Optional arguments:**  

* `-h, --help`  
    Show the help message and exit.


### Build transcription factor binding network: ananse binding

The following command combines genome-wide enhancer intensities (EP300 ChIP-seq, H3K27ac ChIP-seq, ATAC-seq) with sequence features in enhancer peaks to infer cell type-specific TF binding profiles.

Example command:

``` bash
$ ananse binding  -r data/KRT_enhancer.bed \
                  -o results/KRT_binding.txt \
                  -g hg38 -t hg38H3K27ac \
                  -p data/gimme.vertebrate.v5.1.pfm
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-r, --enhancers`  
    The name of the input enhancer peak file. This should be a BED-formatted file with 4 columns. The first column is the chromosome name, the second and third columns are the start and end point of peak. We recommend that all peaks have a size of 200bp. If the peak is not 200bp, ANANSE will change it to 200bp. The fourth column is intensity of the peak, this can be the number of reads, RPKM or equivalent value. You can find the method to generate the enhancer file and an example of an enhancer input file in the section [Input data](input_data/#enhancer-data).  
* `-t, --etype`  
    Enhancer type, hg38H3K27ac, p300, or ATAC.  
* `-o, --output`  
    The name of the output file.

**Optional arguments:**  

* `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. You can find the method to generate genome files in the section [Input data](input_data/#genome). The default genome is `hg38`.   
* `-p, --pfmfile`  
    The input motif file with positional frequence matrices. You can find the definition of the motif file and the default motif files in the section [Input data](input_data/#motif-database).
* `-n, --ncore`  
    Specifies the number of threads to use during analysis.   
* `-h, --help`  
    Show the help message and exit.

### Build gene regulatory network: ananse network

This command infers cell type-specific GRNs based on the predicted TF binding sites using `ananse binding` and the expression levels of both TFs as well as their target genes. TF-gene interaction scores, the edge weights in the network, are calculated based on the predicted TF binding probability, the distance between the enhancer and the target gene, and expression of both TF and the target gene.

Example command:

``` bash
$ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                  -b results/KRT_binding.txt \
                  -o results/KRT_network.txt \
                  -a data/hg38_genes.bed \
                  -g hg38 \
                  --exclude-promoter --include-enhancer
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-b, --binding`  
    The binding network from the `Build binding network` step, generated with `ananse binding`. An example `binding network` can be found [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/binding.txt).  
* `-e, --expression`  
    The expression file(s) of your cell type or tissue of interest. You can supply one or more gene expression file(s). In these files, the 1st column should contain gene name and it should contain a column with the name `tpm`. [***This***](/test/data/KRT_rep1_TPM.txt) is an example of a valid expression file.  
* `-o, --output`  
    The folder to save results. 

**Optional arguments:**  
 
* `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. You can find the method to generate genome files at [Genome](input_data/#genome) part. The default genome is `hg38`.  
* `-a, --annotation`  
    Gene annotation for the genome specified with `-g` as a 12 column BED file. You can find the method to generate this annotation BED file and example BED files in the section [Input data](input_data/#genome).    
* `-c, --corrfiles`  
    A file with gene-gene correlation values. The human gene expression correlation can be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/data/expressioncorrelation.txt). 
* `--include-promoter, --exclude-promoter`
    Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference. By default promoter peaks are **excluded**.
  * `--include-enhancer, --exclude-enhancer`
    Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference. By default enhancer peaks are **included**.
* `-n, --ncore`  
    Specifies the number of threads to use during analysis. 
* `-h, --help`  
    Show the help message and exit.  

### Infer TF influence score: ananse influence

To calculate the influence score for the transition from a *source* cell type (`-s` or `--source`) to a *target* cell type (`t` or `--target`), `ananse influence` uses the GRNs for both cell types, predicted by `ananse network`. For each network, the top 100k interactions are selected, based on the rank of the interaction scores (edge weights). Using the differential GRN, capturing the difference between the two networks, a local network is built for each TF, up to a maximal number of three edges. Using this network, the influence score is calculated based on 1) the edge distance from the TF of interest to the target gene, 2) the predicted interaction score and 3) the change in expression between the source cell type and the target cell type.

Example command: 

``` bash
$ ananse influence  -s results/FB_network.txt \
                    -t results/KRT_network.txt \
                    -d data/FB2KRT_degenes.tsv \
                    -e data/FB_TPM.tsv \
                    -o results/influence.txt \
                    -n 20
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-t, --target`  
    The gene regulatory network of the *target* cell type. It is the result from `Build gene regulatory network` step, generated by `ananase network`. An `network` file can be found [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).  
* `-d, --degenes`  
    The differential expression table between source and target cell type. [***This***](/test/data/FB2KRT_degenes.csv) is an example of differential expression file.  
* `-o, --output`  
    The folder to save results.  

**Optional arguments:**  

* `-s, --source`  
    The gene regulatory network of the *source* cell type (optional but recommended). It is the result from `Build gene regulatory network` step, generated by `ananase network`. An example  `network` file can be found[***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).   
* `-i, --edges`  
    Specifies the number of top edges (interactions) used.  
* `-e, --expression`  
    The gene expression level in the *source* cell type. One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/FB_rep1_TPM.txt) is an example of expression file. This file is used to filter out highly expressed transcription factors in the source cell type.
* `-p, --plot`  
    Plot influence. Default setting is True.  
* `-n, --ncore`  
    Specifies the number of threads to use during analysis. 
* `-h, --help`  
    Show the help message and exit.  
