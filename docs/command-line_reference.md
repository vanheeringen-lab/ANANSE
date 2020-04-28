## Command-line reference

* **All the example input data and result files can be found here:**
  [http://mbdata.science.ru.nl/qxu/ananse/ananse.html](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).

In general, an analysis with ANANSE will consist of the following steps:

1. Generate a binding network for a *target* cell type using `ananse binding`.
2. Generate a gene regulatory network (GRN) for a *target* cell type using `ananse network`, based on the binding network from 1.
3. Generate a binding network for a *source* cell type using `ananse binding`.
4. Generate a GRN for a *source* cell type using `ananse network`, based on the binding network from 3.
5. Run `ananse influence` based on the GRN of the *source* cell type (step 2) and the GRN of the *target* cell type (step 4)

### Build transcription factor binding network: ananse binding

The following command combines genome-wide enhancer intensities (EP300 ChIP-seq, H3K27ac ChIP-seq, ATAC-seq) with sequence features in enhancer peaks to infer cell type-specific TF binding profiles.

Example command:

``` bash
$ ananse binding  -r data/krt_enhancer.bed \
                  -o results/binding.txt \
                  -g hg38 \
                  -p data/gimme.vertebrate.v5.1.pfm
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-r, --enhancers`  
    The name of the input enhancer peak file. This should be a BED-formatted file with 4 columns. The first column is the chromosome name, the second and third columns are the start and end point of peak. We recommend that all peaks have a size of 200bp. If the peak is not 200bp, ANANSE will change it to 200bp. The fourth column is intensity of the peak, this can be the number of reads, RPKM or equivalent value. You can find the method to generate the enhancer file and an example of an enhancer input file in the section [Input data](input_data/#enhancer-data).  
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
                  -b results/binding.txt \
                  -o results/full_features.txt \
                  -a data/hg38_genes.bed \
                  -g hg38
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
$ ananse influence  -t results/full_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -p 
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
