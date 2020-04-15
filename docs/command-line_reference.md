## Command-line reference

* **All the example dataset and result files are able to find at [***http://mbdata.science.ru.nl/qxu/ananse/ananse.html***](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).**

### Build transcription factor binding network

``` bash
$ ananse binding  -r data/krt_enhancer.bed \
                  -o results/binding.txt \
                  -a /data/hg38_genes.bed \
                  -g hg38 \
                  -p /data/gimme.vertebrate.v5.1.pfm
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-r, --enhancers`  
    The name of the input enhancer peak file. This should be a BED-formatted file with 4 columns. The first column is the chromosome name, the second and third columns are the start and end point of peak. We recommend that all peaks have a size of 200bp. If the peak is not 200bp, ANANSE will change it to 200bp. The fourth column is intensity of the peak, this can be the number of reads, RPKM or equivalent value. You can find the method to generate the enhancer file and an example of an enhancer input file in the section [Input data](input_data/#enhancer-data).  
* `-o, --output`  
    The name of the output file.

**Optional arguments:**  

* `-n, --ncore`  
    Specifies the number of threads to use during analysis.  
* `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. You can find the method to generate genome files in the section [Input data](input_data/#genome). The default setting is `hg19`.  
* `-a, --annotation`  
    Gene annotation for the genome specified with `-g` as a 12 column BED file. You can find the method to generate this annotation BED file and example BED files in the section [Input data](input_data/#genome).  
* `-p, --motifs`  
    The input motif file with positional frequence matrices. You can find the definition of the motif file and the default motif files in the section [Input data](input_data/#motif-database).  
* `-f, --filter_promoter`  
    Filter promoters. Default setting is True. If 'True', the function will remove all promoter peaks (+-2k from TSS) from the provided peaks.
* `-h, --help`  
    Show the help message and exit.

### Build gene regulatory network

``` bash
$ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                  -b results/binding.txt \
                  -o results/full_features.txt \
                  -a /data/hg38_genes.bed \
                  -g hg38 \
                  -p ../data/gimme.vertebrate.v5.1.pfm
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-e, --expression`  
    The expression file(s) of your cell type or tissue of interest. You can supply one or more gene expression file(s). In these files, the 1st column should contain gene name and it should contain a column with the name `tpm`. [***This***](/test/data/KRT_rep1_TPM.txt) is an example of a valid expression file.  
* `-b, --binding`  
    The binding network from the `Build binding network` step, generated with `ananse binding`. An example `binding network` can be found [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/binding.txt).  
* `-o, --output`  
    The folder to save results. 

**Optional arguments:**  

* `-n, --ncore`  
    Specifies the number of threads to use during analysis.  
* `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. You can find the method to generate genome files at [Genome](input_data/#genome) part. The default setting is `hg19`.  
* `-a, --annotation`  
    Gene annotation for the genome specified with `-g` as a 12 column BED file. You can find the method to generate this annotation BED file and example BED files in the section [Input data](input_data/#genome).   
* `-p, --motifs`  
    The input motif file with positional frequence matrices. You can find the definition of the motif file and the default motif files in the section [Input data](input_data/#motif-database).  
* `-c, --corrfiles`  
    A file with gene-gene correlation values. The human gene expression correlation can be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/data/expressioncorrelation.txt).  
* `-h, --help`  
    Show the help message and exit.  

### Infer TF influence score

``` bash
$ ananse influence  -a results/full_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -p 
```
!!! tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  

* `-t, --target`  
    The gene regulatory network of the target cell type. It is the result from `Build gene regulatory network` step, generated by `ananase network`. An `network` file can be found [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).  
* `-d, --degenes`  
    The differential expression table between two cells. [***This***](/test/data/FB2KRT_degenes.csv) is an example of differential expression file.  
* `-o, --output`  
    The folder to save results.  

**Optional arguments:**  

* `-n, --ncore`  
    Specifies the number of threads to use during analysis.  
* `-i, --edges`  
    Specifics the number of top edges (interactions) used.  
* `-s, --source`  
    The gene regulatory network of the source cell type (optional but recommended). It is the result from `Build gene regulatory network` step, generated by `ananase network`. An `network` file can be found[***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).  
* `-e, --expression`  
    The gene expression level in the source cell type (optional but recommended). One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/FB_rep1_TPM.txt) is an example of expression file.  
* `-p, --plot`  
    Plot influence. Default setting is True.  
* `-h, --help`  
    Show the help message and exit.  
