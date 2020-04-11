## Command-line reference

* **All the example dataset and result files are able to find at [***http://mbdata.science.ru.nl/qxu/ananse/ananse.html***](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).**

### Build TF binding network

``` bash
$ ananse binding  -r data/krt_enhancer.bed \
                  -o results/binding.txt \
                  -a /data/hg38_genes.bed \
                  -g hg38 \
                  -p /data/gimme.vertebrate.v5.1.pfm
```
''' tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  
  * `-r, --enhancers`  
    The name of the input enhancer peak file. This should be a BED format file, with 4 columns. The first column is chromosome name, the second and third columns are the start and end point of peak. We recommend all peaks have 200bp. If the peak is not 200bp, we will normize it to 200bp. The fourth column is intensity of the peak, it could be RPKM or equivalent value. [***This***](https://github.com/vanheeringen-lab/ANANSE/raw/master/test/data/krt_enhancer.bed) is an example enhancer BED file.
  * `-o, --output`  
    The name of the output file.

**Optional arguments:**  
  * `-n, --ncore`  
    Specifies the number of threads to use during analysis.  
  * `-g, --genome`  
    The genome that is used for the gene annotation and the enhancer location. 
    This can be either the name of a genome installed with [genomepy](https://github.com/vanheeringen-lab/genomepy), for example `hg38`, or the name of a genome FASTA file, for example `/data/genomes/hg38/hg38.fa`. It is recommended to use a genome installed by `genomepy`. The default setting is `hg19`.
  * `-a, --annotation`  
    The input 12 columns BED file with gene annotation in your genome version. [***This***](https://github.com/vanheeringen-lab/ANANSE/raw/master/data/hg38_genes.bed) is an example BED annotation file of human hg38.
  * `-p, --motifs`  
    The input Motif file. [***This***](/data/gimme.vertebrate.v5.1.pfm) is an example Motif file in vertebrate. if provided there should also be a motif2factors.txt file and a factortable.txt file in the same folder. [***This***](/data/gimme.vertebrate.v5.1.motif2factors.txt) is an example of motif2factors file. [***This***](/data/gimme.vertebrate.v5.1.factortable.txt) is an example of factortable file.
  * `-f, --filter_promoter`  
    Filter promoters, True or False, input should be
    either 'True' or 'False'. (Default setting: True; if 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.).
  * `-d, --keep_detail`  
    Keep detail files, True or False, input should be either 'True' or 'False'. (Default setting: True).  
  * `-h, --help`  
    Show the help message and exit.

### Built gene regulatory network

``` bash
$ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                  -b results/binding.txt \
                  -o results/full_features.txt \
                  -a /data/hg38_genes.bed \
                  -g hg38 \
                  -p ../data/gimme.vertebrate.v5.1.pfm
```
''' tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**
  * `-e, --expression`  
    The expression file of your interested cell type or tissue. It could have one or more gene expression file(s). In this file, the 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/KRT_rep1_TPM.txt) is an example of expression file.   

  * `-b, --binding`  
    The binding network from `Build binding network` step. One of the example `binding network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/binding.txt).  
  * `-o, --output`  
    The folder to save results, `-o` is the required arguments. 

**Optional arguments:**
  * `-n, --ncore`  
    Specifies the number of threads to use during analysis. 
  * `-g, --genome`  
    The genome of your data. For example, hg38. The genome is recommended to download by `genomepy`.
  * `-p, --motifs`  
    The input Motif file. [***This***](/data/gimme.vertebrate.v5.1.pfm) is an example Motif file in vertebrate. if provided there should also be a motif2factors.txt file and a factortable.txt file in the same folder. [***This***](/data/gimme.vertebrate.v5.1.motif2factors.txt) is an example of motif2factors file. [***This***](/data/gimme.vertebrate.v5.1.factortable.txt) is an example of factortable file. 
  * `-a, --annotation`  
    The input 12 columns BED file with gene annotation in your genome version. [***This***](/data/hg38_genes.bed) is an example BED annotation file of human hg38.
  * `-f, --filter_promoter`  
    Filter promoters, True or False, input should be
    either 'True' or 'False'. (Default setting: True; if 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.).
  * `-c, --corrfiles`  
    All gene correlation file, the human gene expression correlation can be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/data/expressioncorrelation.txt).
  * `-h, --help`  
    Show the help message and exit.

### Infer TF influence score

``` bash
$ ananse influence  -a results/full_network.txt \
                    -e data/FB_rep1_TPM.txt \
                    -d data/FB2KRT_degenes.csv \
                    -o results/FB2KRT.txt \
                    -p False
```
''' tip
    Please use `-h/--help` for the details of all options.

**Required arguments:**  
  * `-a, --anetwork`  
  The network in second cell. It is the result from `Built GRN` step. One of the example `network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).   
  * `-d, --degenes`  
  The differential expression table between two cells. [***This***](/test/data/FB2KRT_degenes.csv) is an example of differential expression file.  
  * `-o, --output`  
  The folder to save results, `-o` is the required arguments.   

**Optional arguments:**  
  * `-n, --ncore`  
    Specifies the number of threads to use during analysis. 
  * `-s, --edges`  
    Specifics the number of top edges (interactions) used. 
  * `-b, --bnetwork`  
  The network in first cell (optional). It is the result from `Built GRN` step. One of the example `network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).  
  * `-e, --expression`  
  The gene expression in first cell (optional). One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/FB_rep1_TPM.txt) is an example of expression file. 
  * `-p, --plot`  
  Plot influence. True or False, input should be either 'True' or 'False'. (Default setting: True)  
  * `-h, --help`  
  Show the help message and exit.
