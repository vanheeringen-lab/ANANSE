# ANANSE: ANalysis Algorithm for Networks Specified by Enhancers
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/version.svg)](https://anaconda.org/bioconda/ananse)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/license.svg)](https://anaconda.org/bioconda/ananse)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/downloads.svg)](https://anaconda.org/bioconda/ananse)
[![Documentation Status](https://readthedocs.org/projects/anansepy/badge/?version=latest)](https://anansepy.readthedocs.io/en/latest/?badge=latest)

### Prediction of key transcription factors in cell fate determination using enhancer networks
Read **[full ANANSE documentation](https://anansepy.readthedocs.io/en/latest/)** for detailed installation instructions and usage examples.  



![](/pic/Fig2.jpg)
> (A), Data types required and utilized in ANANSE. These data include motif score of all TFs, gene expression data (e.g. RNA-seq) and enhancer data that can be obtained by ATAC-seq, EP300 ChIP-seq, or H3K27ac ChIP-seq from each cell type. The blue and orange peaks represent enhancers in two cell types. The four sequence-logos represent the motif of four TFs. The heatmap represents gene expression intensity in two cell types. (B), The TF binding profiles predicted from enhancer data and TF motif scores in each cell type. Two GRNs below show cell type-specific TF binding profiles in two cell types (source and target cell types). (C), The cell type-specific GRN predicted based on TF-Gene binding and TF/Gene expression. Two networks show cell type-specific GRN in two cell types. The orange circle represents a TF or a gene, and the size of the circle indicates the target gene number of the corresponding TF. The blue arrow indicates regulation between two TFs, and the color intensity represents regulation intensity. (D), The differential GRN between the two cell types. In this step, the interaction specific for the target cell type is kept constant, and if the interaction score of the target cell type is higher than that of the source cell type, the interaction score is further used. (E), The barplot shows the ranked influence score of all TFs calculated from the differential GRN. The influence score is calculated based on gene expression score, distance from the enhancer bound by TF to gene, and the interaction score between TF and gene.

---

## Quick start
* ### **Detail documentation**
  * The **full ANANSE documentation** at [https://anansepy.readthedocs.io](https://anansepy.readthedocs.io). 

* ### **Easy installation**
  * The most straightforward way to install ANANSE is via conda using the bioconda channel.

  * If you have not used bioconda before, first set up the necessary channels (in this order!). You only have to do this once.

    ```
    $ conda config --add channels defaults
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    ```
  * You can now install ANANSE with one command:
    ``` 
    # Create an environment called ananse with all dependencies
    $ conda create -n ananse python=3 ananse

    # Activate the environment
    $ conda activate ananse
    ```
  * Python 3 is the required for ANANSE. Don't forget to activate the environment with conda activate gimme whenever you want to use ANANSE.


* ### **Usage**
  
  * **All the example dataset and result files are able to find at [***http://mbdata.science.ru.nl/qxu/ananse/ananse.html***](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).**
  ---
  > ### ***Build TF binding network***  
  > Predict cell type-specific transcription factor binding with enhancer intensity and motif z-score.

  * Example:
    ```
    $ ananse binding  -r data/krt_enhancer.bed \
                      -o results/binding.txt \
                      -a /data/hg38_genes.bed \
                      -g hg38 \
                      -p /data/gimme.vertebrate.v5.1.pfm
    ```

  * Required arguments:
    * `-r, --enhancers`  
      The name of the input enhancer peak file. This should be a BED format file, with 4 columns. The first column is chromosome name, the second and third columns are the start and end point of peak. We recommend all peaks have 200bp. If the peak is not 200bp, we will normize it to 200bp. The fourth column is intensity of the peak, it could be RPKM or equivalent value. [***This***](https://github.com/vanheeringen-lab/ANANSE/raw/master/test/data/krt_enhancer.bed) is an example enhancer BED file.
    * `-o, --output`  
      The name of the output file.


  * Optional arguments:

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
      Filter promoters. Default setting is True If 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.
    * `-d, --keep_detail`  
      Keep detail files. Default setting is True.  
    * `-h, --help`  
      Show the help message and exit.

  ---
  > ### ***Built gene regulatory network***  
  > Infer cell type-specific gene regulatory network with TF binding and distance to promoter.

  * Example:
    ```
    $ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                      -b results/binding.txt \
                      -o results/full_features.txt \
                      -a /data/hg38_genes.bed \
                      -g hg38 \
                      -p ../data/gimme.vertebrate.v5.1.pfm
    ```

  * Required arguments:
    * `-e, --expression`  
      The expression file of your interested cell type or tissue. It could have one or more gene expression file(s). In this file, the 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/KRT_rep1_TPM.txt) is an example of expression file.   

    * `-b, --binding`  
      The binding network from `Build binding network` step. One of the example `binding network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/binding.txt).  
    * `-o, --output`  
      The folder to save results, `-o` is the required arguments. 

  * Optional arguments:
    * `-n, --ncore`  
      Specifies the number of threads to use during analysis. 
    * `-g, --genome`  
      The genome of your data. For example, hg38. The genome is recommended to download by `genomepy`.
    * `-p, --motifs`  
      The input Motif file. [***This***](/data/gimme.vertebrate.v5.1.pfm) is an example Motif file in vertebrate. if provided there should also be a motif2factors.txt file and a factortable.txt file in the same folder. [***This***](/data/gimme.vertebrate.v5.1.motif2factors.txt) is an example of motif2factors file. [***This***](/data/gimme.vertebrate.v5.1.factortable.txt) is an example of factortable file. 
    * `-a, --annotation`  
      The input 12 columns BED file with gene annotation in your genome version. [***This***](/data/hg38_genes.bed) is an example BED annotation file of human hg38.
    * `-f, --filter_promoter`  
      Filter promoters. Default setting is True. If 'True', the function will filtered all promoter peaks (+-2k from TSS) in provided enhancer peaks.
    * `-c, --corrfiles`  
      All gene correlation file, the human gene expression correlation can be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/data/expressioncorrelation.txt).
    * `-h, --help`  
      Show the help message and exit.

  ---
  > ### ***Infer TF influence score***  
  > Infer key TFs during cell fate determination with TF expression and gene regulatory network.

  * Example:
    ```
    $ ananse influence  -t results/full_network.txt \
                        -e data/FB_rep1_TPM.txt \
                        -d data/FB2KRT_degenes.csv \
                        -o results/FB2KRT.txt \
                        -p 
    ```

  * Required arguments:
  
    * `-t, --target`  
    The network in second cell. It is the result from `Built GRN` step. One of the example `network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).   
    * `-d, --degenes`  
    The differential expression table between two cells. [***This***](/test/data/FB2KRT_degenes.csv) is an example of differential expression file.  
    * `-o, --output`  
    The folder to save results, `-o` is the required arguments.   

  * Optional arguments:

    * `-n, --ncore`  
      Specifies the number of threads to use during analysis. 
    * `-i, --edges`  
      Specifics the number of top edges (interactions) used. 
    * `-s, --source`  
    The network in first cell (optional). It is the result from `Built GRN` step. One of the example `network` could be found at [***here***](http://mbdata.science.ru.nl/qxu/ananse/results/full_network.txt).  
    * `-e, --expression`  
    The gene expression in first cell (optional). One or more gene expression file(s), 1st column should contain gene name, and a column should be named TPM. [***This***](/test/data/FB_rep1_TPM.txt) is an example of expression file. 
    * `-p, --plot`  
    Plot influence. Default setting is True.
    * `-h, --help`  
    Show the help message and exit.


* ### **API documentation**

  * The ***python API documentation*** of this package can be found at here:  
    * [`Binding` class](/docs/api_binding.md).
    * [`Network` class](/docs/api_network.md).
    * [`Influence` class](/docs/api_influence.md).

---
## Help

  * The preferred way to get support is through the Github issues page.

---

## Support

  Reach out to me at one of the following places!

  - Website at <a href="https://github.com/vanheeringen-lab" target="_blank">`vanheeringen-lab`</a>
  - Email to <a href="mailto:qxuchn@gmail.com" target="_blank">`Quan Xu`</a> or <a href="mailto:simon.vanheeringen@gmail.com" target="_blank">`Simon J. van Heeringen`</a>

---

## License

  - **[MIT license](http://opensource.org/licenses/mit-license.php)** [![Anaconda-Server Badge](https://anaconda.org/qxuchn/ananse/badges/license.svg)](https://anaconda.org/qxuchn/ananse)
  - Copyright 2020 © <a href="https://github.com/vanheeringen-lab" target="_blank">vanheeringen-lab</a>.
