# ANANSE: ANalysis Algorithm for Networks Specified by Enhancers
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/version.svg)](https://anaconda.org/bioconda/ananse)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/downloads.svg)](https://anaconda.org/bioconda/ananse)
[![Documentation Status](https://readthedocs.org/projects/anansepy/badge/?version=master)](https://anansepy.readthedocs.io/en/master/?badge=master)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/license.svg)](https://anaconda.org/bioconda/ananse)
### Prediction of key transcription factors in cell fate determination using enhancer networks
ANANSE is a computational approach to infer enhancer-based gene regulatory networks (GRNs) and to use these GRNs to identify the key transcription factors in cell fate determination. You can use it to generate a shortlist of transcription factors for trans-differentiation experiments, but also to study transcription regulation during development and differentiation. It is written in Python and it contains a user-friendly command-line script that includes `ananse binding`, `ananse network`, and `ananse influence`.

![](/pic/Fig2.png)
<!-- > (A), Data types required and utilized in ANANSE. These data include motif score of all TFs, gene expression data (e.g. RNA-seq) and enhancer data that can be obtained by ATAC-seq, EP300 ChIP-seq, or H3K27ac ChIP-seq from each cell type. The blue and orange peaks represent enhancers in two cell types. The four sequence-logos represent the motif of four TFs. The heatmap represents gene expression intensity in two cell types. (B), The TF binding profiles predicted from enhancer data and TF motif scores in each cell type. Two GRNs below show cell type-specific TF binding profiles in two cell types (source and target cell types). (C), The cell type-specific GRN predicted based on TF-Gene binding and TF/Gene expression. Two networks show cell type-specific GRN in two cell types. The orange circle represents a TF or a gene, and the size of the circle indicates the target gene number of the corresponding TF. The blue arrow indicates regulation between two TFs, and the color intensity represents regulation intensity. (D), The differential GRN between the two cell types. In this step, the interaction specific for the target cell type is kept constant, and if the interaction score of the target cell type is higher than that of the source cell type, the interaction score is further used. (E), The barplot shows the ranked influence score of all TFs calculated from the differential GRN. The influence score is calculated based on gene expression score, distance from the enhancer bound by TF to gene, and the interaction score between TF and gene. -->

Read **[full ANANSE documentation](https://anansepy.readthedocs.io/en/master/)** for detailed installation instructions and usage examples. For documentation on the **development version** see [here](https://anansepy.readthedocs.io/en/develop/).


<!-- --- -->

## Quick start
<!-- * ### **Detail documentation**
  * The **full ANANSE documentation** at [https://anansepy.readthedocs.io](https://anansepy.readthedocs.io).  -->

* ### **Easy installation**
  <!-- * The most straightforward way to install ANANSE is via conda using the bioconda channel. -->

  #### ***1. If you have not used bioconda before, first set up the necessary channels (in this order!). You only have to do this once.***

  ```
  $ conda config --add channels defaults
  $ conda config --add channels bioconda
  $ conda config --add channels conda-forge
  ```
  #### ***2. Install ANANSE from bioconda***
  ``` 
  # Create an environment called ananse with all dependencies
  $ conda create -n ananse python=3 ananse

  # Activate the environment
  $ conda activate ananse
  ```
  <!-- * Python 3 is the required for ANANSE. Don't forget to activate the environment with conda activate gimme whenever you want to use ANANSE. -->


* ### **Usage**

  #### ***0. Make enhancer file***
  ```
  $ ananse enhancer -g hg38 -t H3K27ac \
                    -b data/KRT_H3K27ac_rep1.bam \
                    -p data/KRT_H3K27ac_peaks.broadPeak \
                    -o data/KRT_enhancer.bed 
  ```

  * `-t, --etype`. Enhancer type, H3K27ac, p300, or ATAC. **H3K27ac only provide for hg38!** 
  * `-g, --genome`. The genome of the data.
  * `-b, --bam_input`. The H3K27ac or p300 ChIP-seq bam file.
  * `-p, --epeak`. The H3K27ac ChIP-seq broadPeak, or the p300 ChIP-seq / ATAC-seq narrowPeak.
  * `-o, --bed_output`. The output enhancer file.
  * `-h, --help`. Show the help message and exit.
  
  <!-- * **All the example dataset and result files are able to find at [***http://mbdata.science.ru.nl/qxu/ananse/ananse.html***](http://mbdata.science.ru.nl/qxu/ananse/ananse.html).** -->
  <!-- --- -->
  #### ***1. Build TF binding network***  
  <!-- > Predict cell type-specific transcription factor binding with enhancer intensity and motif z-score. -->

  <!-- * Example:  -->
  ```
  $ ananse binding  -r data/KRT_enhancer.bed \
                    -o results/binding.txt \
                    -g hg38 -t H3K27ac
  ```

  * `-r, --enhancers`. The input enhancer peak file. 
  * `-t, --etype`. Enhancer type, H3K27ac, p300, or ATAC. **H3K27ac only provide for hg38!** 
  * `-o, --output`. The output file.
  * `-g, --genome`. The genome of the data.
  * `-h, --help`. Show the help message and exit.

  <!-- --- -->
  #### ***2. Built gene regulatory network***  
  <!-- > Infer cell type-specific gene regulatory network with TF binding and distance to promoter. -->

  <!-- * Example: -->
  ```
  $ ananse network  -e data/KRT_rep1_TPM.txt data/KRT_rep2_TPM.txt \
                    -b results/binding.txt \
                    -o results/KRT_full_features.txt \
                    -g hg38 \
                    --exclude-promoter --include-enhancer
  ```

  <!-- * Required arguments: -->
  * `-e, --expression`. The expression file of your interested cell type or tissue. 
  * `-b, --binding`. The binding network from `Build binding network` step. 
  * `-o, --output`. The output file. 
  * `-g, --genome`. The genome of your data. 
  * `--include-promoter, --exclude-promoter`. Include or exclude promoter peaks (<= TSS +/- 2kb) in network inference. By default promoter peaks are **excluded**.
  * `--include-enhancer, --exclude-enhancer`. Include or exclude enhancer peaks (> TSS +/- 2kb) in network inference. By default enhancer peaks are **included**.
  * `-h, --help`. Show the help message and exit.

  <!-- --- -->
  #### ***3. Infer TF influence score***  
  <!-- > Infer key TFs during cell fate determination with TF expression and gene regulatory network. -->

  <!-- * Example: -->
  ```
  $ ananse influence  -s results/FB_full_network.txt \
                      -t results/KRT_full_network.txt \
                      -d data/FB2KRT_degenes.csv \
                      -o results/FB2KRT.txt 
  ```

  <!-- * Required arguments: -->
  * `-s, --source`. The network in source cell (optional).     
  * `-t, --target`. The network in target cell.  
  * `-d, --degenes`. The differential expression table between two cells.  
  * `-o, --output`. The output file.  
  * `-h, --help`. Show the help message and exit.

<!-- ___ -->
## Citation
  > Xu Q, Georgiou G, Veenstra G J C, et al. ANANSE: An enhancer network-based computational approach for predicting key transcription factors in cell fate determination[J]. [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.06.05.135798v2), 2020.

<!-- --- -->
## Help and Support

  * The preferred way to get support is through the [Github issues page](https://github.com/vanheeringen-lab/ANANSE/issues).

  * Reach out to me at one of the following places!

    - Website at <a href="https://github.com/vanheeringen-lab" target="_blank">`vanheeringen-lab`</a>
    - Email to <a href="mailto:qxuchn@gmail.com" target="_blank">`Quan Xu`</a> or <a href="mailto:simon.vanheeringen@gmail.com" target="_blank">`Simon J. van Heeringen`</a>

<!-- --- -->

## License

  - **[MIT license](http://opensource.org/licenses/mit-license.php)** [![Anaconda-Server Badge](https://anaconda.org/qxuchn/ananse/badges/license.svg)](https://anaconda.org/qxuchn/ananse)
  - Copyright 2020 © <a href="https://github.com/vanheeringen-lab" target="_blank">vanheeringen-lab</a>.
