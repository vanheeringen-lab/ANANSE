# ANANSE: ANalysis Algorithm for Networks Specified by Enhancers
[![bioconda-badge](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/version.svg)](https://anaconda.org/bioconda/ananse)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/downloads.svg)](https://anaconda.org/bioconda/ananse)

[![Documentation Status](https://readthedocs.org/projects/anansepy/badge/?version=master)](https://anansepy.readthedocs.io/en/master/?badge=master)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/ananse/badges/license.svg)](https://anaconda.org/bioconda/ananse)
[![DOI:10.1093/nar/gkab598](http://img.shields.io/badge/DOI-10.1093/nar/gkab598-B31B1B.svg)](https://doi.org/10.1093/nar/gkab598)

[![Maintainability](https://api.codeclimate.com/v1/badges/875df8c40fec66d68b1f/maintainability)](https://codeclimate.com/github/vanheeringen-lab/ANANSE/maintainability)
[![Test Coverage](https://api.codeclimate.com/v1/badges/875df8c40fec66d68b1f/test_coverage)](https://codeclimate.com/github/vanheeringen-lab/ANANSE/test_coverage)

### Prediction of key transcription factors in cell fate determination using enhancer networks
ANANSE is a computational approach to infer enhancer-based gene regulatory networks (GRNs) and to identify key transcription factors between two GRNs. You can use it to study transcription regulation during development and differentiation, or to generate a shortlist of transcription factors for trans-differentiation experiments. 

ANANSE is written in Python and comes with a command-line interface that includes 3 main commands: `ananse binding`, `ananse network`, and `ananse influence`. A graphical overview of the tools is shown below.

![](docs/img/Fig2.png)

Check out the **[ANANSE documentation](https://anansepy.readthedocs.io/en/master/)** for 
* [installation instructions](https://anansepy.readthedocs.io/en/master/installation/)
* [command explanations](https://anansepy.readthedocs.io/en/master/command-line_reference/)
* [input explanations and examples](https://anansepy.readthedocs.io/en/master/input_data/)
* [usage examples](https://anansepy.readthedocs.io/en/master/examples/)
* [FAQ](https://anansepy.readthedocs.io/en/master/faq/)
* and more!
 
For documentation on the **development version** see [here](https://anansepy.readthedocs.io/en/develop/).

## Citation

> ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination 
> Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Gert Jan C Veenstra, Huiqing Zhou, Simon J van Heeringen
> Nucleic Acids Research, gkab598, https://doi.org/10.1093/nar/gkab598


## Help and Support

* The preferred way to get support is through the [Github issues page](https://github.com/vanheeringen-lab/ANANSE/issues).

## scANANSE: Gene regulatory network and motif analysis of single-cell clusters

Recently a pipeline was developed to run ANANSE using single-cell RNA- sequencing data and single-cell ATAC-sequencing data. It consists of packages to export single-cell cluster data from Seurat or Scanpy objects, and a snakemake workflow. Afterwards, results can be imported back into your single-cell object.

For more info on this implementation  check out the
* [scANANSE workflow](https://doi.org/10.12688/f1000research.130530.1)
* [Python package for Scanpy objects](https://github.com/Arts-of-coding/AnanseScanpy)
* [R package for Seurat objects](https://github.com/JGASmits/AnanseSeurat)
* [ANANSE snakemake workflow](https://github.com/vanheeringen-lab/anansnake)

## License

  - **[MIT license](http://opensource.org/licenses/mit-license.php)** [![Anaconda-Server Badge](https://anaconda.org/qxuchn/ananse/badges/license.svg)](https://anaconda.org/qxuchn/ananse)
  - Copyright 2020 © <a href="https://github.com/vanheeringen-lab" target="_blank">vanheeringen-lab</a>.
