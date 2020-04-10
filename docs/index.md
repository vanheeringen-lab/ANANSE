## **ANANSE**: **AN**alysis **A**lgorithm for **N**etworks **S**pecified by **E**nhancers

### What is ANANSE Network?
ANANSE is an analysis framework for key transcription factor prediction during cell fate switch written in Python. It contains command-line scripts to predict all TF binding network (`binding`), predict gene regulatory network (`network`) and infer influence score between two cell types(tissues) (`influence`). In addition, all this functionality is available from a Python API.


### Citation
ANANSE is free and open source research software. If you find it useful please cite our paper:
> Xu, Q., Georgiou, G., Veenstra, G.J.C., Zhou, H., and van Heeringen, S.J. (2020). **ANANSE: An enhancer network-based computational approach for predicting key transcription factors in cell fate determination.** [https://www.biorxiv.org/](https://www.biorxiv.org/)

### Getting started
* The easiest way to install ANANSE is using [bioconda](https://bioconda.github.io/) on Linux (include [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) ) or Mac. 
* Have a look at these simple [examples](examples.md) to get a taste of what is possible.
* Check out the more detailed [command-line usage tutorials](command-line_reference.md).
* There’s also an [python API documentation](API_documentation.md) for python users.

### Get help
* First, check the FAQ for common issues.
* The preferred way to get support is through the [Github issues page](https://github.com/vanheeringen-lab/ANANSE/issues).
* Finally, you can reach me by email to <a href="mailto:qxuchn@gmail.com" target="_blank">Quan Xu</a> or <a href="mailto:simon.vanheeringen@gmail.com" target="_blank">Simon J. van Heeringen</a>.
* You can also visit our website at <a href="https://github.com/vanheeringen-lab" target="_blank">vanheeringen-lab</a>

### Full contents
* [Installation](installation.md)
    - [The easiest way to install](installation/#the-easiest-way-to-install)
    - [Alternative installation](installation/#alternative-installation)
* [Input data](input_data.md)
    - [Genome](input_data/#genome)
    - [Motif database](input_data/#motif-database)
    - [Enhancer data](input_data/#enhancer-data)
    - [Expression data](input_data/#enhancer-data)
    - [Differential expression data](input_data/#enhancer-data)
* [Command-line reference](command-line_reference.md)
    - [Command: ananse binding](command-line_reference/#motif-database)
    - [Command: ananse network](command-line_reference/#motif-database)
    - [Command: ananse influence](command-line_reference/#motif-database)
* [API documentation](API_documentation.md)
    - [Working with `binding` Class](API_documentation/#working-with-binding-class)
    - [Working with `network` Class](API_documentation/#working-with-network-class)
    - [Working with `influence` Class](API_documentation/#working-with-influence-class)
* [Examples](examples.md)
    - [Command: ananse binding](examples/#motif-database)
    - [Command: ananse network](examples/#motif-database)
    - [Command: ananse influence](examples/#motif-database)
* FAQ
    - Sorry, motif prediction tool [X] is not supported
    - I get motifs that have differential scores in gimme maelstrom, however, the number is not different across clusters
    - I have upgraded GimmeMotifs and now it doesn’t find my genome
    - I cannot run gimme index anymore
    - I get ‘RuntimeError: Invalid DISPLAY variable’
    - I get a KeyError when running gimme maelstrom
* Acknowledgments
