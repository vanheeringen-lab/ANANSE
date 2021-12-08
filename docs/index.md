## **ANANSE**: **AN**alysis **A**lgorithm for **N**etworks **S**pecified by **E**nhancers

### What is ANANSE?

ANANSE is a computational approach to infer enhancer-based gene regulatory networks (GRNs) and to use these GRNs to identify the key transcription factors in cell fate determination. You can use it to study transcription regulation during development and differentiation, or to generate a shortlist of transcription factors for trans-differentiation experiments. It is written in Python and it contains command-line script that includes the following three tools:

| Command           | Function                                                       |
| ----------------- | -------------------------------------------------------------- |
|  ananse binding   | predict cell type-specific transcription factor binding profiles   |
|  ananse network   | infer a cell type-specific gene regulatory network           |
|  ananse influence | infer key transcription factors during cell fate determination |   
|  ananse plot      | plot influence results in a dotplot and optionally a GRN of the Top TFs |   
|  ananse view      | inspect the output of `ananse binding` |


All functionality is also available through a Python API.

ANANSE is free and open source research software. If you find it useful please cite our preprint:

!!! note "Citation"
    > Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Gert Jan C Veenstra, Huiqing Zhou, Simon J van Heeringen, **ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination**, Nucleic Acids Research, Volume 49, Issue 14, 20 August 2021, Pages 7966–7985, [https://doi.org/10.1093/nar/gkab598](https://doi.org/10.1093/nar/gkab598)

### Getting started

* Install ANANSE on Linux or Mac, see the [Installation](installation.md) page for details.
* Have a look at these simple [examples](examples.md) to get a taste of what is possible.
* Check out the [command-line reference](command-line_reference.md) to get going with your own data.
* There’s also an [python API documentation](API_documentation.md) for Python users.

### Get help

* First, check the [FAQ](faq.md) for common issues.
* The preferred way to get support is through the [Github issues page](https://github.com/vanheeringen-lab/ANANSE/issues).
* Finally, you can reach us by email to <a href="mailto:simon.vanheeringen@gmail.com" target="_blank">Simon J. van Heeringen</a> or <a href="mailto:qxuchn@gmail.com" target="_blank">Quan Xu</a>.

### Full contents

* [Model description](model_description.md)
    - [Overview of ANANSE](model_description/#overview_of_ANANSE)
    - [Prediction of transcription factor binding](model_description/#prediction_of_transcription_factor_binding)
    - [Inference of gene regulatory networks](model_description/#inference_of_gene_regulatory_networks)
    - [Calculation of influence score](model_description/#calculation_of_influence_score)
* [Installation](installation.md)
    - [The easiest way to install](installation/#the-easiest-way-to-install)
    - [Alternative installation](installation/#alternative-installation)
* [Input data](input_data.md)
    - [Genome](input_data/#genome)
    - [Motif database](input_data/#motif-database)
    - [Enhancer regions](input_data/#enhancer-regions)
    - [Enhancer activity](input_data/#enhancer-activity)
    - [Expression data](input_data/#expression-data)
    - [Differential expression data](input_data/#differential-expression-data)
* [Command-line reference](command-line_reference.md)
    - [Command: ananse binding](command-line_reference/#ananse-binding)
    - [Command: ananse network](command-line_reference/#ananse-network)
    - [Command: ananse influence](command-line_reference/#ananse-influence)
    - [Command: ananse plot](command-line_reference/#ananse-plot)
    - [Command: ananse view](command-line_reference/#ananse-view)
* [API documentation](API_documentation.md)
    - [Working with the `binding` Class](API_documentation/#working-with-binding-class)
    - [Working with the `network` Class](API_documentation/#working-with-network-class)
    - [Working with the `influence` Class](API_documentation/#working-with-influence-class)
* [Examples](examples.md)
    - [Ananse setup](examples/#prepare-code-and-dataset)
    - [Example of ananse binding](examples/#build-tf-binding-network)
    - [Example of ananse network](examples/#built-gene-regulatory-network)
    - [Example of ananse influence](examples/#infer-tf-influence-score)
* [FAQ](faq.md)
    - [I can apply ANANSE to what kind of enhancer data](faq/#I_can_apply_ANANSE_to_what_kind_of_enhancer_data)
    - [I can apply ANANSE to what bioogy process](faq/#I_can_apply_ANANSE_to_what_bioogy_process)
    - [I can apply ANANSE to which species data](faq/#I_can_apply_ANANSE_to_which_species_data)
    - [Do I need provide motif data](faq/#Do_I_need_provide_motif_data)
* [ANANSE news](ananse_news.md)
    - [ANANSE v0.3.0](ananse_news/#ananse-v030)
    - [ANANSE v0.1.5](ananse_news/#ananse-v015)
    - [ANANSE v0.1.4](ananse_news/#ananse-v014)
    - [ANANSE v0.1.3](ananse_news/#ananse-v013)
    - [ANANSE v0.1.2](ananse_news/#ananse-v012)
    - [ANANSE 0.1.1](ananse_news/#ananse-v011)
    - [ANANSE 0.1.0](ananse_news/#ananse-v010)
* [Acknowledgments](acknowledgments.md)
