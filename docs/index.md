## **ANANSE**: **AN**alysis **A**lgorithm for **N**etworks **S**pecified by **E**nhancers

### What is ANANSE?

ANANSE is a computational approach to infer enhancer-based gene regulatory networks (GRNs) and to identify key transcription factors between two GRNs. You can use it to study transcription regulation during development and differentiation, or to generate a shortlist of transcription factors for trans-differentiation experiments. 

ANANSE is written in Python and comes with a command-line interface that includes the following commands:

| Command           | Function                                                       |
| ----------------- | -------------------------------------------------------------- |
|  ananse binding   | predict transcription factor binding profiles                  |
|  ananse network   | infer a gene regulatory network                                |
|  ananse influence | infer key transcription factors between two networks           |   
|  ananse plot      | plot influence results in a dotplot and optionally a GRN of the Top TFs |   
|  ananse view      | inspect the output of `ananse binding`                         |

All functionality is also available through a [Python API](API_documentation.md).

ANANSE is free and open source research software. If you find it useful please cite it:

!!! note "Citation"
    > Quan Xu, Georgios Georgiou, Siebren Frölich, Maarten van der Sande, Gert Jan C Veenstra, Huiqing Zhou, Simon J van Heeringen, **ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination**, Nucleic Acids Research, Volume 49, Issue 14, 20 August 2021, Pages 7966–7985, [https://doi.org/10.1093/nar/gkab598](https://doi.org/10.1093/nar/gkab598)
    
ANANSE can also be useful when using CAGE-seq data. If you used this tool with CAGE-seq data, please cite:

!!! note "Citation"
    > Heuts BMH, Arza-Apalategi S, Frölich S, Bergevoet SM, van den Oever SN, van Heeringen SJ, et al. **Identification of transcription factors dictating blood cell development using a bidirectional transcription network-based computational framework**. Scientific Reports 2022 12:1 [Internet]. 2022 Nov 4 [cited 2022 Dec 6];12(1):1–12. Available from: https://www.nature.com/articles/s41598-022-21148-w

### Getting started

* Install ANANSE on Linux or Mac, see the [Installation](installation.md) page for details.
* Have a look at these simple [examples](examples.md) to get a taste of what is possible.
* Check out the [command-line reference](command-line_reference.md) to get going with your own data.

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
    - [CAGE-seq](input_data/#CAGE-seq)
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
* [Acknowledgments](acknowledgments.md)
