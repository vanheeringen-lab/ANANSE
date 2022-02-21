## FAQ

### What kind of enhancer data can I use as input for ANANSE?

In training the regression model for ANANSE, we used H3K27ac ChIP-seq and/or p300 ChIP-seq peaks as measure of enhancer activity. You can use either one, or both! If you use H3K27ac ChIP-seq data, make sure that you use another source of information to determine enhancer locations! The peaks from H3K27ac are not precise enough to yield informative regions for motif analysis! The detailed description on how to generate enhancer data as input can be found in the section [Enhancer data](https://anansepy.readthedocs.io/en/master/input_data/#enhancer-activity).

### In what kind of experiments can I apply ANANSE?

ANANSE is a method to predict key transcription factors between enhancer networks (e.g. cell fate determination). It is suitable in any type of experimemt where you expect transcription factors to be a main component of the transcriptional regulation. This includes development, (*in vitro*) differentiation, trans-differentiation, before and after treatment and so on. Please let us know if it worked (or if it didn't work) for your type of data.

### To which species can I apply ANANSE?

The model has been trained and evaluated on human data. However, the whole approach is, by design, agnostic to the species. As long as you have a genome, and set of gene annotation that maps to gene symbols you can run it. If you have non-human/mouse gene symbols or gene names, you will need a motif database for your species of interest mapping motifs to gene symbols (also see the question below).

### Do I need provide transcription factor motif data?

By default ANANSE uses a non-redundant, clustered database of known vertebrate motifs: `gimme.vertebrate.v5.0`. These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources. Large-scale benchmarks using ChIP-seq peaks show that this database shows good performance and should be a good default choice. This motif database should be fine for human or mouse data. Currently, you have to provide the motif database for other species. The detailed description of the required format of the motif database can be found in this section: [Motif database](https://anansepy.readthedocs.io/en/master/input_data/#motif-database).
