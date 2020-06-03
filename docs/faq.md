## FAQ

### What kind of enhancer data can I use as input for ANANSE?

In training the regression model for ANANSE we used EP300 ChIP-seq peaks as measure of enhancer activity. In our examples, we also used ANANSE on H3K27ac ChIP-seq signal in ATAC-seq peaks. You probably can also use ATAC-seq signal directly. This will have a weaker correlation with enhancer activity, but it will likely work as well. This is something that we are currently testing. If you use H3K27ac ChIP-seq data, make sure that you use another source of information to determine enhancer locations! The peaks from H3K27ac are not precise enough to yield informative regions for motif analysis! The detailed description on how to generate enhancer data as input can be found in the section [Enhancer data](https://anansepy.readthedocs.io/en/latest/input_data/#enhancer-data).

### In what kind of cell types or experiments can I apply ANANSE?

ANANSE is a method to predict key transcription factors in cell fate determination using enhancer networks. It is suitable in any type of experimemt where you expect transcription factors to be a main component of the transcriptional regulation. This includes development, (*in vitro*) differentiation, trans-differentiation, before and after treatment and so on. Please let us know if it worked (or if it didn't work) for your type of data.

### To which species can Iapply ANANSE?

The model has been trained and evaluated on human data. However, the whole approach is, by design, agnostic to the species. As long as you have a genome, and set of gene annotation that maps to human gene symbols you can run it straightaway. If you have different gene symbols or gene names, you would need a motif database for your species of interest mapping motifs to gene symbols (also see the question below).

### Do I need provide transcription factor motif data?

By default ANANSE uses a non-redundant, clustered database of known vertebrate motifs: `gimme.vertebrate.v5.0`. These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources. Large-scale benchmarks using ChIP-seq peaks show that this database shows good performance and should be a good default choice. This motif database should be fine for human or mouse data. Currently, you have to provide the motif database for other species. The detailed description of the required format of the motif database can be found in this section: [Motif database](https://anansepy.readthedocs.io/en/latest/input_data/#motif-database).

