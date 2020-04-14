## FAQ

### I can apply ANANSE to what kind of enhancer data
In our tranning model, we used EP300 ChIP-seq peak as enhancer. In our examples, we also test to use ATAC-seq together with H3K27ac ChIP-seq as enhancer. The detail discription how to generate enhancer data as input can be found at [Enhancer data](https://anansepy.readthedocs.io/en/latest/input_data/#enhancer-data) part.

### I can apply ANANSE to what bioogy process
ANANSE is a method to predict key transcription factors in cell fate determination using enhancer networks. We firmly believe it can be used in any cell fate determination data, which include cell proliferation, differentiation, cellular movement, trans-differentiation and so on.

### I can apply ANANSE to which species data
All the tranning and test data are from human, but we firmly believe you could apply ANANSE to other species.

### Do I need provide motif data
By default ANANSE uses a non-redundant, clustered database of known vertebrate motifs: `gimme.vertebrate.v5.0`. These motifs come from CIS-BP (http://cisbp.ccbr.utoronto.ca/) and other sources. Large-scale benchmarks using ChIP-seq peaks show that this database shows good performance and should be a good default choice. The detail discription about motif database can be found at [Motif database](https://anansepy.readthedocs.io/en/latest/input_data/#motif-database) part.

