## Model Description

### Overview of ANANSE
An overview of the workflow for the ANANSE method to generate key TFs for cell conversion.

![](img/Fig2.jpg)

* (A), Data types required and utilized in ANANSE. These data include motif score of all TFs, gene expression data (e.g. RNA-seq) and enhancer data that can be obtained by ATAC-seq, EP300 ChIP-seq, or H3K27ac ChIP-seq from each cell type. The blue and orange peaks represent enhancers in two cell types. The four sequence-logos represent the motif of four TFs. The heatmap represents gene expression intensity in two cell types. 
* (B), The TF binding profiles predicted from enhancer data and TF motif scores in each cell type. Two GRNs below show cell type-specific TF binding profiles in two cell types (source and target cell types). 
* (C), The cell type-specific GRN predicted based on TF-Gene binding and TF/Gene expression. Two networks show cell type-specific GRN in two cell types. The orange circle represents a TF or a gene, and the size of the circle indicates the target gene number of the corresponding TF. The blue arrow indicates regulation between two TFs, and the color intensity represents regulation intensity. 
* (D), The differential GRN between the two cell types. In this step, the interaction specific for the target cell type is kept constant, and if the interaction score of the target cell type is higher than that of the source cell type, the interaction score is further used. 
* (E), The barplot shows the ranked influence score of all TFs calculated from the differential GRN. The influence score is calculated based on gene expression score, distance from the enhancer bound by TF to gene, and the interaction score between TF and gene.
