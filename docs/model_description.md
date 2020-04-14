## Model description

### Overview of ANANSE

An overview of the workflow that ANANSE uses to prioritize transcription factors for cellular fate changes. ANANSE infers gene regulatory networks (GRNs) to predict key transcription factors (TFs) that drive cellular transitions.

![](img/Fig2.jpg)

* **(A)** An overview of the data types required by ANANSE. These data include motif score of all TFs and gene expression data (e.g. RNA-seq) and enhancer activity data for each cell type. The enhancer data can be obtained by ATAC-seq, EP300 ChIP-seq or H3K27ac ChIP-seq. The blue and orange peaks represent enhancers in two different cell types. The four sequence logos represent the motif of four TFs. The heatmap represents gene expression intensity in the two cell types. 
* **(B)** The TF binding profiles predicted from the enhancer data and TF motif scores in each cell type. Two GRNs below show cell type-specific TF binding profiles in two cell types (source and target cell types). 
* **(C)** The cell type-specific GRN predicted based on TF-gene binding, TF expression and target gene expression. The two networks show cell type-specific GRNs in two cell types. The orange circle represents a TF or a gene, and the size of the circle indicates the target gene number of the corresponding TF. The blue arrow indicates regulation between two TFs, and the color intensity represents regulation intensity as the edge weight.
* **(D)** The differential GRN between the two cell types. In this step, the interaction specific for the target cell type is kept constant, and if the interaction score of the target cell type is higher than that of the source cell type, the interaction score is further used. 
* **(E)** The barplot shows the ranked influence score of all TFs calculated from the differential GRN. The influence score is calculated based on gene expression score, distance from the enhancer bound by TF to gene, and the interaction score between TF and gene.


#### Prediction of transcription factor binding

The enhancer intensity is combined with sequence features in enhancer peaks to infer cell type-specific TF binding profiles. The enhancer locations can be obtained from ChIP-seq analyses of the transcriptional co-activator EP300, chromatin accessibility data such as ATAC-seq or a combination of TF ChIP peaks. Basically any type that gives sharp peaks would be usable. To determine the enhancer activity, we recommend to use either EP300 ChIP-seq or H3K27ac ChIP-seq, as both of these have been shown to correlate with enhancer activity (Creyghton et al., 2010; Rada-Iglesias et al., 2011). The enhancer activity is combined with TF motif scores in the enhancer sequences usin logistic regression. The motif analysis is performed using [GimmeMotifs](https://gimmemotifs.readthedocs.org).

#### Inference of gene regulatory networks

ANANSE infers cell type-specific GRNs based on the predicted TF binding sites and the expression levels both TFs as well as their target genes. TF-gene interaction scores, the edge weights in the network, are calculated based on the predicted TF binding probability, the distance between the enhancer and the target gene, and expression of both TF and the target gene. By integrating these data, ANANSE determines the interaction score of each TF-gene pair. 

The weighted sum of TF predicted enhancer intensity within 100kb around TSS is defined as the TF-gene binding score (Eq. 1). 

<!-- \begin{equation*} -->
$$ B_{x,r} = \sum_{k} w_k s_k $$
<!-- \end{equation*} -->

where $ B_{x,r} $ is the binding score between TF $x$ and target gene $r$, $w_k$ is the weighted distance between an enhancer and the target gene and where $s_k$ is predicted binding intensity at genomic position $k$ of TF $x$. 

The distance weight is based on a linear genomic distance between the enhancer and the TSS of a gene according to: 

\begin{equation}
w_k =
  \begin{quad}
    0               & \quad k \in (\text{0kb,2kb}]\\
    1               & \quad k \in (\text{2kb,5kb}]\\
    \frac{2e^{-\mu|k-t_r|}}{1+e^{-\mu|k-t_r|}}  & \quad k \in (\text{5kb,100kb}]
  \end{quad}
\end{equation}

where $t_r$ is the genomic position of the TSS of gene $r$ and the parameter $\mu$, which determines the decay rate as a function of distance from the TSS, is set such that an enhancer 10 kb from the TSS contributes one-half of that at the TSS. This distance weight calculation is similar to the method previously described in Wang et al., 2016, except that only signal in enhancers is used, enhancers within 2kb around TSS are removed and the weight of enhancers within 2kb to 5kb is set to 1.

We scaled the expression level of the TF and the target gene, expressed as transcripts per million (TPM), and the TF-gene binding score $B_{x,r}$  we calculated in the first step from 0 to 1, with 1 being the highest and 0 the lowest. Combining the TF-gene binding score and TF and target expression scores by taking the mean, we obtained a TF-gene interaction score.

#### Calculation of influence score

**TODO**