# CT analysis

- Circuit Topology approach in Cancer research and Immune cells

- TCGA data of RNA-seq, WGS, scRNA-seq and Hi-C, scHi-c, etc.

## Steps

1 We extract the topology via [Extract_CT_scHiC.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/Extract_CT_scHiC.ipynb), and combined the gene abandency with [bin_RNAseq.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/bin_RNAseq.ipynb) and [CoupleRNAseq_CT.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/CoupleRNAseq_CT.ipynb) scripts alsa gene annotation and gene expression data (recommend the raw counts).

2 We analyze the single-cell RNA-seq data by [scRNAseq_analysis.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNAseq_analysis.ipynb), which can combine the multiple data set annotate the cell types. Then we use the [scRNA_TF_acti_infer.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNA_TF_acti_infer.ipynb) to infer the transcription factors activity and [scRNAseq_pthway_infer.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNAseq_pthway_infer.ipynb) to deduce pathway activity.



