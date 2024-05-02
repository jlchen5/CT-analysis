# CT analysis

- Circuit Topology approach in Cancer research and Immune cells

- TCGA data of RNA-seq, WGS, scRNA-seq and Hi-C, scHi-c, etc.

## Steps

- We extract the circuit topology matrix by [Extract_CT_scHiC.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/Extract_CT_scHiC.ipynb) or use [Extract_topology.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/Extract_topology.ipynb) from `.pdb` files, and combined the gene abundance with [bin_RNAseq.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/bin_RNAseq.ipynb) and [CoupleRNAseq_CT.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/CoupleRNAseq_CT.ipynb) scripts alsa gene annotation and gene expression data (recommend the raw counts).

- We analyze the single-cell RNA-seq data by [scRNAseq_analysis.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNAseq_analysis.ipynb), which can combine the multiple data set annotate the cell types. Then we use the [scRNA_TF_acti_infer.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNA_TF_acti_infer.ipynb) to infer the transcription factors activity and [scRNAseq_pthway_infer.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/scRNAseq_pthway_infer.ipynb) to deduce pathway activity.

- We analyze the single cell persudotime via the scRNA-seq data.
 
- The scRNA-seq data was analyzed with scHi-C data.

- We seperate the hic contact files to 2 files which `_1.bed` and `_2.bed`.The motif analysis is performed by [run_homer.sh](https://github.com/jlchen5/CT-analysis/blob/main/run_homer.sh).
 
     ```bash
     bash run_homer.sh
     ```

- We use [CTforCumulativeAnalysis](https://github.com/jlchen5/CT-analysis/blob/main/CTforCumulativeAnalysis.py) to generate the csv files and as [Cumulative_analysis.ipynb](https://github.com/jlchen5/CT-analysis/blob/main/Cumulative_analysis.ipynb) input files.


 

