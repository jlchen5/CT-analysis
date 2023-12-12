# CT analysis

- Circuit Topology approach in Cancer research and Immune cells

- TCGA data of RNA-seq, WGS, scRNA-seq and Hi-C, scHi-c, etc.

#### - Hi-C data processing

1. 双端测序的单细胞Hi-C处理package：[NucProcess](https://github.com/TheLaueLab/nuc_processing)。
2. 这个软件 ([NucProcess](https://github.com/TheLaueLab/nuc_processing)) 需要原始的 *FASTQ* 测序数据, 参考基因组和实验过程 (限制性内切酶的种类和DNA片段大小的范围) . 步骤如下:
   - 查找限制性内切酶位点（restriction enzyme cut sites）
   - 根据酶切位点分割数据
   - 比对到参考基因组（建议使用Bowtie2），同时不允许错配（allowed for no base pair mismatches）
   - 比对结果分为两部分：1）双端测序结果都是唯一匹配；2）其中一条或者两条比对结果含糊不清
   - The mapped read pairs were allocated to their originating chromosomal restriction enzyme fragment regions (i.e. one for each end) to identify the specific RE1-RE1 ligation junction involved. They were then filtered to remove pairs that represent obviously aberrant molecular events or contain no useful spatial information. The filtering stages are:
     - Removal of pairs that map to the same RE1 fragment: inward facing reads carry no ligation junction, useful or otherwise, and outward facing reads represent circularisation.
     - Removal of pairs that map to sequentially adjacent, inward facing RE1 fragments: these represent re-ligation and carry no useful structural information. 
     - Removal of pairs that map to RE1 fragments which are so close that they could be observed even if there is no ligation junction between them. 
     - Removal or pairs where the combined sequence separation from the read positions to the estimated ligation junction does not match the size of the fragments in the sequenced DNA library: i.e. the predicted length is too short and the sequenced molecule must contain unknown intermediate sequence, or the predicted length is too long and the sequenced molecule has somehow been truncated.
   - The redundancy in the restriction fragment pairs, where a particular contact is observed more than once due to sequence amplification, is removed. At the same time the data is filtered so that at least two separate, albeit often identical, molecules must be paired-end sequenced to confirm a ligation junction. Accordingly, the data is output in two separate groups: 1) redundant/supported pairs that are observed more than once, which are used in the structure calculations, and 2) unique/unsupported pairs that are only observed once, and which are not used in structure calculation. The latter have been shown to represent proportionately more noise. 
   - 由于研究的 G1 单倍体细胞中每条染色体只有一个拷贝，因此任何杂乱的配对（即限制性片段末端参与了一次以上的接触）都会被删除（这样删除的contacts非常少）。
   - 输出结果为`ncc`格式，`bin_size`默认值为5 Mb。
