# 😀 🌁 2023.12.4

## 一、principle of Analysis  

### Hi-C data processing

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



### 3D genome structure calculations and analysis 

3D genome structures were calculated using simulated annealing, of an initially random conformation of a particle-on-a-string representation of the chromosomes, to generate structures that were compatible with the experimental distance restraints derived from the Hi-C contacts. 

The software written to perform this task (which includes both a command line version and Jupyter notebook) is available at: https://github.com/TheLaueLab/nuc_dynamics. 



## 二、代码部分

~~~shell
# 文章中的单细胞Hi-C的ncc格式文件利用nuc_tools可视化为contact map
(base) jialechen@Jiale-MacBook-Pro ~/Desktop/PhD/Pang_2022_GenomeBiol_3D/ncc $ nuc_tools contact_map  GSM6081054_Cell_ID_1.ncc -o Cell_ID_1.svg -sc # -sc表示为single cell

	INFO: Making PDF contact map for GSM6081054_Cell_ID_1.ncc
	INFO: Loading NCC format contact data
	INFO: Reading GSM6081054_Cell_ID_1.ncc
	INFO:  .. found 16,238 contacts
	
	INFO: Min. contig size not specified, using 10.0% of largest: 24,814,690 bp
	INFO: Considering 23 chromosomes/contigs
	INFO: Skipped 1 small chromosomes/contigs < 24,814,690 bp
	INFO: Binning contact lists
	INFO:  .. chrX - chrX : 20
	INFO: Full contact map size 1488 x 1488
	INFO:  .. making map Chromosome chrX (dpi=380)
	INFO: Written Cell_ID_1.pdf
~~~

<img src="/Users/jialechen/Desktop/PhD/Procedure_Circuit_Topology.assets/截屏2023-12-04 16.56.39.png" alt="截屏2023-12-04 16.56.39" style="zoom:33%;" />





# 😉 🌧️ 2023.12.5

## 一、代码部分 （文件格式转换）

~~~shell
# 将单细胞的ncc格式文件可视化代码
$ nuc_tools contact_map  GSM6081054_Cell_ID_1.ncc -o Cell_ID_1.svg -sc

# sample name 
Jiale-MacBook-Pro:ncc $ cat sample.txt
GSM6081054_Cell_ID_1
GSM6081055_Cell_ID_2
GSM6081056_Cell_ID_3
GSM6081057_Cell_ID_4
GSM6081058_Cell_ID_5
GSM6081059_Cell_ID_6
GSM6081060_Cell_ID_7
GSM6081061_Cell_ID_8
GSM6081062_Cell_ID_9
GSM6081063_Cell_ID_10
GSM6081064_Cell_ID_11
GSM6081065_Cell_ID_12
GSM6081066_Cell_ID_13
GSM6081067_Cell_ID_14
GSM6081068_Cell_ID_15
GSM6081069_Cell_ID_16
GSM6081070_Cell_ID_17
GSM6081071_Cell_ID_18
GSM6081072_Cell_ID_19
GSM6081073_Cell_ID_20
GSM6081074_Cell_ID_21
GSM6081075_Cell_ID_22
GSM6081076_Cell_ID_23
GSM6081077_Cell_ID_24
GSM6081078_Cell_ID_43
GSM6081079_Cell_ID_44
GSM6081080_Cell_ID_49
GSM6081081_Cell_ID_56
GSM6081082_Cell_ID_84
GSM6081083_Cell_ID_86
GSM6081084_Cell_ID_89
GSM6081085_Cell_ID_90
GSM6081086_Cell_ID_93

# 提取ncc格式文件为contact txt格式
$ cat sample.txt|while read line ;do (nohup awk  '{print $1"\t"$2"\t"$7"\t"$8}' $line.ncc > $line.contact.txt &);done

-rw-r--r--@ 1 jialechen  staff   207K 12  5 11:13 GSM6081068_Cell_ID_15.contact.txt
-rw-r--r--@ 1 jialechen  staff   465K 12  5 11:13 GSM6081054_Cell_ID_1.contact.txt
-rw-r--r--@ 1 jialechen  staff   512K 12  5 11:13 GSM6081058_Cell_ID_5.contact.txt
-rw-r--r--@ 1 jialechen  staff   384K 12  5 11:13 GSM6081070_Cell_ID_17.contact.txt
-rw-r--r--@ 1 jialechen  staff   652K 12  5 11:13 GSM6081067_Cell_ID_14.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.3M 12  5 11:13 GSM6081061_Cell_ID_8.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.3M 12  5 11:13 GSM6081069_Cell_ID_16.contact.txt
-rw-r--r--@ 1 jialechen  staff   833K 12  5 11:13 GSM6081073_Cell_ID_20.contact.txt
-rw-r--r--@ 1 jialechen  staff   861K 12  5 11:13 GSM6081084_Cell_ID_89.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.7M 12  5 11:13 GSM6081056_Cell_ID_3.contact.txt
-rw-r--r--@ 1 jialechen  staff   872K 12  5 11:13 GSM6081081_Cell_ID_56.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.0M 12  5 11:13 GSM6081078_Cell_ID_43.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.0M 12  5 11:13 GSM6081083_Cell_ID_86.contact.txt
-rw-r--r--@ 1 jialechen  staff   803K 12  5 11:13 GSM6081086_Cell_ID_93.contact.txt
-rw-r--r--@ 1 jialechen  staff   783K 12  5 11:13 GSM6081085_Cell_ID_90.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.1M 12  5 11:13 GSM6081059_Cell_ID_6.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.4M 12  5 11:13 GSM6081079_Cell_ID_44.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.2M 12  5 11:13 GSM6081075_Cell_ID_22.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.4M 12  5 11:13 GSM6081080_Cell_ID_49.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.5M 12  5 11:13 GSM6081064_Cell_ID_11.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.3M 12  5 11:13 GSM6081060_Cell_ID_7.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.7M 12  5 11:13 GSM6081066_Cell_ID_13.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.7M 12  5 11:13 GSM6081082_Cell_ID_84.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.7M 12  5 11:13 GSM6081062_Cell_ID_9.contact.txt
-rw-r--r--@ 1 jialechen  staff   1.9M 12  5 11:13 GSM6081077_Cell_ID_24.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.5M 12  5 11:13 GSM6081063_Cell_ID_10.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.9M 12  5 11:13 GSM6081074_Cell_ID_21.contact.txt
-rw-r--r--@ 1 jialechen  staff   2.7M 12  5 11:13 GSM6081071_Cell_ID_18.contact.txt
-rw-r--r--@ 1 jialechen  staff   3.3M 12  5 11:13 GSM6081072_Cell_ID_19.contact.txt
-rw-r--r--@ 1 jialechen  staff   4.1M 12  5 11:13 GSM6081057_Cell_ID_4.contact.txt
-rw-r--r--@ 1 jialechen  staff   4.3M 12  5 11:13 GSM6081055_Cell_ID_2.contact.txt
-rw-r--r--@ 1 jialechen  staff   3.9M 12  5 11:13 GSM6081065_Cell_ID_12.contact.txt
-rw-r--r--@ 1 jialechen  staff   3.9M 12  5 11:13 GSM6081076_Cell_ID_23.contact.txt


# 给文件添加首行 / chr_A	pos_A	chr_B	pos_B

$ cat sample.txt|while read line;do echo -e "chr_A\tpos_A\tchr_B\tpos_B"|cat - $line.contact.txt > $line.contact_pair.txt && mv $line.contact_pair.txt $line.contact.txt ;done

# 按照单细胞在/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/results/matrices路径下新建文件夹
$ cat sample.txt|sed 's/_/\t/' |cut -f 2 |while read line ;do mkdir -p ../results/matrices/$line ;done

~~~



# 🥱 ☁️ 2023.12.6

## 一、代码分析

- 使用ExtractTopology_singlecellsHiC.ipynb代码从单细胞数据中提取CT结果，地址为：http://localhost:8888/notebooks/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/Hi-C%20map%20analysis/ExtractTopology_singlecellsHiC.ipynb。



# 😐 🌁 2022.12.7

## 一、代码部分

### 1. 使用nuc_dynamics将ncc文件生成基因组结构pdb文件：

~~~shell
# Running NucDynamics
# Typical use, generating 10 conformational models in PDB format:
$ nuc_dynamics example_chromo_data/Cell_1_contacts.ncc -m 10 -f pdb

# Specifying the particle sizes (8 Mb, 2 Mb, 1 Mb, 500 kb) and an output file
name:
$ nuc_dynamics example_chromo_data/Cell_1_contacts.ncc -m 10 -f pdb -o Cell_1.pdb -s 8 2 1 0.5
~~~



但是！！文章给的ncc文件格式太老，现有的工具不支持，无语。还是继续看Barbara的代码。



# 😄 ☁️ 2023.12.11

## 一、condensate project (preparation)

1. phase-separated biomolecular condesates
2. membraneless compartments (made of IDPs)
3. CellHesion device (atomic force microscope)
4. oppositely charged macromolecules: PLL:H and PLL:HS
5. laminin and Fus protein use various crowding agents and also various buffers



## 二、analysis data

分析卵巢癌（ovarian cancer cell）Hi-C数据中的： - P*「E」* / - S*「U」*/ - X*「E」*



# 🥱 ☁️ 2023.12.12



- 把hic contact txt提取Circuit topology的代码（/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/Hi-C map analysis/ExtractTopology_singlecellsHiC.ipynb）中加入提取染色体信息的PSX。



# 🥱 ☁️ 2023.12.13

- **laminin** structure and droplet formation 
- Laminins are large cell-adhesive glycoproteins
- Calculate the net charge





#  😉 ☁️ 2023.12.15

- Generate the hic.txt from ncc file 

~~~shell
$ cat sample.txt|while read line ;do (nohup awk  '{print $1"\t"$2"\t"$3"\t"$7"\t"$8"\t"$9"\t"$15}' $line.ncc > $line.hic.txt &);done

$ cat sample.txt|while read line;do echo -e "#chr_A\tstart_A\tend_A\tchr_B\tstart_B\tend_b\tnum_obs"|cat - $line.hic.txt > $line.hic_hder.txt  ;done

$ rm  *hic.txt

$ cat sample.txt|while read line;do mv $line.hic_hder.txt $line.hic.txt;done 
~~~

生成的结果不尽人意，RNAseq数据是tpm，是标准化的，但是源代码是用的counts值，所以出来的结果大部分是Nan。



# 😃 ☀️2024.02.12 大年初三

办公室门锁坏了，只能在旁边的会议室办公。

修复了原来coupleRNAseq_CT.ipynb的bug，只是加了the codes below：

~~~shell
# Define the directory path
output_dir = 'results/correlation rna seq/{}/'.format(resolution_string)

# Create the directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
~~~



# 😄 ☀️ 2024.02.13

- Pang_2022_GenomeBiol_3D
  - Cell ID-1~24: single-cell Hi-C profiling of ovarian cancer cell line PEO1 and HEYA8
  - Cell ID-43/44/49/56/84/86/89/90/93: single-cell Hi-C profiling of OVCA429 shGRHL2 Tet-inducible, an ovarian cancer cell line 
- `contact map` transfer to `sparse matrix`

~~~python
## install dependent packages
!pip install numpy scipy scikit-image tqdm matplotlib

## 使用skimage库中的io模块打开tiff文件，并将其转换为科学计算库scipy中的稀疏矩阵格式csr_matrix	
from tqdm import tqdm
from skimage import io
from scipy.sparse import csr_matrix

# set path	
filename = "path/to/large_file.tiff"

# open and read tif 
image = io.imread(filename)

# 将tif文件矩阵转换为稀疏矩阵
sparse_matrix = csr_matrix(image)

## 如有必要，可以选择进行降采样以减少稀疏矩阵的大小。可以使用skimage库中的transform模块进行降采样：
from skimage import transform

# 设置降采样倍数
downsample = 5

# 对稀疏矩阵进行降采样
downsampled_matrix = transform.downscale_local_mean(sparse_matrix.toarray(), (downsample, downsample))
~~~



# 😅 ☔️ 2024.02.14 

- Use `Local_matrix_analysis.ipynb` script to analyze the chromosome segments, and get the `.csv` file.
- This code can extract the `.tif` matrix and analyze the L-loops.



# 2024.02.15

开完start up meeting，老板讲了很多物理知识，包括固体/流体/牛顿流体/viscoelastic等。

安排我尽快找到带有HBD的laminin的蛋白（带荧光）然后买，heparin binding domain位于laminin的alpha chain LG4 domains。



# 2024.02.16

- 准备做周一的ppt，汇报一下结果：首先是背景：为什么要分析宫颈癌单细胞Hi-C的数据，结果（hic contact和ct_map）。



- 循环读取hic.txt文件代码：

~~~python
import os

data_dir = '/Users/jialechen/Desktop/PhD/CT/Pang_2022_GenomeBiol_3D/data/'
file_names = [f for f in os.listdir(data_dir) if f.startswith('GSM') and f.endswith('.hic.txt')]


for file in file_names:
    contacts = pd.read_csv(os.path.join(data_dir, file), sep='\t')
~~~



# 2024.02.19

单细胞数据分析HGSOC，推断通路，拟时序分析

- 牛顿流体
- Hook solid



# 2024.02.20

- create the STAR genome index of repeats

  ~~~shell
  $ STAR --runThreadN 10 --runMode genomeGenerate \
   			 --genomeDir  ./ #folder of index results \
         --genomeFastaFiles ~/Desktop/Analysis_work/reference/Hsapiens/hg38.fa  \
         --sjdbGTFfile ~/Desktop/Analysis_work/reference/Hsapiens/repeat.hg38.gtf \
         --limitSjdbInsertNsj 5005834
  ~~~
  



# 2024.02.21 learn the concentration

learn the concentration with Dirk-jan and calcalate the different units about gram, Liter, etc.

| 英文名称     | 缩写  | 中文名称     | 备注                                                         |
| ------------ | ----- | ------------ | ------------------------------------------------------------ |
| gram         | g.    | 克           | 适应于表示重量的情况，如：this chip is 10g by weight.        |
| centigram    | cg    | 厘克         | 重量单位                                                     |
| milligram    | mmg   | 毫克         | 重量单位                                                     |
| kilogram     | kg.   | 千克（公斤） | 通常，kilo这个词根表示“千”，如kilogram（千克），kilometer(千米)。此处，也是一个重量单位。 |
| metric ton   | m.t.  | 公吨         | 重量单位                                                     |
| long ton     | l.t.  | 长吨         | 重量单位                                                     |
| short ton    | sh.t. | 短吨         | 重量单位                                                     |
| pound        | lb.   | 磅           | 重量单位1lb.=0.454kg                                         |
| ounce        | oz    | 盎司         | 重量单位1 oz=31.1035g                                        |
| meter        | m     | 米           | 长度单位                                                     |
| kilometer    | km    | 千米（公里） | 长度单位                                                     |
| decimeter    |       | 分米（公寸） | 长度单位                                                     |
| centimeter   | cm    | 厘米（公分） | 长度单位                                                     |
| millimeter   | mm    | 毫米         | 长度单位                                                     |
| micrometer   | μm    | 微米         | 长度单位                                                     |
| nanometer    | nm    | 纳米         | 长度单位                                                     |
| yard         | yd    | 码           | 长度单位                                                     |
| foot         | ft    | 英尺         | 长度单位1ft=12in=30.48cm                                     |
| inch         | in.   | 英寸         | 长度单位                                                     |
| square meter | sq.m  | 平方米       | 面积单位                                                     |
| square foot  | sq.ft | 平方英尺     | 面积单位                                                     |
| square inch  | sq.in | 平方英寸     | 面积单位                                                     |
| square yard  | sq.yd | 平方码       | 面积单位                                                     |
| liter        | l.    | 升           | 容积单位                                                     |
| milliliter   | ml    | 毫升         | 容积单位                                                     |
| gallon       | gal   | 加仑         | 容积单位                                                     |
| pint         |       | 品脱         | 容积单位                                                     |
| bushel       | bu    | 蒲式耳       | 容积单位                                                     |
| cubic meter  | cu.m  | 立方米       | 体积单位                                                     |
| cubic foot   | cu.ft | 立方英尺     | 体积单位                                                     |
| piece        | pc    | 只，件，块   | 数量单位                                                     |
| package      | pkg   | 件，包       | 数量单位                                                     |
| pair         |       | 双，对       | 数量单位                                                     |
| set          |       | 台，套，架   | 数量单位                                                     |
| dozen        | doz.  | 打           | 数量单位                                                     |
| gross        | gr.   | 罗           | 数量单位                                                     |
| great gross  | g.gr. | 大罗         | 数量单位                                                     |
| ream         | rm    | 令           | 数量单位                                                     |
| roll         |       | 卷           | 数量单位                                                     |
| unit         |       | 件，辆       | 数量单位                                                     |
| head         |       | 头           | 数量单位                                                     |
| barrel       |       | 桶           | 数量单位                                                     |
| bag          |       | 袋           | 数量单位                                                     |
| lot          |       | 批           | 数量单位                                                     |
| bar          |       | 条           | 数量单位                                                     |
| batch        |       | 批           | 数量单位                                                     |
| bolt         |       | 匹           | 数量单位                                                     |
| sheet        |       | 张           | 数量单位                                                     |
| bunch        |       | 串，束       | 数量单位                                                     |
| pile         |       | 堆           | 数量单位                                                     |
| group        |       | 组，套       | 数量单位                                                     |
| portion      |       | 份           | 数量单位                                                     |
| set          |       | 套           | 数量单位                                                     |
| hour         | h     | 小时         | 时间单位                                                     |
| minute       | min   | 分钟         | 时间单位                                                     |
| second       | sec   | 秒           | 时间单位                                                     |
| ampere       | A     | 安培         | 电流单位                                                     |
| milliampere  | mA    | 毫安         | 电流单位                                                     |
| microampere  | μA    | 微安         | 电流单位                                                     |
| volt         | V     | 伏特         | 电压单位                                                     |
| millivolt    | mV    | 毫微         | 电压单位                                                     |
| microvolt    | μV    | 微伏         | 电压单位                                                     |
| Ohm          | Ω     | 欧姆         | 电阻单位                                                     |
| milliohm     | mΩ    | 毫欧         | 电阻单位                                                     |
| microohm     | μΩ    | 微欧         | 电阻单位                                                     |
| watt         | W     | 瓦特         | 功率单位                                                     |
| kilowatt     | KW    | 千瓦         | 功率单位                                                     |
| Kilocalorie  | kcal  | 千卡，大卡   | 热量单位                                                     |
| calorie      | cal   | 卡路里       | 热量单位                                                     |



# 2024.02.22

- download the reference genome:https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/ with the gene name

- create STAR genome index 

  ~~~
  # coding_genes
  STAR --runThreadN 6 --runMode genomeGenerate --genomeDir coding_genes/ --genomeFastaFiles ~/Desktop/Analysis_work/reference/Hsapiens/hg38.fa --sjdbGTFfile ~/Desktop/Analysis_work/reference/Hsapiens/hg38.refGene.gtf --sjdbOverhang 149
  
  # repeats
  STAR --runThreadN 6 --runMode genomeGenerate --genomeDir repeats/ --genomeFastaFiles ~/Desktop/Analysis_work/reference/Hsapiens/hg38.fa --sjdbGTFfile ~/Desktop/Analysis_work/reference/Hsapiens/repeat.hg38.gtf --sjdbOverhang 149
  ~~~

  



# 2024.03.05

- remote connect lab computer and I can use it now : )
- I have downloaded the raw fastq data of the single cell data from Pang 2022 Genome Biology.
- try to use the `nuc_process` to process the raw data

~~~shell
# nuc_process

# step 0: change the fastq name from _1.fastq to _r_1.fastq
# The splitFastqBarcodes.py script is provided to split FASTQ files that represent many cells, each with a different barcode sequence, into separate paired read files. The script can be run as follows, specifying the names of the two paired, multiplex FASTQ read files after the script:
 
(base) msbb@msbb-System-Product-Name:~/jiale/projects/Pang_2022_GenomeBiol_3D/data/fastq$ cat /home/msbb/jiale/projects/Pang_2022_GenomeBiol_3D/data/fastq/sample.txt  |while read line ;do split_fastq_barcodes ${line}_r_1.fastq.gz ${line}_r_2.fastq.gz -s 3 ;done

   # INFO: No -b barcode file specified.
   # INFO: No FASTQ files will be written.
   # INFO: Only analysis will be performed.


# step 1: creates a genome index and RE-digest files:
(jiale) msbb@msbb-System-Product-Name:~/jiale/projects/Pang_2022_GenomeBiol_3D/data/fastq$ nuc_process -f ~/jiale/library/chromosome/*.fa -o cell_id_1 -v -a -k -re1 MboI -re2 AluI -s 150-2000 -n 12 -g ~/jiale/tools/nuc_processing/hg38/  FASTQ_FILE ./Cell_ID-1_r_1.fastq.gz ./Cell_ID-1_r_2.fastq.gz

~~~



# 2024.03.06 🌁 & ☀️

➡️ make condensate with certain concentration: 

- HNa2O4P  100mM
- KCl              200mM
- PEG             230mg/ml
- BSA             30mg/ml

➡️ the density of the condensate is high.









# 2024.03.11

## 绘制condensate不同浓度结果图

~~~python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

# 你的数据
data = {
    'Type': ['one', 'one', 'two', 'two', 'two', 'one', 'one', 'two', 'two', 'two', 'one', 'one', 'two', 'two'],
    'BSA_Concentration': [15, 15, 15, 15, 15, 30, 30, 30, 30, 30, 60, 60, 60, 60],
    'PEG_Concentration': [50, 100, 150, 200, 230, 50, 100, 150, 200, 230, 25, 50, 100, 150]
}

df = pd.DataFrame(data)

# 分离两个相区的数据
one_phase = df[df['Type'] == 'one']
two_phase = df[df['Type'] == 'two']

# 计算每个BSA浓度级别的"one"和"two"相区PEG浓度的平均值的中点
bsa_levels = sorted(df['BSA_Concentration'].unique())
middle_points = []

for level in bsa_levels:
    one_avg = one_phase[one_phase['BSA_Concentration'] == level]['PEG_Concentration'].mean()
    two_avg = two_phase[two_phase['BSA_Concentration'] == level]['PEG_Concentration'].mean()
    middle_points.append((level, (one_avg + two_avg) / 2))

middle_points = np.array(middle_points)

# 线性插值来生成平滑的分界线
f_linear = interp1d(middle_points[:,0], middle_points[:,1], kind='linear', fill_value='extrapolate')
xnew = np.linspace(np.min(bsa_levels), np.max(bsa_levels), num=100, endpoint=True)

# 绘制原始数据点
plt.scatter(one_phase['BSA_Concentration'], one_phase['PEG_Concentration'], color='blue', label='One Phase')
plt.scatter(two_phase['BSA_Concentration'], two_phase['PEG_Concentration'], color='red', label='Two Phases')

# 绘制平滑的分界线
plt.plot(xnew, f_linear(xnew), 'k--', label='Phase Boundary')

# 填充两个相区域
plt.fill_between(xnew, f_linear(xnew), df['PEG_Concentration'].min() - 10, color='lightblue', alpha=0.2,label='One Phase')
plt.fill_between(xnew, f_linear(xnew), df['PEG_Concentration'].max() + 10, color='salmon', alpha=0.2,label='Two Phases')

plt.text(one_phase['BSA_Concentration'].mean(), one_phase['PEG_Concentration'].mean(), 'One Phase', color='blue', horizontalalignment='center', verticalalignment='bottom')
plt.text(two_phase['BSA_Concentration'].mean(), two_phase['PEG_Concentration'].mean(), 'Two Phases', color='red', horizontalalignment='center', verticalalignment='top')


plt.title('PEG-BSA Phase Diagram')
plt.xlabel('BSA  (mg/mL)')
plt.ylabel('PEG  (mg/mL)')
#plt.legend()
plt.show()

~~~



## 订购peptides and DNA （2024-Science Adv）

- 文章推荐：We used ***poly(dT) sequences*** featuring different numbers of nucleotides ranging from 20 to 200 nt (dT20, dT40, dT90, and dT200).

| System        | 𝐸 ஺(RT) | 𝐸 ஺(RT) | 𝐸 ஺(RT) | 𝐸 ஺േ ∆𝐸 ஺(RT) |
| ------------- | ------- | ------- | ------- | ------------- |
| [RPRPP]5-dT40 | 10.35   | 8.44    | 8.72    | 9±1           |
| [RGRGG]5-dT40 | 15.97   | 18.87   | 17.23   | 17 ± 1        |
| [RGYGG]3-dT40 | 27.01   | 22.07   | 21.58   | 24 ± 3        |



# 2024.3.12 ☁️

- peptides and DNA：文章中使用最多的就是[RGRGG]5-dT40，所以本实验也可以使用同上序列。
- 用nuc_process将fastq转为ncc（new）：

~~~shell
# github
$ nuc_process  /data/hi-c/Cell1.r[12].fq.gz -o Cell1 -v -a -k -re1 MboI -re2 AluI -s 50-5000 -n 8 -f /data/genome/mm10_chromosomes/*.fa -g /data/genome/mm10 -cn /data/genome/mm10_chromo_names.tsv

#my server
### Typical first-time use, which creates a genome index and RE-digest files:
$ nuc_process  -f /home/msbb/jiale/library/hg38.fa -o cell -v -a -k -re1 MboI -re2 AluI -s 150-2000 -n 12 -g /home/msbb/jiale/library/Bowtie2Index/hg38  ~/jiale/projects/Pang_2022_GenomeBiol_3D/data/fastq/Cell_ID-10_r[12].fastq

$ nuc_process Cell_ID-11_r_[12].fastq.gz  -o Cell1 -v -a -k -re1 MboI -re2 AluI -s 50-5000 -n 8  -f /home/msbb/jiale/library/hg38.fa  -g /home/msbb/jiale/library/Bowtie2Index/hg38  -cn ~/jiale/tools/nuc_processing/hg38_chromo_names.tsv
~~~



- Activate lion env: 

~~~
source venv/bin/activate
conda activate jiale 
~~~



# 2024.3.13

~~~bash
# use `nuc_sequence_names` to generate the tsv file

$ nuc_sequence_names  hg38_chromo_names.tsv j.chen@lacdr.leidenuniv.nl ~/jiale/library/hg38.fa
~~~



# 2024.3.14 ☀️

- Nuc_process运行成功！！！

```bash
(jiale) msbb@msbb-System-Product-Name:/media/msbb/J.Chen_MSBB/PhD/CT/Pang_2022_GenomeBiol_3D/data/fastq$ nuc_process Cell_ID-10_r_[12].fastq -o Cell1 -v -a -k -re1 MboI -re2 AluI -s 50-5000 -n 8 -f /home/msbb/jiale/library/hg38.fa -g /home/msbb/jiale/library/Bowtie2Index/hg38 -cn ~/jiale/tools/nuc_processing/hg38_chromo_names.tsv


```

- change a loop command now

```bash
cat sample.txt|while read line ;do (nohup nuc_process ${line}_r_[12].fastq.gz -o $line -v -a -k -re1 MboI -re2 AluI -s 50-5000 -n 8 -f /home/msbb/jiale/library/hg38.fa -g /home/msbb/jiale/library/Bowtie2Index/hg38 -cn ~/jiale/tools/nuc_processing/hg38_chromo_names.tsv &);done 

```



- Use `nuc_dynmics` to transfer data from `ncc` to `pdb`

```bash
## Typical use, generating 10 conformational models in PDB format:

$ nuc_dynamics Cell_ID-10_nuc_pair.ncc -m 10 -f pdb

## Specifying the particle sizes (8 Mb, 2 Mb, 1 Mb, 500 kb) and an output file
name:

$ nuc_dynamics example_chromo_data/Cell_1_contacts.ncc -m 10 -f pdb -o Cell_1.pdb -s 8 2 1 0.5
```





# 2024.3.20

## transfer the hic contact files to bed files (same chromosome).



- ~~~bash
  $ head GSM6081054_Cell_ID_01.contact.txt                                               
  chr_A	pos_A	chr_B	pos_B
  chr10	100160223	chr2	67288359
  chr10	100201623	chr3	127643636
  chr10	100288384	chr2	225501213
  chr10	100288397	chr2	225501213
  chr10	10031568	chr2	24049144
  chr10	100322067	chr6	90732008
  chr10	100322078	chr6	90732008
  chr10	100331209	chr12	31235800
  chr10	100427285	chr5	63196830
  
  $ head Cell_ID_01.bed                                                                
  chr1	629508	629357	chr1	634058	634209
  chr1	2045342	2045455	chr1	113465200	113465280
  chr1	2045342	2045455	chr1	113465223	113465373
  chr1	2054571	2054420	chr1	45397481	45397526
  chr1	2168591	2168678	chr1	23909210	23909059
  chr1	2796648	2796786	chr1	41629119	41629083
  chr1	3449613	3449763	chr1	187898227	187898378
  chr1	3655493	3655342	chr1	196643496	196643345
  chr1	7753440	7753289	chr1	221959932	221960081
  ~~~

- Code: 

~~~bash
$ sed '1d' GSM6081054_Cell_ID_01.hic.txt |awk '{if ($1 == $4)print $0}' |sort -k1,1 -k2,2n |cut -f 1-6 >Cell_ID_01.bed

# make a loop
$ ls *hic.txt|sed 's/.hic.txt//' |while read line;do sed '1d' $line.hic.txt |awk '{if ($1 == $4)print $0}' |sort -k1,1 -k2,2n |cut -f 1-6 >$line.bed ;done 
~~~





# 2024.3.21

## bedtools intersect

~~~bash
(base) [jialechen@Jiale-MacBook-Pro hic_contact (main ✗)]$ cut -f 1-3 Cell_ID_01.bed >Cell_ID_01_chrA.bed
(base) [jialechen@Jiale-MacBook-Pro hic_contact (main ✗)]$ bedtools intersect -a Cell_ID_01_chrA_new.bed -b /Users/jialechen/Desktop/Analysis_work/reference/Hsapiens/ENCFF420VPZ.bed
~~~





 # 2024.3.27 contact motif analysis

## find the motifs regions with hic_contact

~~~bash
# step 1: seperate the bed file to 2 files
$  cut -f 1-3  Cell_ID_12.bed   > Cell_ID_12_1.bed 
$  cut -f 4-6  Cell_ID_12.bed   > Cell_ID_12_2.bed 
# step 2: install genome 
$ sudo perl /Users/jialechen/Desktop/Analysis_work/App/homer/.//configureHomer.pl -install hg38
# step 3: use homer
$ findMotifsGenome.pl Cell_ID_12_1.bed hg38 ../motif/ -size 200 -mask

~~~



# 2024.4.2 Urea condition in BSA condensate

- Analyze the adding the Urea of different concertration to BSA condensate and explore the appropriate condition.

![0321_statis_signi](/Users/jialechen/Desktop/PhD/Condensate/0321_statis_signi.png)

# 2024.4.18

- [Building a custom data service - An example using single-cell Hi-C data](https://github.com/nucleome/Tutorial-SingeCellHiC?tab=readme-ov-file#start_demo)



- [scHiCTools](https://github.com/liu-bioinfo-lab/scHiCTools)

A computational toolbox for analyzing single cell Hi-C (high-throughput sequencing for 3C) data which includes functions for:

1. Loading single-cell HiC datasets
2. Screening valid single-cell data.
3. Smoothing the contact maps with linear convolution, random walk or network enhancing
4. Calculating pairwise similarity using measures include InnerProduct, HiCRep and Selfish
5. Calculating embeddings for single cell HiC datasets efficiently with MDS, t-SNE and PHATE
6. Clustering the cells using scHiCluster, k-means and spectral clustering.
7. Visualizing embedded cells via 2-D or 3-D scatter plot.



# 2024.04.26

- Hic_contact file format transform

~~~bash
java -Xmx48000m  -Djava.awt.headless=true -jar   ~/Desktop/Analysis_work/App/juicer-1.6/scripts/common/juicer_tools.jar pre --threads 16 GSM6081061_Cell_ID_08.hic.txt Cell_ID_08.hic hg38.genome
~~~





