
# CRISPR notebook


记录`CRISPR`相关的资料，包括概念、实验策略、分析算法等

参考资料：

+ [awesome-CRISPR](https://github.com/davidliwei/awesome-CRISPR)
+ [chronos](https://github.com/broadinstitute/chronos)
+ [High-content CRISPR screening](https://www.nature.com/articles/s43586-021-00093-4)


## 文件夹说明

一些自己归纳、分析整理的文件初稿


### assay01

基于已发表文章和数据，进行分析结果重现和研究思路学习

+ GSE126486

	CRISPR在合成致死领域的应用，利用CRISPR技术发现合成致死基因： [Pre-clinical activity of combined LSD1 and mTORC1 inhibition in MLL-translocated acute myeloid leukaemia](https://www.nature.com/articles/s41375-019-0659-6)

+ GSE182185

	`drugZ`算法分析CRIPR数据，用于drug resistance研究
	
	[CRISPR screens uncover protective effect of PSTK as a regulator of chemotherapy-induced ferroptosis in hepatocellular carcinoma](https://molecular-cancer.biomedcentral.com/articles/10.1186/s12943-021-01466-9)
	
	** drugz **： https://github.com/hart-lab/drugz
	[Identifying chemogenetic interactions from CRISPR screens with drugZ](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-019-0665-3)
	
	
+ 。。。。待续



### assay02

比较 `MAGECK` 和 `Chronos` 在 gene fitness 评价上的一致性和相关性

+ DLD1 细胞系

+ HT29 细胞系

+ MCF7 细胞系

+ MDMAB231 细胞系

+ U2OS 细胞系

### chronos

chronos 分析脚本： [chronos](https://github.com/broadinstitute/chronos)

具体分析实例，可以参考 `assay02/` 或 `chronos.ipynb`


### chronos_vs_mageck

chronos_vs_mageck 比较分析结果和报告


### CRISPR_papers

CRISPR 相关已发表文章


### essential_noessential

CRISPR 分析中的阳性对照基因


