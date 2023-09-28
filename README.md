# MAGIC

<div align=center>
<img src="https://github.com/sulab-wmu/MAGIC/blob/main/pic/MAGIC.jpg" width="350" height="140"/>
</div>

The **MAGIC** (*Myopia Associated Genetics and Intervention Consortium*) toolkits are a collection of Python scripts for rare variant-related analysis. The toolkits are composed of four subcommands, they are **vcf2mat** (convert the GATK vcf file into a variant by sample matrix), **RVAS** (fisher test, permutation test, SKAT and ACAT analysis), **logistic**, and **scPBS**, respectively. The **scPBS** (*single-cell Polygenic Burden Score*) is a method to evaluate polygenic burden enrichment of rare variants in individual cells of scRNA-seq data to consider heterogeneity within each cell type. <br />

![Image text](https://github.com/sulab-wmu/MAGIC/blob/main/pic/scPBS.jpg)

**scPBS** assesses whether a given cell has an excess of rare PTVs among EM cases in cell-specific highly expressed genes derived from scRNA-seq using an appropriately matched empirical null distribution to estimate well-calibrated P values. scPBS consists of the following three steps. 

**_(1)._** scPBS constructs a set of cell-specific highly expressed genes, the most specific genes in each individual cell from high cell-to-cell variable features in the dataset. <br />
**_(2)._** scPBS selects qualifying variants (QV) that pass all filters of a specific model and builds the geneset(cell)-by-individual indicator matrix used for the polygenic burden test. <br /> 
**_(3)._** scPBS quantifies the EM-related polygenic burden coefficient (odds ratio) of highly expressed genes in each cell to generate single-cell specific polygenic burden scores (scPBS). Following previous studies, we correlate scPBS with the expression level of each gene across cells and prioritize the trait-relevant genes by ranking the Pearson correlation coefficients (PCCs). <br />

Finally, a rare variant trait-relevant score (rvTRS) of each cell is computed by averaging the expression level of the top 10% trait-relevant genes based on ranked PCCs and subtracting the random control cell score via the cell-scoring method used in Seurat.

**[We also provide a step-by-step tutorial for you to walk through the scPBS pipeline](https://github.com/sulab-wmu/MAGIC/blob/main/tutorial.md)**

_If you have any questions using our toolkits, please reach out via <likai@ucas.ac.cn> and <yuanjian0415@gmail.com>._
