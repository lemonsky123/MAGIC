# INSTALLATION
We recommend installing scPBS inside a conda environment.
```
conda create -n scPBS python=3.9
conda activate scPBS
git clone https://github.com/sulab-wmu/MAGIC.git
cd MAGIC
pip install . --ignore-installed
```
# USAGE
## scPBS
Before using scPBS, we suggest you walk through the [Seurat tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html) since scPBS takes the Seurat object as input. After clustering cells in Seurat, please save the Seurat object to a rds file. Then we need to convert the rsd file to h5ad format. Note that you must also provide the path to the Rscript executable with Seurat and dplyr installed.
```
magic rds2h5ad -r seurat_obj.rds -o seurat_obj.h5ad -rscript /usr/bin/Rscript
```
If you have some trouble using rds2h5ad, you could also convert by yourself. [This article](https://zqfang.github.io/2020-04-28-seurat2scanpy/) explains how to do this process in detail. 
```
library(Seurat)
library(dplyr)

# MuDataSeurat
seurat_obj <- readRDS("seurat_obj.rds")

seurat_obj %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 5000) -> seurat_obj
# MUST run DietSeurat before WriteH5AD
seu <- DietSeurat(
  seurat_obj,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  #features = rownames(seurat_obj), # export all genes, not just top highly variable genes
  features = seurat_obj@assays$RNA@var.features, # export highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

MuDataSeurat::WriteH5AD(seu, "seurat_obj.h5ad", assay="RNA")
```
After converting, run scPBS with the following commands.
```
magic scPBS -p rare_lof.txt -g sample_EM_class.txt -b seurat_obj.h5ad -o seurat_obj_corr_odds.txt
```
- -p rare_lof.txt  A snp x sample matrix, each row stands for an snp and each column is a sample, 0 is "0/0", 1 is "0/1" and 2 is "1/1". Other columns are "symbol", "ENSG", "ENST", "sample_stat", "variant_stat", "info". You can get this file by "magic vcf2mat".
- -g sample_EM_class.txt A two-column file, the first column is sample names and the second is sample types. The headers must be "sample" and "type". The "type" must be "case" or "control". Other words are not accepted.
- -b seurat_obj.h5ad The converted seurat object.
- -o Path to put the correlation results between odds ratio and gene expresssion.
After getting the correlation result. You could run the magic seurat_plot command to generate the Feature plot.
```
magic seurat_plot -r seurat_obj.rds -c seurat_obj_corr_odds.txt -o seurat_obj_feature.pdf -rscript /usr/bin/Rscript
```
