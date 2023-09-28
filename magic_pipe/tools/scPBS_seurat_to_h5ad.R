suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

all_args <- commandArgs()
seurat_obj_path <- all_args[6]
h5ad_path <- all_args[7]

seurat_obj <- readRDS(seurat_obj_path)
# export top 5000 highly variable genes
seurat_obj %>% FindVariableFeatures(., selection.method = "vst", nfeatures = 5000) -> seurat_vst
# MUST run DietSeurat before WriteH5AD
seu <- DietSeurat(
  seurat_vst,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
  #features = rownames(seurat_vst), # export all genes, not just top highly variable genes
  features = seurat_vst@assays$RNA@var.features, # export highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap"),
  graphs = c("RNA_nn", "RNA_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

MuDataSeurat::WriteH5AD(seu, h5ad_path, assay="RNA")
# sceasy, be sure Python is in your PATH and the Python package anndata is installed. Better is to run this code in a conda env
# sceasy::convertFormat(seurat_obj, from="seurat", to="anndata", outFile=h5ad_path)
