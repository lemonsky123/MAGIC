suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

all_args <- commandArgs()
seurat_obj_path <- all_args[6]
cor_top500_path <- all_args[7]
pdf_path <- all_args[8]


pcw_new <- readRDS(seurat_obj_path)
cor_top500_file <- fread(cor_top500_path, header = T)
cor_top500 <- cor_top500_file[1:500, -c(2, 3)]
cor_top500 <- as.list(cor_top500)

pcw_new <- AddModuleScore(pcw_new, features = cor_top500 ,name = "rvTRS")
pcw_new$rvTRS1[which(pcw_new$rvTRS1<0)]<-0

pdf(pdf_path)
FeaturePlot(pcw_new, features = "rvTRS1", cols =c("lightgrey", "#DC143C")) + ggtitle("cor_top500")
dev.off()
# ggsave(file = pdf_path)
