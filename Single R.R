library("Seurat")
library("RcppAnnoy")
library("sceasy")
library(reticulate)
# SINGLE R ----------------------------------------------------------------
#Original Guide: https://bioconductor.org/packages/release/bioc/vignettes/SingleR/inst/doc/SingleR.html
# i use Melnoma_clusters resulted from Seurat to use single R 
Melnoma = Melnoma_clusters
# Predict using the MonacoImmuneData from celldex
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celldex")
library("celldex")
library(devtools)




monaco.ref <- celldex::MonacoImmuneData()

# What is the structure of mimd.se?
monaco.ref
View(monaco.ref)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
library("SingleR")

Querydata <- readRDS("D:/Multiple Myeloma project/ANNOTATION STEP BY TREEARCHES TOOL/garnet/final_Downsample_RCPA.rds")

#Save the handeled reference and use it 
Ref <- readRDS("D:/Multiple Myeloma project/ANNOTATION STEP BY TREEARCHES TOOL/local.rds")
# reference modification gene names
#install.packages("Matrix")
library(Matrix)
symbols <- Ref@assays[["RNA"]]@meta.features[["feature_name"]]
matrix <- Matrix(Ref@assays$RNA@counts, sparse = TRUE)
library(Seurat)
rownames(matrix) <- symbols
seurat_ref <- CreateSeuratObject(counts =matrix, meta.data = Ref@meta.data)

SeuratREF <- saveRDS(seurat_ref, "D:/Multiple Myeloma project/ANNOTATION STEP BY TREEARCHES TOOL/SeuratREF.rds")
#Convert Seurat object to Single cell experiment
Melnoma.sce <- as.SingleCellExperiment(DietSeurat(Querydata))
Ref.sce <- as.SingleCellExperiment(DietSeurat(seurat_ref))

test_assay <- GetAssayData(Querydata)
ref_assay <- GetAssayData(seurat_ref)


# Predict using the Main-grained labels
Melnoma.main <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
Myloma <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = Ref.sce,labels = Ref.sce$Cell_label)
#Melnoma <- SingleR(test = test_assay ,ref = ref_assay ,labels = Ref_local$Cell_label)
install.packages("viridis")
library("viridis")
install.packages("pheatmap")
library("pheatmap")

plotScoreHeatmap(Myloma,
                 max.labels = 21,
                 clusters = Querydata$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
table(Myloma$pruned.labels)

 plotScoreHeatmap(Melnoma.main,
                max.labels = 21,
                clusters = Querydata$seurat_clusters,
                order.by = "clusters",
                show_colnames = FALSE)
 

# Predict using the Fine-grained labels
Myloma.fine <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
plotScoreHeatmap(Myloma.fine,
                 max.labels = 21,
                 clusters = Querydata$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
# hpca.main <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.main)
# hpca.fine <- SingleR(test = sce,assay.type.test = 1,ref = hpca.ref,labels = hpca.ref$label.fine)
# dice.main <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.main)
# dice.fine <- SingleR(test = sce,assay.type.test = 1,ref = dice.ref,labels = dice.ref$label.fine)
saveRDS(Melnoma.main, file = "/home/balqees/Multiple-Myeloma/Single R/Melnoma.main.rds")

table(Melnoma.main$pruned.labels)
# B cells       Basophils    CD4+ T cells    CD8+ T cells 
# 12410              67            2158            3736 
# Dendritic cells       Monocytes     Neutrophils        NK cells 
# 1961            3479             107            3120 
# Progenitors         T cells 
# 964            1993
saveRDS(Melnoma.fine, file = "/home/balqees/Multiple-Myeloma/Single R/Melnoma.fine.rds")
table(is.na(Melnoma.fine$pruned.labels))
# FALSE  TRUE 
# 29990    10
table(is.na(Melnoma.main$pruned.labels))
# FALSE  TRUE 
# 29995     5 

