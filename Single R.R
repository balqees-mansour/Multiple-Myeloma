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

monaco.ref <- celldex::MonacoImmuneData()

# What is the structure of mimd.se?
monaco.ref
View(monaco.ref)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
library("SingleR")

Squery.data <- readRDS("/home/balqees/Multiple-Myeloma/mm-project/myeloma_subset.rds")

#Save the handeled reference and use it 
SeuratREF <- saveRDS(seurat_ref, "/home/balqees/Multiple-Myeloma/Single R/SeuratREF.rds")
#Convert Seurat object to Single cell experiment
Melnoma.sce <- as.SingleCellExperiment(DietSeurat(Squery.data))
Ref.sce <- as.SingleCellExperiment(DietSeurat(seurat_ref))

test_assay <- GetAssayData(Squery.data)
ref_assay <- GetAssayData(Ref_local)


# Predict using the Main-grained labels
Melnoma.main <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.main)
#Melnoma <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = Ref.sce,labels = Ref.sce$Cell_label)
#Melnoma <- SingleR(test = test_assay ,ref = ref_assay ,labels = Ref_local$Cell_label)
install.packages("viridis")
library("viridis")
install.packages("pheatmap")
library("pheatmap")
plotScoreHeatmap(Melnoma.main,
                 max.labels = 21,
                 clusters = Squery.data$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)
# Predict using the Fine-grained labels
Melnoma.fine <- SingleR(test = Melnoma.sce ,assay.type.test = 1,ref = monaco.ref,labels = monaco.ref$label.fine)
plotScoreHeatmap(Melnoma.fine,
                 max.labels = 21,
                 clusters = Squery.data$seurat_clusters,
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

#Let’s add the annotations to the Seurat object metadata so we can use them:
#Melnoma@meta.data$monaco.main <- monaco.main$pruned.labels
#Melnoma@meta.data$monaco.fine <- monaco.fine$pruned.labels

# Combine the labels for improved accuracy
pred.comb <- combineCommonResults(
  list(
    "Broad" = Melnoma.main,
    "Fine" = Melnoma.fine
  )
)

plotScoreHeatmap(pred.comb,
                 max.labels = 14,
                 clusters = Squery.data$seurat_clusters,
                 order.by = "clusters",
                 show_colnames = FALSE)

table(is.na(pred.comb$pruned.labels))
# FALSE  TRUE 
# 29996     4  


# Add the predicted labels and compare to clusters
# You are looking for clusters which have unambiguous labels
Squery.data$predicted_id <- pred.comb$pruned.labels
table(Squery.data$predicted_id, Squery.data$seurat_clusters)
saveRDS(Squery.data, file = "/home/balqees/Multiple-Myeloma/Single R/SMelnoma.combined.rds")
DimPlot(Squery.data)
DimPlot(Squery.data, group.by = "predicted_id", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()

#DimPlot(raw.data, label.size = 5, label = TRUE) 
# Select final IDs
# new.ids <- c("CD8+ T", "CD4+ T", "B", "progintor", "CD4+ T", "Macrophages", 
#              "progintor","progintor", "progintor", "progintor", "progintor",
#              "Keratinocytes", "Fibroplasts", "NK", "B", "CD8+ T",
#              "Endothelial", "progintor", "Endothelial.", "B", "DC")
# new.ids
# names(new.ids) <- levels(Melnoma)
# names(new.ids)
# levels(Melnoma)
# Melnoma <- RenameIdents(Melnoma, new.ids)
# Melnoma$cell_type <- Idents(Melnoma)
# #Visualize
# DimPlot(Melnoma)
# DimPlot(Melnoma, label.size = 5, label = TRUE) 
# #Compare Seurat clusters versus single R final annotation
# DimPlot(Melnoma, group.by = "seurat_clusters") + DimPlot(Melnoma)
# saveRDS(Melnoma, file = "D:/CELLCHAT2/Dr Abdelrahman TASK/Editting/Melnoma.single.r.rds")
# Melnoma <- Melnoma.single.r
#Convert from Single cell expriment format to seurat format
# library(Seurat)
# library("scater")
# library("loomR")
# library(SeuratDisk)
# library(SeuratData)
# library(patchwork)
# library(dplyr)
# library(Signac)
# Melnoma.seurat <- as.Seurat(Melnoma, assay = "NULL", counts =  "counts ", data = "logcounts")