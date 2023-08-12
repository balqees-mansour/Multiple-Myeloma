# packages installation
update.packages(oldPkgs = c("withr", "rlang"))
install.packages("SeuratObject")
remotes::install_github("hhoeflin/hdf5r")
remotes::install_github("mojaveazure/seurat-disk")
remotes::install_github('satijalab/azimuth', ref = 'master')
install.packages(c( "‘hdf5r’, ‘SeuratDisk’"))
library(Azimuth)
library(hdf5r)
library(SeuratDisk)

##==================================================
library(Seurat)
library(Azimuth)
library(tidyverse)

# #plot the reference data
rownames(seurat_ref)

view(seurat_ref@meta.data)

DimPlot(seurat_ref, group.by ="Cell_label", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend()

DimPlot(lungref, group.by ="ann_finest_level", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend()

raw.data #query data 

view(data_query@meta.data)

# The RunAzimuth function can take a Seurat object as input
mm.query <- RunAzimuth(raw.data, 
                       reference = "/home/balqees/Multiple-Myeloma/reference")

view (mm.query@meta.data)

DimPlot(mm.query, group.by = "predicted.celltype.l1", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()

DimPlot(mm.query, group.by = "predicted.celltype.l2", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()
mm.query@meta.data$predicted.celltype.l1
mm.query@meta.data$predicted.celltype.l2
# Set Idents
Idents(mm.query)

DimPlot(mm.query)

Idents(mm.query) <- "predicted.celltype.l2"

DimPlot(mm.query)
DimPlot(ref_z)

# Azimuth normalizes data before mapping, but does not return the results
# normalize the data here before visualization.
mm.query <- NormalizeData(mm.query)

FeaturePlot(mm.query, features = c("SDC1", "SLAMF7", "PTP4A3", "XBP1"))

# Save RDS
saveRDS(mm.query, "../Desktop/Video_Tutorials/Data/lungquery.rds")

ref_z <- readRDS("/home/balqees/Multiple-Myeloma/ref.Rds")






