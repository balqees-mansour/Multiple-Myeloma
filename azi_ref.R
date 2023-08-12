install.packages("Annoy")
library(Annoy)
library(Azimuth)
library(hdf5r)
library(SeuratDisk)
devtools::install_version("RcppAnnoy", version = "0.0.16")
# generating Azimuth reference from seurat object 
# https://github.com/satijalab/azimuth/wiki/Generating-an-Azimuth-Reference
seurat_ref@meta.data$Cell_label
seurat_ref@meta.data$Cluster_ID
seurat_ref@meta.data$cell_type

az_ref <- seurat_ref
#Step 1: Perform SCTransform normalization
az_ref <- SCTransform(az_ref)
# Store the normalized data in the SCT assay
az_ref[["SCT"]] <- az_ref
ref <- AzimuthReference(
  object = az_ref,
  refUMAP = "umap",
  refDR = "pca",
  metadata = c("Cell_label", "Cluster_ID"),
  dims = 1:50,
  k.param = 31
)
saveRDS(ref, file = "ref_Z.Rds")

#Saving the annoy index
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = "idx.annoy")

#==========================================================================================
# The RunAzimuth function can take a Seurat object as input
DimPlot(ref, group.by ="Cell_label", raster=FALSE,
        label = TRUE, repel = TRUE) + NoLegend() 
ref@meta.data$Cell_label

mm.query <- RunAzimuth(raw.data, 
                       reference = "/home/balqees/Multiple-Myeloma/azimuth_ref")

View(mm.query@meta.data)
mm.query@meta.data$predicted.Cell_label
mm.query@meta.data$predicted_id
DimPlot(mm.query, group.by = "predicted.Cell_label", 
        label = TRUE, label.size = 3, repel = TRUE ) + NoLegend()

# Set Idents
Idents(mm.query)

DimPlot(mm.query)

Idents(mm.query) <- "predicted.Cell_label"

DimPlot(mm.query)
DimPlot(ref_z)

# Azimuth normalizes data before mapping, but does not return the results
# normalize the data here before visualization.
mm.query <- NormalizeData(mm.query)

FeaturePlot(mm.query, features = c ("NARF", "SDC1,", "PTP4A3", "COL1A2"))

# Save RDS
saveRDS(mm.query, "../Desktop/Video_Tutorials/Data/lungquery.rds")

ref_z <- readRDS("/home/balqees/Multiple-Myeloma/ref.Rds")



