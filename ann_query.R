library(Matrix)

data.ref <- readRDS("/home/balqees/Multiple-Myeloma/mm-project/local.rds")
colnames(data.ref)
length(unique(colnames(data.ref)))
colnames(data.ref@assays$RNA@counts)
colnames(data.ref@meta.data)
rownames(data.ref)


symbols <- data.ref@assays[["RNA"]]@meta.features[["feature_name"]]
matrix <- Matrix(data.ref@assays$RNA@counts, sparse = TRUE)
rownames(matrix) <- symbols
seurat_ref <- CreateSeuratObject(counts =matrix, meta.data = data.ref@meta.data)
rownames(seurat_ref)
data.ref@assays[["RNA"]]@meta.features
colnames(seurat_ref)
# determine the ident by cell type 
Idents(seurat_ref) <- "Cell_label"

# How many cells are in each cluster
table(Idents(seurat_ref))

# What proportion of cells are in each cluster?
prop.table(table(Idents(seurat_ref)))

# What are the cell names of all NK cells?
WhichCells(seurat_ref, idents = "Small pre-B cell")


#============================== preprocessing ==================================

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ref[["percent.mt"]] <- PercentageFeatureSet(seurat_ref, pattern = "^MT-")

# Visualize QC metrics as a violin plot
#VlnPlot(seurat_ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

#plot1 <- FeatureScatter(seurat_ref, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(seurat_ref, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 
plot2
# subset 
#seurat_ref <- subset(seurat_ref, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

seurat_ref <- NormalizeData(seurat_ref)
# FindVariableFeatures
seurat_ref <- FindVariableFeatures(seurat_ref , selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_ref), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_ref)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# Scaling the reference 
seurat_ref <- ScaleData(seurat_ref)
seurat_ref <- RunPCA(seurat_ref, features = VariableFeatures(object = seurat_ref))
# Examine and visualize PCA results a few different ways
print(seurat_ref[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(seurat_ref, reduction = "pca")
DimHeatmap(seurat_ref, dims = 1, cells = 500, balanced = TRUE)
# find clusters 
seurat_ref <- FindNeighbors(seurat_ref, dims = 1:10)
seurat_ref <- FindClusters(seurat_ref, resolution = 0.5)
seurat_ref <- RunUMAP(seurat_ref, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_ref, reduction = "umap")

DimPlot(object = seurat_ref,reduction = "umap", group.by = "Cell_label", label = TRUE, label.size = 3, repel = TRUE) + NoLegend()
DimPlot(object = seurat_ref,reduction = "umap", group.by = "Cell_label", label = FALSE, label.size = 3, repel = TRUE) 



