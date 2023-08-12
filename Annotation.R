install.packages("Seurat")
install.packages("SeuratData")
install.packages("Matrix")
install.packages("remotes")
install.packages("hdf5r")
install.packages("reticulate")
install.packages("anndata")
remotes::install_github("mojaveazure/seurat-disk")
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(reticulate)
library(rhdf5)
library("anndata")

# read the query data with 30000  cells and 2000 genes
#raw.data <- read_h5ad("/home/balqees/Downloads/downsample_2000_raw.h5ad")
#raw.data<- CreateSeuratObject(counts = t(as.matrix(raw.data$X)), meta.data = raw.data$obs)

raw.data <- readRDS("/home/balqees/Multiple-Myeloma/mm-project/myeloma_subset.rds")
raw.data@assays$RNA@counts
#======================= preprocess the query data ==================================

raw.data@meta.data$nCount_RNA[2000:3000] # data needs preprocessing and filtering 

which(raw.data@meta.data$nCount_RNA == 0)
dim(raw.data)

# QC and selecting cells for further analysis
raw.data[["percent.mt"]] <- PercentageFeatureSet(raw.data, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(raw.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#ata_que!y <- subset(raw.data, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5)

# Normalizing the query data 
data_query <- NormalizeData(raw.data, normalization.method = "LogNormalize", scale.factor = 10000)


data_query <- FindVariableFeatures(data_query , selection.method = "vst", nfeatures = 2000)
 
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_query), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_query)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
# scaling query data 
data_query <- ScaleData(data_query)

# perform dimentionality reduction 
data_query@meta.data
data_query <- RunPCA(data_query, features = VariableFeatures(object = data_query))
DimPlot(data_query, reduction = "pca")

data_query <- RunUMAP(data_query, features = VariableFeatures(object = data_query))
DimPlot(data_query, reduction = "umap")


#Transfer labels
features = intersect(rownames(seurat_ref@assays$RNA@counts),rownames(data_query@assays$RNA@counts))

length(features)
query.anchors = FindTransferAnchors(reference = seurat_ref,
                                    query = data_query, 
                                    dims = 1:20,
                                    features = features,
                                    reference.reduction = 'pca')

predictions <- TransferData(anchorset = query.anchors,
                            k.weight = 10 ,
                            refdata = seurat_ref@meta.data$Cell_label, 
                            dims = 1:20)
data_query = AddMetaData(data_query, metadata = predictions)
data_query@meta.data$prediction.score.CD8..central.memory.T.cells

ref-labels = DimPlot(object = seurat_ref,reduction = "umap", group.by = "Cell_label", label = FALSE, label.size = 1.5, repel = TRUE) 
ref-labels

query_labels= DimPlot(object = data_query, reduction = "pca", group.by = "predicted.id", label = FALSE, label.size =1.5,   repel = TRUE) 
query_labels

FeaturePlot(data_query, features = c("Early erythroid progenitor", "Pro-B cells"),  reduction = "umap", cols = c("lightgrey", "darkred"), ncol = 3) & theme(plot.title = element_text(size = 10))
raw.data@meta.data$predicted.id

seurat_ref$prediction.match <- data_query@meta.data$predicted.id == seurat_ref$Cell_label
table(seurat_ref$prediction.match)

raw.data@meta.data$predicted.id

#seurat_ref <- RunUMAP(seurat_ref, dims = 1:30, reduction = "pca", return.model = TRUE)
data_query <- MapQuery(anchorset = query.anchors, reference = seurat_ref, query = data_query,
                     refdata = list(celltype = "Cell_label"), reference.reduction = "pca", reduction.model = "umap")



p1 <- DimPlot(seurat_ref, reduction = "umap", group.by = "Cell_label", label = TRUE, label.size =3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(data_query, reduction = "umap", group.by = "predicted.id", label = TRUE,
              label.size =3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

#merge reference and query
seurat_ref$id <- 'reference'
data_query$id <- 'query'
refquery <- merge(seurat_ref, data_query)
refquery[["pca"]] <- merge(seurat_ref[["pca"]], data_query[["pca"]])
refquery <- RunUMAP(refquery, reduction = 'pca', dims = 1:50)
DimPlot(refquery, group.by = 'id', shuffle = TRUE)









