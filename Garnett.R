# First install Bioconductor and Monocle
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install()
BiocManager::install(c("monocle"))

# Next install a few more dependencies
BiocManager::install(c('DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db'))
library(monocle)

install.packages("devtools")
devtools::install_github("cole-trapnell-lab/garnett", force = TRUE)

library(garnett)

# devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")

#To check garnett version 
packageDescription("garnett")$Version

library(Seurat)
library(monocle3)
install.packages(monocle3)

# First install Bioconductor and Monocle 3
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")

BiocManager::install()

# Next install a few more dependencies
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

#read query
Querydata <- readRDS("D:/Multiple Myeloma project/ANNOTATION STEP BY TREEARCHES TOOL/garnet/final_Downsample_RCPA.rds")
#Integrated assay has zero couunts so i switch to RNA assay
DefaultAssay(object = Querydata) <- "RNA"

# Create cds object 
qdata <- GetAssayData(object = Querydata, slot = "counts")
celldata <- as.data.frame(Querydata@meta.data)
genedata <- as.data.frame(x = row.names(Querydata), row.names = row.names(Querydata))
colnames(genedata) <- "gene_short_name"

cds <- new_cell_data_set(qdata, cell_metadata = celldata, gene_metadata = genedata)

# cds1 <- new_cell_data_set(as(qdata, "dgCMatrix"),
#                           cell_metadata = celldata,
#                           gene_metadata = genedata)

# generate size factors for normalization later - Monocole 2
#MM_cds <- estimateSizeFactors(cds)

classifier <- readRDS("D:/Multiple Myeloma project/ANNOTATION STEP BY TREEARCHES TOOL/garnet/ceWhole_20191017.RDS")

library(org.Hs.eg.db)
#classify your cells 
MM_cds <- classify_cells(cds, classifier,
                           db = org.Hs.eg.db,
                           cluster_extend = FALSE,
                           cds_gene_id_type = "SYMBOL")
### Error: is(object = cds, class2 = "CellDataSet") is not TRUE
library(scran)
cds1 <- as.CellDataSet(cds)
class(cds)

update.packages()
#cds <- updateCDS(cds) 


