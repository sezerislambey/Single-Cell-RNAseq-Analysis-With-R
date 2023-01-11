# working with the Smart-Seq2 dataset

print ("hello workshop")

library(Seurat)
library(tidyverse)

# set your working directory
setwd("/home/klaus/Documents/Bioinformatics/workshop/qcb_single_cell")

# load the dataset
mydata <- read.table("GSE102130_K27Mproject.RSEM.vh20170621.txt", header = T, row.names = 1)

# look at the first 6 rows
head(mydata)

dim(mydata)

# create the Seurat object
glioma <- CreateSeuratObject(counts = mydata, min.cells = 300, min.features = 200)

# the "counts" slot
GetAssayData(object = glioma, slot = "counts")

# check the class of the object
class(glioma)

# Violin plots - based on nGenes and nUMIS
VlnPlot(object = glioma, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# make the GenePlot - number of genes vs number of UMIS I
FeatureScatter(object = glioma, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# check the dimension of the data
dim(glioma@assays$RNA)

# a QC step - filter cells having less than 500 and more than 8000 genes
glioma <- subset(x = glioma, subset = nFeature_RNA > 500 & nFeature_RNA < 8000)

# check again the dimension of the data
dim(glioma@assays$RNA)

# scale the data
all.genes <- rownames(x = glioma)
glioma <- ScaleData(object = glioma, features = all.genes)

# glioma@assays$RNA@counts
# glioma@assays$RNA@data

# do some checks
expr1 <- GetAssay(glioma, slot="counts")[,"MGH66. P08.BØ5"]
# first, let's see the raw expression value of the gene A2M of the cell with name MGH66.P08.B05
expr1["A2M",] # should be equal to 18.78
log1p(expr1["A2M" ,] * 10000 / sum(expr1)) 

glioma <- NormalizeData(object = glioma, normalization.method = "LogNormalize", scale.factor = 10000)

# Now, the normalized data is stored in the slot counts
# check what is the normalized expression of the same gene in the same cell
expr2 <- GetAssay(glioma, slot="counts")[, "MGH66.P08.B05"]
expr2["A2M",] # should be equal to 0.1721429
# that is how it was computed!
# should give the same result

# glioma[, "MGH66.P08.B05"]

# finding variable genes
# Expmean function: Calculate mean of logged values in non-log space (return answer in log-space)
# LogVMR function: Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
glioma <- FindVariableFeatures(object = glioma, selection.method = "vst", nfeatures = 2000)
top10 <- head(x = VariableFeatures(object = glioma), 10)

# VariableFeatures(glioma)

plot1 <- VariableFeaturePlot(object = glioma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
LabelPoints(plot = plot1, points = top10, repel = TRUE) 



# Dimensionality reduction step by using the Principal Component Analysis
glioma <- RunPCA(object = glioma, features = VariableFeatures(object = glioma))

print(x = glioma[["pca"]], dims = 1:15, nfeatures = 20)

# The results of dimensionality reduction step are stored in the slot dr
glioma@reductions$pca

# glioma@reductions$pca@cell.embeddings
# glioma@reductions$pca@cell.embeddings[, c(1, 2)]
# glioma@reductions$pca@feature.loadings

# this is how to check the "embeddings" of the cells into the PCA 2D space (ba)
# the new coordinates)
head(glioma@reductions$pca@cell.embeddings)

# this is how to check the standard deviation of the PCs
# basically, it is the standard deviation of the first "new" coordinate across the cell
glioma@reductions$pca@stdev # check the standard deviation
plot(glioma@reductions$pca@stdev)

# plot(glioma@reductions$pca@feature.loadings)

# vizualise the results of the PCA step
VizDimLoadings(object = glioma, dims = 1:2, reduction = "pca")

# Draws a heatmap focusing on a principal component.
# Both cells and genes are sorted by theirprincipal component scores.
# Allows for nice visualization of sources of heterogeneity in the dataset
# do.balanced = TRUE -> Return an equal number of cells with both + and - PC scores
DimHeatmap(object = glioma, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(object = glioma, dims = 1:6, cells = 500, balanced = TRUE)

# check the top cells by the first PC
TopCells(glioma)
TopCells(glioma, dim = 2)

# Coordinate of a particular cell (the "top" cell)
Embeddings(glioma)["BCH836.P01.D12",]

# indeed, it's first PC coordinate is the largest by the absolute value
max(abs(Embeddings(glioma)[, "PC_2"]))
# top PC genes
TopFeatures(glioma)
TopFeatures(glioma, dim = 2)
TopFeatures(glioma, 50)

# PCA loadings of a gene, i.e. how much the gene contributes to each component
Loadings (glioma)["OLR1",]

# JackStraw procedure - find the significant PCs
#glioma <- JackStraw(object = glioma, num.replicate = 100)
#glioma <- ScoreJackStraw(object = glioma, dims = 1:20)
#JackStrawPlot(object = glioma, dims = 1:20)


ElbowPlot(glioma)
ElbowPlot(glioma, ndim=50)

# Find the clusters - "cell types"
glioma <- FindNeighbors(object = glioma, dims = 1:8)
glioma <- FindClusters(object = glioma, resolution = 0.01)
head(x = Idents(object = glioma), 10)
Idents(glioma)
table(Idents(object = glioma)) 

# a helper function to get the cell of a specific identity
WhichCells(object = glioma, ident = 1) 


# or, select the cells which have at least 7600 genes
WhichCells(object = glioma, subset.name = c("nGene"), accept.low = 7600)

# Find average expression in each of the cluster
AverageExpression(object = glioma) 


# build the cluster tree
glioma <- BuildClusterTree(glioma)
# plot cluster tree
PlotClusterTree(glioma)


# how many cells of each cell type do we have?
table(Idents(glioma))

# run the tSNE visualization
glioma <- RunTSNE(object = glioma,  dims.use = 1:8,  do.fast = TRUE, seed.use = 26)

# and see the SNE plot I
TSNEPlot(object = glioma) 
DimPlot(glioma, reduction = "pca") 
# DimPlot(glioma, reduction = "pca", dims = c(3, 5))

  
# Differential expression analysis
# Find the marker genes for each cluster vs all the others
markers <- FindAllMarkers(object = glioma,  only.pos = TRUE,  min.pct = 0.25,  thresh.use = 0.25)
  
head(markers)
dim(markers)
table(markers$cluster)
markers[markers$cluster == 1,]

# Top markers for each cluster
# topGenes <- markers %>% group_by(cluster) %>% top_n(1, wt=avg_log2FC) %>% select(gene)
VlnPlot(object = glioma, features = c("A2M", "KCNQ2"))
VlnPlot(object = glioma, features = c("PTPRZ1"))
VlnPlot(object = glioma, features = c("SLAIN1"))
rownames(glioma)
VariableFeatures(glioma)
rownames(markers)
setdiff(VariableFeatures(glioma), rownames(markers))
VlnPlot(object = glioma, features = c("SKA3"))

# Colors single cells on a dimensional reduction plot
# according to a 'feature' (i.e. gene expression, PC scores, number of genes detected, etc.)
FeaturePlot(glioma,  c("BCAN"),  cols = c("yellow", "red", "blue"),  ncol = 3)
FeaturePlot(glioma,  c("SKA3"),  cols = c("yellow", "red", "blue"),  ncol = 3)
FeaturePlot(glioma,  c("PTPRZ1"),  cols = c("yellow", "red", "blue"),  ncol = 3)


# Read10X for 10X data set
mydata <- Read10X("filtered_gene_bc_matrices/GRCh38/")
glioma <- CreateSeuratObject(counts = mydata)
glioma <- subset(x = glioma, subset = nFeature_RNA > 500 & nFeature_RNA < 3000)
all.genes <- rownames(x = glioma)
glioma <- ScaleData(object = glioma, features = all.genes)
glioma <- NormalizeData(object = glioma, normalization.method = "LogNormalize", scale.factor = 10000)
glioma <- FindVariableFeatures(object = glioma, selection.method = "vst", nfeatures = 2000)
plot1 <- VariableFeaturePlot(object = glioma)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
LabelPoints(plot = plot1, points = top10, repel = TRUE)
glioma <- RunPCA(object = glioma, features = VariableFeatures(object = glioma))
print(x = glioma[["pca"]], dims = 1:15, nfeatures = 20)
VizDimLoadings(object = glioma, dims = 1:4, reduction = "pca")
glioma <- FindNeighbors(object = glioma, dims = 1:8)
glioma <- FindClusters(object = glioma, resolution = 0.2)
glioma <- RunTSNE(object = glioma,  dims.use = 1:8,  do.fast = TRUE, seed.use = 26)
TSNEPlot(object = glioma)
DimPlot(glioma, reduction = "pca") 
markers <- FindAllMarkers(object = glioma,  only.pos = TRUE,  min.pct = 0.25,  thresh.use = 0.25)
FeaturePlot(glioma,  c("IL7R"),  cols = c("yellow", "red", "blue"),  ncol = 3)
FeaturePlot(glioma,  c("MS4A1"),  cols = c("yellow", "red", "blue"),  ncol = 3)





  



