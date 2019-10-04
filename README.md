# RNA_SEQ_ANALYSIS
Tutorial for RNA_seq analysis
Step by step flow on processing of data

getwd()
library(dplyr)
library(Seurat)

##loading the data
pbmc.data <- Read10X(data.dir = 'C:/Users/Tony Isebe/Desktop/Seurat Tutorial_Single_Cell RNA Seq/filtered_gene_bc_matrices/hg19')
View(pbmc.data)
?write.table
write.table(pbmc.data, file = 'pbmc.txt')
getwd()

#starting the Seurat object with the raw data
pbmc <- CreateSeuratObject(counts = pbmc.data, project = 'pbmc3k', min.cells = 3, 
                           min.features = 200)
pbmc

## adding columns to object metadata and stash QC stats

pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern = '^MT-')

#Visualize QC metrics as a violin plot

VlnPlot(pbmc, features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)

# Feature scatter to visualize feature-feature relationship

plot1 <- FeatureScatter(pbmc, feature1 = 'nCount_RNA', feature2 = 'percent.mt')
plot2 <- FeatureScatter(pbmc, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA')
CombinePlots(plots = list(plot1,plot2))

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt <5)
pbmc

##Data normalization, global scaling approach adopted, LogNormalize normalizes feature expression measurements for each cell by the total expression, multiplies this by a scale facror (10,000 by default) and log-transforms the result. Normalized values are stored in pbmc [['RNA]]@DATA

pbmc <- NormalizeData(pbmc, normalization.method = 'LogNormalize', scale.factor = 10000)
 #alternatively

pbmc <- NormalizeData(pbmc)

##Identifying highly variable features

pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)


##Identify the 10 most highly variable genes

top10 <- head(VariableFeatures(pbmc), 10)
top10


##plot variable features with and without labels

plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

##scaling data prior to running a PCA, shifts expression of each gene so that the mean expression across cells is 0, ensures highly expressed genes do not dominate

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- ScaleData(pbmc)

##removing unwanted sources of variation 

pbmc <- ScaleData(pbmc, vars.to.regress = 'percent.mt')

##performing linear dimensional reduction, PCA

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

##examining and visualizing PCA results a few different ways

print(pbmc[['pca']], dims = 1:10, nfeatures = 10)


#visualizing using VizDimLoadings

VizDimLoadings(pbmc, dims = 1:5, reduction = 'pca')

DimPlot(pbmc, reduction = 'pca')

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)


DimHeatmap(pbmc, dims = 1:5, cells = 500, reduction = 'pca')

## clustering  cells based on their PC scores

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:20)

ElbowPlot(pbmc)
## cluster of cells

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

##cluster IDs of first five cells

head(Idents(pbmc),5)


pbmc <- RunUMAP(pbmc, dims = 1:10)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("umap-learn")


saveRDS(pbmc, file = 'C:/Users/Tony Isebe/Desktop/Seurat Tutorial_Single_Cell RNA Seq/filtered_gene_bc_matrices/hg19_tutorial')
View(saveRDS())
##find differentially expressed features

cluster1.markers <- FindAllMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n=5)

cluster2.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster2.markers, n=5)
saveRDS(pbmc, file = 'C:/Users/Tony Isebe/Desktop/Seurat Tutorial_Single_Cell RNA Seq/filtered_gene_bc_matrices/hg19_tutorial_more')
