################################################################################
#                                                                              #
#                     Seurat - Guided Clustering Tutorial                      #
#                                                                              #
################################################################################


### Install Package ###

# install.packages("Seurat")

library(dplyr)
library(Seurat)
library(patchwork)

# Load PBMC dataset
pbmc.data <- Read10X(data.dir = "C:\\Users\\Ryan\\Documents\\filtered_gene_bc_matrices\\hg19\\")

# Initialize Seurat object with raw (non-normalized) data
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc


### Standard pre-processing workflow ###
## QC = Quality Control 

# [[ operator adds columns to object metadata -- used to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature-Feature relationships 
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


### Normalizing Data ###

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# pbmc <- NormalizeData(pbmc) can produce the same results


### Identification of Highly Variable Features (Feature Selection) ###

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# Plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


### Scaling Data ###

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


### Perform Linear Dimensional Reduction ###
## PCA = Principal Components Analysis (Dimensionality Reduction)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and Visualize PCA results

print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


### Determine Dimensionality of Dataset ###

pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)


### Cluster the Cells ###

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)


### Run Non-Linear Dimensional Reduction (UMAP/tSNE) ###

# Install UMAP
reticulate::py_install(packages = "umap-learn")

pbmc <- RunUMAP(pbmc, dims = 1:10)

# Individual Clusters
DimPlot(pbmc, reduction = "umap")
saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")


### Finding Differentially Expressed Features (Cluster Biomarkers) ###

# Find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# Find markers for ever cluster compared to all remainign cells, report only posiive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# ROC test (classification power)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Visualizing marker expression
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
FeaturePlot(pbmc, features = c("MS4A1","GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# Expression heatmap for cells and features
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()


### Assigning Cell Type Identity to Clusters ###

new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
savdRDS(pbmc, file = "../output/pbmc3k_final.rds")






