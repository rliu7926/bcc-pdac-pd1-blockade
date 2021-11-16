################################################################################
#                                                                              #
#       Defining Epidermal Basal Cell States During Skin Homeostasis and       #
#              Wound Healing Using Single-Cell Transcriptomics                 #
#                                                                              #
################################################################################

library(dplyr)
library(Seurat)
library(patchwork)

GSE_data <- Read10X(data.dir = "C:\\Users\\Ryan\\Downloads\\GSE142471_RAW\\Un-Wounded_1")
unwounded_1 <- CreateSeuratObject(counts = GSE_data, project = "GSE142471", min.cells = 3, min.features = 200)
unwounded_1

unwounded_1[c("Krt14", "Mki67", "Krt1", "Krt17", "Cd34", "Postn", "Lhx2", "Lgr5"), 1:30]

unwounded_1[["percent.mt"]] <- PercentageFeatureSet(GSE_data, pattern = "^MT-")

VlnPlot(unwounded_1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot <- FeatureScatter(unwounded_1, feature1 = "nCount_RNA", "nFeature_RNA")
plot

unwounded_1 <- NormalizeData(unwounded_1, normalization.method = "LogNormalize", scale.factor = 10000)

unwounded_1 <- FindVariableFeatures(unwounded_1, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(unwounded_1), 10)

plot1 <- VariableFeaturePlot(unwounded_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(unwounded_1)
unwounded_1 <- ScaleData(unwounded_1, features = all.genes)

unwounded_1 <- RunPCA(unwounded_1, features = VariableFeatures(object = unwounded_1))
print(unwounded_1[["pca"]], dims = 1:10, nfeatures = 10)

VizDimLoadings(unwounded_1, dims = 1:2, reduction = "pca")
DimPlot(unwounded_1, reduction = "pca")
DimHeatmap(unwounded_1, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(unwounded_1, dims = 1:10, cells = 500, balanced = TRUE)

unwounded_1 <- JackStraw(unwounded_1, num.replicate = 100)
unwounded_1 <- ScoreJackStraw(unwounded_1, dims = 1:20)
JackStrawPlot(unwounded_1, dims = 1:20)
ElbowPlot(unwounded_1)

unwounded_1 <- FindNeighbors(unwounded_1, dims = 1:10)
unwounded_1 <- FindClusters(unwounded_1, resolution = 0.6)
head(Idents(unwounded_1), 5)

unwounded_1 <- RunUMAP(pbmc, dims = 1:10)
DimPlot(unwounded_1, reduction= "umap")

unwounded_1 <- RunUMAP(unwounded_1, dims = 1:10)
DimPlot(unwounded_1, reduction = "umap")

cluster1.markers <- FindMarkers(unwounded_1, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

cluster5.markers <- FindMarkers(unwounded_1, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n=5)

unwounded_1.markers <- FindAllMarkers(unwounded_1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
unwounded_1.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

unwounded_1.markers <- FindMarkers(unwounded_1, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(unwounded_1, features = c("Atf3"))
VlnPlot(unwounded_1, features = c("Atf3"), slot = "counts", log = TRUE)
FeaturePlot(unwounded_1, features = c())

top10 <- unwounded_1.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(unwounded_1, features = top10$gene) + NoLegned()

