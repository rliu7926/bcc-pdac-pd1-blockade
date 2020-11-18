preprocessing <- function(dataset, normalize = TRUE, scale = TRUE, pca = TRUE, heatmap = TRUE) {
  
  ## PARAMETERS
  ## dataset = Seurat object
  ## normalize = NormalizeData T/F
  ## scale = ScaleData T/F
  ## pca = RunPCA T/F
  ## heatmap = DoHeatmap T/F
  
  dataset <- dataset
  
  if (isTRUE(normalize)) {
    dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  }
  
  if (isTRUE(scale)) {
    dataset <- FindVariableFeatures(dataset, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(dataset)
    dataset <- ScaleData(dataset, features = all.genes)
  }
  
  if (isTRUE(pca)) {
    dataset <- RunPCA(dataset, features = VariableFeatures(object = dataset))
  }
  
  if (isTRUE(heatmap)) {
    DimHeatmap(dataset, dims = 1:24, cells = 500, balanced = TRUE)
  }
  
  return(dataset)
}



determine_PCs <- function(dataset) {
  ElbowPlot(dataset, ndims=50)
  DimHeatmap(dataset, dims=1:24, cells = 500, balanced=TRUE)
  DimHeatmap(dataset, dims=25:48, cells = 500, balanced=TRUE)
}



clustering <- function(dataset, PCs=20, resolution=0.3) {
  ## Find neighbors
  dataset <- FindNeighbors(dataset, dims = 1:PCs)
  
  ## Do clustering based on resolution 
  dataset <- FindClusters(dataset, resolution = resolution, algorithm = 3)
  dataset.markers <- FindAllMarkers(dataset, test.use = 'wilcox', min.pct = 0.25, logfc.threshold = 0.50, max.cells.per.ident = 500)
  top5 <- dataset.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
  DoHeatmap(dataset, features = top5$gene, size = 3)
  
  ## Get list of genes
  dataset.top25 <- dataset.markers %>% group_by(cluster) %>% top_n(n=25, wt=avg_logFC)
  vec.dataset.top25 <- dplyr::pull(dataset.top25, 7)
  View(vec.dataset.top25)
  
  ## Create UMAP
  dataset <- RunUMAP(dataset, dims = 1:20)
}



analyze <- function(dataset, normalize = TRUE, scale = TRUE, pca = TRUE, heatmap = FALSE, PCs=20, resolution=0.3) {
  
  preprocessing(dataset=dataset, normalize = normalize, scale = scale, pca = pca, heatmap = heatmap)
  clustering(dataset=dataset, PCs=PCs, resolution=resolution)
  
}