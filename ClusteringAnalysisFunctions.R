dist_table <- function(dataset, uniqueclusters1, uniqueclusters2, loc1, loc2) {
  
  ## Datset should be a dataframe 
  ## Uniqueclusters 1 and uniqueclustesr2 are matrices
  ## loc1 and loc2 describe position of cluster info to be compared in dataset
  
  distribution.table <- data.frame(matrix(0, nrow=nrow(uniqueclusters1), ncol=nrow(uniqueclusters2)))
  for (cluster1 in 1:nrow(uniqueclusters1)) {
    for (cluster2 in 1:nrow(uniqueclusters2)) {
      for (row in 1:nrow(dataset)) {
        if ((dataset[row,loc1] == uniqueclusters1[cluster1,]) & (dataset[row,loc2] == uniqueclusters2[cluster2,])) {
          distribution.table[cluster1, cluster2] <- distribution.table[cluster1, cluster2] + 1
        }
      }
    }
  }
  return(distribution.table)
}

prop_table <- function(distribution.table, uniqueclusters1, uniqueclusters2){
  distribution.table <- as.matrix(distribution.table)
  
  rownames(distribution.table) <- uniqueclusters1
  colnames(distribution.table) <- uniqueclusters2
  
  proportions.table <- prop.table(distribution.table, 1)
  proportions.table <- as.data.frame(proportions.table)
  proportions.table <- unlist(proportions.table)
}



ggplot_heatmap <- function(proportions.table, uniqueclusters1, uniqueclusters2, xlab, ylab) {
  x <- uniqueclusters1
  y <- uniqueclusters2
  data <- expand.grid(X=x, Y=y)
  data$Z <- proportions.table
  ggplot(data, aes(X, Y, fill= Z)) + geom_tile() + labs(x=xlab, y=ylab) + theme(axis.text.x = element_text(angle=45, hjust=1))
}


ggplot_barcharts <- function(dataset.batch, dataset.clusters, colors="Spectral") {
  
  dataframe <- data.frame("Batch" = dataset.batch, "Clusters" = dataset.clusters)
  tibble <- as_tibble(dataframe)
  
  summary <- tibble %>% group_by(Clusters, Batch) %>% count()
  
  normalized_dataset <- summary %>% group_by(Batch) %>% mutate(normalize = n/sum(n))

  colorCount = length(unique(dataset.batch))
  getPalette = colorRampPalette(brewer.pal(9, colors))
  
  ggplot(normalized_dataset, aes(x=Clusters, y=normalize, fill=Batch)) + 
      geom_col(position = 'fill') + 
      scale_fill_manual(values= getPalette(colorCount)) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) + 
      ylab("Cell Proportion Normalized by Batch")
}


bioc_markers <- function(dataset, celltype, group, ident, pCutoff=0.05, FCcutoff=0.5) {
  
  ## Find differentially expressed markers for Enhanced Volcano plot
  
  celltype.markers <- FindMarkers(dataset, subset.ident=celltype, group.by=group, ident.1=ident)
  celltype.markers_tibble <- as_tibble(celltype.markers, rownames='GeneNames')
  celltype.markers_tibble <- arrange(celltype.markers_tibble, avg_logFC)
  
  return(celltype.markers_tibble)
}


bioc_volcano <- function(dataset, celltype, group, ident, ident2, pCutoff=0.05, FCcutoff=0.5) {
  
  ## Create Enhanced Volcano plots from Bioconductor
  ## celltype: cluster (ex. macrophages)
  ## group: how to divide cluster (ex. responders vs. nonresponders)
  ## ident: subcluster to be compared to rest of group (ex. responders)
  ## ident2: rest of group (just used to generate title)
  
  celltype.markers <- FindMarkers(dataset, subset.ident=celltype, group.by=group, ident.1=ident)
  celltype.markers_tibble <- as_tibble(celltype.markers, rownames='GeneNames')
  celltype.markers_tibble <- arrange(celltype.markers_tibble, avg_logFC)
  
  celltype.volcano <- EnhancedVolcano(celltype.markers_tibble,
                                      lab = celltype.markers_tibble$GeneNames,
                                      x = 'avg_logFC', y = 'p_val_adj',
                                      pCutoff = pCutoff, FCcutoff = FCcutoff,
                                      title = paste(celltype, "DE:", ident, "vs.", ident2))
  
  return(celltype.volcano)
  
}
