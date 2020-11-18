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



ggplot_heatmap <- function(proportions.table, uniqueclusters1, uniqueclusters2, palette, xlab, ylab) {
  x <- uniqueclusters1
  y <- uniqueclusters2
  data <- expand.grid(X=x, Y=y)
  data$Z <- proportions.table
  ggplot(data, aes(X, Y, fill= Z)) + geom_tile() + scale_fill_distiller(palette=palette) + labs(x=xlab, y=ylab)
}


