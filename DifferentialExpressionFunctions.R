bioc_markers <- function(dataset, celltype, group, ident, pCutoff=0.05, FCcutoff=0.5) {
  
  ## Find differentially expressed markers for Enhanced Volcano plot
  
  celltype.markers <- FindMarkers(dataset, subset.ident=celltype, group.by=group, ident.1=ident)
  celltype.markers_tibble <- as_tibble(celltype.markers, rownames='GeneNames')
  celltype.markers_tibble <- arrange(celltype.markers_tibble, avg_log2FC)
  
  celltype.markers.de <- celltype.markers_tibble[celltype.markers_tibble$p_val_adj < pCutoff & celltype.markers_tibble$avg_log2FC < -(FCcutoff) | celltype.markers_tibble$avg_log2FC > FCcutoff, ]
  
  return(celltype.markers.de)
}


bioc_volcano <- function(dataset, celltype, group, ident, ident2, pCutoff=0.05, FCcutoff=0.5) {
  
  ## Create Enhanced Volcano plots from Bioconductor
  ## celltype: cluster (ex. macrophages)
  ## group: how to divide cluster (ex. responders vs. nonresponders)
  ## ident: subcluster to be compared to rest of group (ex. responders)
  ## ident2: rest of group (just used to generate title)
  
  celltype.markers <- FindMarkers(dataset, subset.ident=celltype, group.by=group, ident.1=ident)
  celltype.markers_tibble <- as_tibble(celltype.markers, rownames='GeneNames')
  celltype.markers_tibble <- arrange(celltype.markers_tibble, avg_log2FC)
  
  celltype.volcano <- EnhancedVolcano(celltype.markers_tibble,
                                      lab = celltype.markers_tibble$GeneNames,
                                      x = 'avg_logFC', y = 'p_val_adj',
                                      pCutoff = pCutoff, FCcutoff = FCcutoff,
                                      title = paste(celltype, "DE:", ident, "vs.", ident2))
  
  return(celltype.volcano)
  
}


anova_test <- function(data.frame, gene, category, confidence = 0.95, xlab = "Comparisons", ylab = "Difference in Mean Expression") {

  anova <- aov(gene~category, data=data.frame)
  tukey <- TukeyHSD(anova, conf.level = confidence)
  plot <- tuk_plot(tukey, xlab = xlab, ylab = ylab, title = paste0("ANOVA Test: Mean Expression of ", gene, ": ", confidence, " Confidence Level"))
  
  return(plot)
  
}


tuk_plot <- function (x, xlab="Difference in Mean Expression", ylab="Comparisons", ylabels = NULL, title="ANOVA Test", ...) {
  
  ## Plot results of ANOVA / Tukey Test
  
  for (i in seq_along(x)) {
    xi <- x[[i]][, -4L, drop = FALSE]
    yvals <- nrow(xi):1L
    dev.hold()
    on.exit(dev.flush())
    plot(c(xi[, "lwr"], xi[, "upr"]), rep.int(yvals, 2L), 
         type = "n", axes = FALSE, xlab = "", ylab = "", main = NULL, 
         ...)
    axis(1, ...)
    
    if (is.null(ylabels)) ylabels <- dimnames(xi)[[1L]]
    
    axis(2, at = nrow(xi):1, labels = ylabels, 
         srt = 0, ...)
    abline(h = yvals, lty = 1, lwd = 0.5, col = "lightgray")
    abline(v = 0, lty = 2, lwd = 0.5, ...)
    segments(xi[, "lwr"], yvals, xi[, "upr"], yvals, ...)
    segments(as.vector(xi), rep.int(yvals - 0.1, 3L), as.vector(xi), 
             rep.int(yvals + 0.1, 3L), ...)
    title(main = title,
          xlab = xlab, ylab = ylab)
    
    box()
    dev.flush()
    on.exit()
  }
}

