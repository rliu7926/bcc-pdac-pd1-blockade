bioc_markers <- function(dataset, celltype, group, ident, pCutoff=0.05, FCcutoff=0.5) {
  
  ## Find differentially expressed markers for Enhanced Volcano plot
  
  celltype.markers <- FindMarkers(dataset, subset.ident=celltype, group.by=group, ident.1=ident)
  celltype.markers_tibble <- as_tibble(celltype.markers, rownames='GeneNames')
  celltype.markers_tibble <- arrange(celltype.markers_tibble, avg_log2FC)
  
  celltype.markers.de <- celltype.markers_tibble[celltype.markers_tibble$p_val_adj < pCutoff & celltype.markers_tibble$avg_log2FC < -(FCcutoff) | celltype.markers_tibble$avg_log2FC > FCcutoff, ]
  
  return(celltype.markers.de)
}


bioc_volcano <- function(dataset, celltype, group, ident, ident2, pCutoff=0.05, FCcutoff=0.5, title) {
  
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
                                      title = paste(title, "DE:", ident, "vs.", ident2))
  
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

fourway_diffexpression_MHCandHSP <- function(celltype, celltypename) {
  
  Idents(BCC.responder.pre) <- BCC.responder.pre$seurat_clusters
  Idents(BCC.responder.post) <- BCC.responder.post$seurat_clusters
  
  Idents(BCC.nonresponder.pre) <- BCC.nonresponder.pre$seurat_clusters
  Idents(BCC.nonresponder.post) <- BCC.nonresponder.post$seurat_clusters
  
  BCC.celltype.responder.pre <- subset(BCC.responder.pre, idents=celltype)
  BCC.celltype.responder.post <- subset(BCC.responder.post, idents=celltype)
  
  BCC.celltype.nonresponder.pre <- subset(BCC.nonresponder.pre, idents=celltype)
  BCC.celltype.nonresponder.post <- subset(BCC.nonresponder.post, idents=celltype)
  
  ## Set idents to analyze by patient
  
  Idents(BCC.celltype.responder.pre) <- BCC.celltype.responder.pre$patient
  Idents(BCC.celltype.responder.post) <- BCC.celltype.responder.post$patient
  
  Idents(BCC.celltype.nonresponder.pre) <- BCC.celltype.nonresponder.pre$patient
  Idents(BCC.celltype.nonresponder.post) <- BCC.celltype.nonresponder.post$patient
  
  ## Determine average expression of genes
  
  BCC.celltype.responder.pre.averages <- AverageExpression(BCC.celltype.responder.pre)
  BCC.celltype.responder.post.averages <- AverageExpression(BCC.celltype.responder.post)
  BCC.celltype.nonresponder.pre.averages <- AverageExpression(BCC.celltype.nonresponder.pre)
  BCC.celltype.nonresponder.post.averages <- AverageExpression(BCC.celltype.nonresponder.post)
  
  BCC.celltype.responder.pre.averages <- BCC.celltype.responder.pre.averages$RNA
  BCC.celltype.responder.post.averages <- BCC.celltype.responder.post.averages$RNA
  BCC.celltype.nonresponder.pre.averages <- BCC.celltype.nonresponder.pre.averages$RNA
  BCC.celltype.nonresponder.post.averages <- BCC.celltype.nonresponder.post.averages$RNA
  
  ## Filter out important genes
  
  MHC1.genes <- c('HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G')
  MHC2.genes <- c('HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQA2', 'HLA-DQB1', 'HLA-DQB1-AS1', 'HLA-DQB2', 'HLA-DRA', 'HLA-DRB1', 'HLA-DRB5')
  HSP.genes <- c('HSP90AA1', 'HSP90AB1', 'HSP90B1', 'HPSA12A', 'HSPA12B', 'HSPA13', 'HSPA14', 'HPSA1A', 'HSPA1B', 'HSPA1L', 'HSPA2', 'HSPA4', 'HSPA4L', 'HSPA5', 'HSPA6', 'HSPA8', 'HSPA9', 'HSPB1', 'HSPB11', 'HSPB2', 'HSPB3', 'HSPB6', 'HSPB7', 'HSPB8', 'HSPB9', 'HSPBAP1', 'HSPBP1', 'HSPD1', 'HSPE1', 'HSPE1-MOB4', 'HSPG2', 'HSPH1')
  
  BCC.celltype.responder.pre.averages.mhc1 <- BCC.celltype.responder.pre.averages[MHC1.genes,]
  BCC.celltype.responder.pre.averages.mhc2 <- BCC.celltype.responder.pre.averages[MHC2.genes,]
  BCC.celltype.responder.pre.averages.hsp <- BCC.celltype.responder.pre.averages[HSP.genes,]
  
  BCC.celltype.responder.post.averages.mhc1 <- BCC.celltype.responder.post.averages[MHC1.genes,]
  BCC.celltype.responder.post.averages.mhc2 <- BCC.celltype.responder.post.averages[MHC2.genes,]
  BCC.celltype.responder.post.averages.hsp <- BCC.celltype.responder.post.averages[HSP.genes,]
  
  BCC.celltype.nonresponder.pre.averages.mhc1 <- BCC.celltype.nonresponder.pre.averages[MHC1.genes,]
  BCC.celltype.nonresponder.pre.averages.mhc2 <- BCC.celltype.nonresponder.pre.averages[MHC2.genes,]
  BCC.celltype.nonresponder.pre.averages.hsp <- BCC.celltype.nonresponder.pre.averages[HSP.genes,]
  
  BCC.celltype.nonresponder.post.averages.mhc1 <- BCC.celltype.nonresponder.post.averages[MHC1.genes,]
  BCC.celltype.nonresponder.post.averages.mhc2 <- BCC.celltype.nonresponder.post.averages[MHC2.genes,]
  BCC.celltype.nonresponder.post.averages.hsp <- BCC.celltype.nonresponder.post.averages[HSP.genes,]
  
  ## Calculate scores (arithmetic mean)
  
  BCC.celltype.responder.pre.averages.mhc1$gene <- NULL
  BCC.celltype.responder.pre.averages.mhc2$gene <- NULL
  BCC.celltype.responder.pre.averages.hsp$gene <- NULL
  
  BCC.celltype.responder.post.averages.mhc1$gene <- NULL
  BCC.celltype.responder.post.averages.mhc2$gene <- NULL
  BCC.celltype.responder.post.averages.hsp$gene <- NULL
  
  BCC.celltype.nonresponder.pre.averages.mhc1$gene <- NULL
  BCC.celltype.nonresponder.pre.averages.mhc2$gene <- NULL
  BCC.celltype.nonresponder.pre.averages.hsp$gene <- NULL
  
  BCC.celltype.nonresponder.post.averages.mhc1$gene <- NULL
  BCC.celltype.nonresponder.post.averages.mhc2$gene <- NULL
  BCC.celltype.nonresponder.post.averages.hsp$gene <- NULL
  
  BCC.celltype.responder.pre.averages.mhc1.score <- as.data.frame(colMeans(BCC.celltype.responder.pre.averages.mhc1, na.rm=TRUE))
  BCC.celltype.responder.pre.averages.mhc2.score <- as.data.frame(colMeans(BCC.celltype.responder.pre.averages.mhc2, na.rm=TRUE))
  BCC.celltype.responder.pre.averages.hsp.score <- as.data.frame(colMeans(BCC.celltype.responder.pre.averages.hsp, na.rm=TRUE))
  
  BCC.celltype.responder.post.averages.mhc1.score <- as.data.frame(colMeans(BCC.celltype.responder.post.averages.mhc1, na.rm=TRUE))
  BCC.celltype.responder.post.averages.mhc2.score <- as.data.frame(colMeans(BCC.celltype.responder.post.averages.mhc2, na.rm=TRUE))
  BCC.celltype.responder.post.averages.hsp.score <- as.data.frame(colMeans(BCC.celltype.responder.post.averages.hsp, na.rm=TRUE))
  
  BCC.celltype.nonresponder.pre.averages.mhc1.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.pre.averages.mhc1, na.rm=TRUE))
  BCC.celltype.nonresponder.pre.averages.mhc2.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.pre.averages.mhc2, na.rm=TRUE))
  BCC.celltype.nonresponder.pre.averages.hsp.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.pre.averages.hsp, na.rm=TRUE))
  
  BCC.celltype.nonresponder.post.averages.mhc1.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.post.averages.mhc1, na.rm=TRUE))
  BCC.celltype.nonresponder.post.averages.mhc2.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.post.averages.mhc2, na.rm=TRUE))
  BCC.celltype.nonresponder.post.averages.hsp.score <- as.data.frame(colMeans(BCC.celltype.nonresponder.post.averages.hsp, na.rm=TRUE))
  
  ## Create dataframe to store all scores
  
  x <- nrow(BCC.celltype.responder.pre.averages.mhc1.score)
  y <- nrow(BCC.celltype.responder.post.averages.mhc1.score)
  z <- nrow(BCC.celltype.nonresponder.pre.averages.mhc1.score)
  w <- nrow(BCC.celltype.nonresponder.post.averages.mhc1.score)
  
  BCC.celltype.averages.score <- data.frame()
  
  BCC.celltype.averages.score[1:x,1] <- BCC.celltype.responder.pre.averages.mhc1.score
  BCC.celltype.averages.score[1:x,2] <- BCC.celltype.responder.pre.averages.mhc2.score
  BCC.celltype.averages.score[1:x,3] <- BCC.celltype.responder.pre.averages.hsp.score
  BCC.celltype.averages.score[(x+1):(x+y),1] <- BCC.celltype.responder.post.averages.mhc1.score
  BCC.celltype.averages.score[(x+1):(x+y),2] <- BCC.celltype.responder.post.averages.mhc2.score
  BCC.celltype.averages.score[(x+1):(x+y),3] <- BCC.celltype.responder.post.averages.hsp.score
  
  BCC.celltype.averages.score[(x+y+1):(x+y+z),1] <- BCC.celltype.nonresponder.pre.averages.mhc1.score
  BCC.celltype.averages.score[(x+y+1):(x+y+z),2] <- BCC.celltype.nonresponder.pre.averages.mhc2.score
  BCC.celltype.averages.score[(x+y+1):(x+y+z),3] <- BCC.celltype.nonresponder.pre.averages.hsp.score
  BCC.celltype.averages.score[(x+y+z+1):(x+y+z+w),1] <- BCC.celltype.nonresponder.post.averages.mhc1.score
  BCC.celltype.averages.score[(x+y+z+1):(x+y+z+w),2] <- BCC.celltype.nonresponder.post.averages.mhc2.score
  BCC.celltype.averages.score[(x+y+z+1):(x+y+z+w),3] <- BCC.celltype.nonresponder.post.averages.hsp.score
  
  colnames(BCC.celltype.averages.score) <- c("MHC1", "MHC2", "HSP")
  
  
  BCC.celltype.averages.score$Response[1:(x+y+z+w)] <- c("placeholder")
  
  BCC.celltype.averages.score$Response[1:x] <- c("Responsive.Pre")
  BCC.celltype.averages.score$Response[(x+1):(x+y)] <- c("Responsive.Post")
  BCC.celltype.averages.score$Response[(x+y+1):(x+y+z)] <- c("Nonresponsive.Pre")
  BCC.celltype.averages.score$Response[(x+y+z+1):(x+y+z+w)] <- c("Nonresponsive.Post")
  
  
  BCC.celltype.averages.score$Patient <- rownames(BCC.celltype.averages.score)
  
  BCC.celltype.averages.score.melt <- melt(BCC.celltype.averages.score, id=c("Response", "Patient"))
  
  BCC.celltype.averages.score.melt.MHC1 <- BCC.celltype.averages.score.melt[c(1:(x+y+z+w)),]
  BCC.celltype.averages.score.melt.MHC2 <- BCC.celltype.averages.score.melt[c((x+y+z+w+1):(2*(x+y+z+w))),]
  BCC.celltype.averages.score.melt.HSP <- BCC.celltype.averages.score.melt[c((2*(x+y+z+w)+1):(3*(x+y+z+w))),]
  
  ## Plot scores
  
  comparisons <- list(c("Responsive.Pre", "Responsive.Post"), c("Nonresponsive.Pre", "Nonresponsive.Post"), c("Responsive.Pre", "Nonresponsive.Pre"), c("Responsive.Post", "Nonresponsive.Post"))
  
  BCC.celltype.averages.score.graph.MHC1 <- ggboxplot(BCC.celltype.averages.score.melt.MHC1, x="Response", y="value", color="Response", title = "MHC Class I", add = "jitter") + stat_compare_means(comparisons = comparisons)
  BCC.celltype.averages.score.graph.MHC1 <- ggpar(p=BCC.celltype.averages.score.graph.MHC1, font.main=c(16, "bold", "black")) + theme(plot.title = element_text(hjust = 0.5))
  
  BCC.celltype.averages.score.graph.MHC2 <- ggboxplot(BCC.celltype.averages.score.melt.MHC2, x="Response", y="value", color="Response", title = "MHC Class II", add = "jitter") + stat_compare_means(comparisons = comparisons) 
  BCC.celltype.averages.score.graph.MHC2 <- ggpar(p=BCC.celltype.averages.score.graph.MHC2, font.main=c(16, "bold", "black")) + theme(plot.title = element_text(hjust = 0.5))
  
  BCC.celltype.averages.score.graph.HSP <- ggboxplot(BCC.celltype.averages.score.melt.HSP, x="Response", y="value", color="Response", title = "HSP Genes", add = "jitter") + stat_compare_means(comparisons = comparisons)
  BCC.celltype.averages.score.graph.HSP <- ggpar(p=BCC.celltype.averages.score.graph.HSP, font.main=c(16, "bold", "black")) + theme(plot.title = element_text(hjust = 0.5))
  
  graph <- grid.arrange(BCC.celltype.averages.score.graph.MHC1,
                        BCC.celltype.averages.score.graph.MHC2,
                        BCC.celltype.averages.score.graph.HSP,
                        nrow = 1,
                        top = textGrob(paste0("Differential Expression of MHC and HSP Genes in the ", celltypename, " of BCC Patients"), gp=gpar(fontsize=30)))
  
  return(graph)
  
}

BCC_diffexp_dumbbell <- function(df1, df2, filterx=0.5, filtery=0.05, title) {
  
  df1.data <- df1$data %>% select(GeneNames, xvals, yvals)
  df2.data <- df2$data %>% select(GeneNames, xvals, yvals)
  
  filterxneg <- -1 * filterx
  
  df1.data <- filter(df1.data, xvals < filterxneg | xvals > filterx, yvals < filtery)
  df2.data <- filter(df2.data, xvals < filterxneg | xvals > filterx, yvals < filtery)
  
  df1.data$xvals[df1.data$xvals < -2] <- -2
  df1.data$xvals[df1.data$xvals > 2] <- 2
  
  df2.data$xvals[df2.data$xvals < -2] <- -2
  df2.data$xvals[df2.data$xvals > 2] <- 2
  
  df1.data$Treatment <- "Pre-treatment"
  df2.data$Treatment <- "Post-treatment"
  
  data <- rbind(df1.data, df2.data)
  data <- data %>% group_by(GeneNames) %>% filter(n()>1)
  
  data <- data[order(data$GeneNames),]
  data$paired <- 0
  for (i in 1:nrow(data)) {
    data$paired[i] <- ceiling(i/2)
  }
  
  data$yvals <- -log10(data$yvals)
  
  graph <- data %>%
    ggplot(aes(x = xvals, y = GeneNames)) +
    labs(title = title,
         x = expression("Log"[2]*" Fold Change"),
         y = "Gene") +
    geom_line(aes(group = paired), color="black", size=1) + 
    geom_point(aes(color = Treatment, size = yvals)) +
    guides(size = FALSE) +
    theme_light() +
    theme(legend.position="top", 
          plot.title = element_text(color="black", size=24, face="bold", hjust = 0.5)) +
    xlim(-2,2)
  
  return(graph)
 
}

BCC_score_dumbbell <- function(df, patient, xlim=50) {
  
  patient.id.char <- BCC.patients.char[patient]
  
  if (df == BCC.MHC1.scores.na) {
    title <- c("MHC-1")
  } else if (df == BCC.MHC2.scores.na) {
    title <- c("MHC-2")
  } else {
    title <- c("HSP")
  }
  
  colnames(df)[patient] <- "xvals"
  
  graph <- df %>%
    ggplot(aes(x = xvals, y = Cluster)) +
    labs(title = paste0(title, " Gene Score: ", patients[patient]), x = "Score", y = "Cluster") +
    geom_line(aes(group = paired), color="black", size=1) + 
    geom_point(aes(color = Treatment, size = 3)) +
    guides(size = FALSE) +
    theme_light() +
    theme(legend.position="top", 
          plot.title = element_text(color="black", size=18, face="bold", hjust = 0.5))
  
  colnames(df)[patient] <- BCC.patients.char[patient]
  
  return(graph)
  
}
