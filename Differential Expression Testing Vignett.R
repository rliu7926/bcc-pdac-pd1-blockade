################################################################################
#                                                                              #
#                     Seurat - Guided Clustering Tutorial                      #
#                                                                              #
################################################################################

library(Seurat)
pbmc <- readRDS(file = "C:\\Users\\Ryan\\Downloads\\pbmc3k_final.rds")


### Perform Default Differential Expression Tests ###

# List options for groups to perform idfferential expression on
levels(pbmc)

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
# View results
head(monocyte.de.markers)

# Find differentially expressed features between CD+14 monocytes and all other cells
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = NULL, only.pos = TRUE)
# View results
head(monocyte.de.markers)


### Prefilter Features or Cells to Increase the Speed of DE Testing ###

# Pre-filter features that are detected at <50% frequency in either CD14+ or FCGR3A+ monocytes
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.pct = 0.5))

# Pre-filter features that have less than a two-fold change between the average expression fo CD14+ monocytes vs. FCGR3A+ monocytes
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))

# Pre-filter features whose detection percentages across the two groups are similar (within 0.25) 
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.diff.pct = 0.25))

#Subsample each group to a maximum of 200 cells
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", max.cells.per.ident = 200))


### Perform DE Analysis Using Alternative Tests ###

# Test for DE features using the MAST package

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Mast")

head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "MAST"))

# Test for DE features using the DESeq2 package

BiocManager::install("DESeq2")

head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "DESeq2", max.cells.per.ident = 50))
