library(Seurat)
library(dplyr)
setwd("/Volumes/BigHDD/Seurat")
lung <- load("datasets/facs_Lung_seurat_tiss.Robj", verbose = T)
# PCA analysis
tissue <- RunPCA(object = tissue, pc.genes = tissue@var.genes,
                 do.print = TRUE, pcs.print = 1:5, genes.print = 5)
PCAPlot(object = tissue, dim.1 = 1, dim.2 = 2)

# Find significant clusters
PCElbowPlot(object = tissue)
# elbow around PC 18
# Will use clusters 1-18
tissue <- FindClusters(object = tissue, reduction.type = "pca", dims.use = 1:18, 
                     resolution = 0.6, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = tissue)
# do tsne
tissue <- RunTSNE(object = tissue, dims.use = 1:18, do.fast = TRUE)
TSNEPlot(object = tissue)
# find CD4 t cells
VlnPlot(object = tissue, features.plot = c("Cd4", "Foxp3"))
# Cd4 Tcells are cluster 8. It seems that there are 3 Tregs
# Saved FeaturePlot as "Tabula_Muris_Lung_Finding_Rora_Treg.pdf"
FeaturePlot(object = tissue, features.plot = c(
  "Cd3e","Cd8a","Cd4","Foxp3","Rora", "Areg", "Gata3","Irf4", "Melk"), cols.use = c("grey", "red"), 
  reduction.use = "tsne")
# Trying to find subpopulations.
# Take only most variable PCs and increase resolution
tissue <- FindClusters(object = tissue, reduction.type = "pca", dims.use = 1:10, 
                       resolution = 1, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = tissue)
# do tsne again
tissue <- RunTSNE(object = tissue, dims.use = 1:10, do.fast = TRUE)
TSNEPlot(object = tissue)
# Find Tregs
FeaturePlot(object = tissue, features.plot = c(
  "Cd3e","Cd8a","Cd4","Foxp3","Rora", "Areg", "Gata3","Irf4", "Cd44"), cols.use = c("grey", "red"), 
            reduction.use = "tsne")
# Saved FeaturePlot as "Tabula_Muris_Lung_alternative_clust_Finding_Rora_Treg.pdf"
# Alternative find markers for clusters
lung.markers <- FindAllMarkers(object = tissue, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
lung.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# Alternative find markers for clusters
cluster1.markers <- FindMarkers(object = tissue, ident.1 = 1, min.pct = 0.25)
names(cluster.markers) <- paste("cluster.marker", seq(1,15), sep = "")
print(x = head(x = cluster1.markers, n = 5))
for (i in seq(1,15)) {
  cluster.markers[[i]] <- FindMarkers(object = tissue, ident.1 = i, min.pct = 0.25)
}

