install.packages("caTools")
library(caTools)



BiocManager::install('multtest')
library(multtest)


install.packages('Seurat')

library(Seurat)

install.packages('patchwork')

library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)


data <- Read10X(data.dir = "./")

# min.cells = 3 means gene must be present in atleast 3 cells to keep
seurat_ob <- CreateSeuratObject(counts = data, project = "E1_S1", min.cells = 3, min.features = 200)
seurat_ob




# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob@meta.data, 5)


# Visualize QC metrics as a violin plot: N_features: no. of genes expressed/cell and N_count= no. of transcripts/cell, percent.mt = %age of mt transcripts
VlnPlot(seurat_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


seurat_ob_filter <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

seurat_ob_norm <- NormalizeData(seurat_ob_filter, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_ob_norm  <- NormalizeData(seurat_ob_norm )

seurat_ob_var <- FindVariableFeatures(seurat_ob_norm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_ob_var), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seurat_ob_var)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1

plot2

#plot1 + plot2

all.genes <- rownames(seurat_ob_var)
seurat_ob_scaled <- ScaleData(seurat_ob_var, features = all.genes)

seurat_ob_pca <- RunPCA(seurat_ob_scaled, features = VariableFeatures(object = seurat_ob_scaled))

# Examine and visualize PCA results a few different ways
print(seurat_ob_pca[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(seurat_ob_pca, dims = 1:2, reduction = "pca")

DimPlot(seurat_ob_pca, reduction = "pca")

DimHeatmap(seurat_ob_pca, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(seurat_ob_pca, dims = 1:15, cells = 500, balanced = TRUE)


# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
seurat_ob_n <- JackStraw(seurat_ob_pca, num.replicate = 100)
seurat_ob_n1 <- ScoreJackStraw(seurat_ob_n, dims = 1:20)

#jackstraw plot
JackStrawPlot(seurat_ob_n1, dims = 1:15)

#Elbow plot
ElbowPlot(seurat_ob_n1)



seurat_ob_neigbours <- FindNeighbors(seurat_ob_n1, dims = 1:10)
seurat_ob_cluster <- FindClusters(seurat_ob_neigbours, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_ob_cluster), 5)

# UMAP
seurat_ob_umap <- RunUMAP(seurat_ob_cluster, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(seurat_ob_umap, reduction = "umap")

# find all markers of cluster 6
cluster6.markers <- FindMarkers(seurat_ob_umap, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 5)


# find all markers distinguishing cluster 6 from clusters 7 and 1
cluster6_1.markers <- FindMarkers(seurat_ob_umap, ident.1 = 6, ident.2 = c(7, 1), min.pct = 0.25)
head(cluster6_1.markers, n = 5)


# find markers for every cluster compared to all remaining cells, report only the positive
# ones
all_clusters_markers <- FindAllMarkers(seurat_ob_umap, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

all_clusters_markers %>% group_by(cluster) %>%  slice_max(n = 2, order_by = avg_log2FC)

all_clusters_markers <- as.data.frame(all_clusters_markers)

#write into file
tail(all_clusters_markers, 20)


#cluster0.markers <- FindMarkers(seurat_ob_umap, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(seurat_ob_umap, features = c("Usp48",  "Hdgfl2" ))

# you can plot raw counts as well
VlnPlot(seurat_ob_umap, features = c("Usp48",  "Hdgfl2"), slot = "counts", log = TRUE)


all_clusters_markers %>% group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(seurat_ob_umap, features = top5$gene) + NoLegend()



######### types of cells
new.cluster.ids <- c("Cell_1", "Cell_2", "Cell_3", "Cell_4", "Cell_5", "Cell_6", "Cell_7", "Cell_8")
names(new.cluster.ids) <- levels(seurat_ob_umap)
seurat_ob_type<- RenameIdents(seurat_ob_umap, new.cluster.ids)
DimPlot(seurat_ob_type, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

