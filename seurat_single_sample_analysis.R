library(Seurat)
library(cowplot)
library(harmony)
library(patchwork)
library(dplyr)
library(multtest)
library(metap)


data <- Read10X(data.dir = "pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")

seurat_ob <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

#QC
seurat_ob[["percent.mt"]] <- PercentageFeatureSet(seurat_ob, pattern = "^MT-")

#Visualize QC metrics as a violin plot
VlnPlot(seurat_ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

data_filter <- subset(seurat_ob, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

VlnPlot(data_filter , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


data_norm <- NormalizeData(data_filter, normalization.method = "LogNormalize", scale.factor = 10000)

all.genes <- rownames(data_norm)
data_scaled <- ScaleData(data_norm, features = all.genes)

data_hv <- FindVariableFeatures(data_scaled, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(data_hv), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(data_hv)
plot1
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2


data_pca <- RunPCA(data_hv)


DimPlot(data_pca, reduction = "pca")

DimHeatmap(data_pca, dims = 1:15, cells = 500, balanced = TRUE)


data_p <- JackStraw(data_pca, num.replicate = 100)
data_score <- ScoreJackStraw(data_p, dims = 1:20)
JackStrawPlot(data_score, dims = 1:20)


ElbowPlot(data_score)

data_umap <- RunUMAP(data_pca, dims = 1:10)


#find neighbours
data_nn <- FindNeighbors(data_umap, dims = 1:10) ## taken first 10 pca components

#find clusters
data_cluster <- FindClusters(data_nn, resolution = 0.5)

#Visualize clusters
DimPlot(data_cluster, reduction = "umap")


data.markers <- FindAllMarkers(data_cluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Filter top 25 markers
Top25_markers <- data.markers %>% group_by(cluster) %>% slice_max(n = 25, order_by = avg_log2FC)

VlnPlot(data_cluster, features = c("CCR7", "LYZ" , "AQP3", "CD79A" , "GZMK", "LST1", "GNLY", "HLA-DPB1", "PF4"))

FeaturePlot(data_cluster, features = c("CCR7", "LYZ" , "AQP3", "CD79A" , "GZMK", "LST1", "GNLY", "HLA-DPB1", "PF4"))

top25_markers

library(clustermole)

Type <- c()
clusters <- c()
for (i in seq(0,8, by=1)) {
  markers_filt <- as.data.frame(Top25_markers %>% filter(cluster == as.character(i)))
  my_overlaps <- clustermole_overlaps(genes = markers_filt$gene, species = "hs")
  type <- as.character(my_overlaps[1,5])
  Type <- c(Type, type)
  clusters <- c(clusters, i)
}


clust_map_pbmc <- as.data.frame(cbind(clusters,Type))
head(clust_map_pbmc,20)



write.table(clust_map,file="clust_map_cell_types_based_on_clustermole.txt", sep='\t',  quote = F,row.names = FALSE)

new.cluster.ids <- c("Naive T cells", "Neuro_ep Dendritic Cells", "Naive T cells","B cell",
                     "NK-NKT cells","Neuro_ep Dendritic cells","NK Cells", "BM-Dendritic Cells","Platelets")
names(new.cluster.ids) <- levels(data_cluster)

#Rename clusters with cell types
data_cluster_new <- RenameIdents(data_cluster, new.cluster.ids)


#type based

DimPlot(data_cluster_new, reduction = "umap", label = TRUE, pt.size = 0.7) + NoLegend()
