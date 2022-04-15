library(Seurat)
library(cowplot)
library(harmony)
library(patchwork)
library(dplyr)
library(multtest)
library(metap)


ctrl.data <- read.table(file = "immune_control_expression_matrix_copy.txt", sep = "\t")
stim.data <- read.table(file = "immune_stimulated_expression_matrix_copy.txt", sep = "\t")

#Set up stimulated object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)

stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)


#Add column with sample information
ctrl$stim <- "CTRL"
stim$stim <- "STIM"


#mitocondrial 
ctrl[["percent.mt"]] <- PercentageFeatureSet(ctrl, pattern = "^MT-")
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")


#Visualize QC metrics as a violin plot
VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Remove low quality cells
ctrl_filter <- subset(ctrl, subset = nFeature_RNA > 250)
stim_filter <- subset(stim, subset = nFeature_RNA > 250)


#Visualize QC metrics as a violin plot
#VlnPlot(ctrl_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#VlnPlot(stim_filter, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Normalize data
ctrl_norm <- NormalizeData(ctrl_filter, verbose = FALSE)
stim_norm <- NormalizeData(stim_filter, verbose = FALSE)

#Select higher variable features
ctrl_hv <- FindVariableFeatures(ctrl_norm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_ctrl <- head(VariableFeatures(ctrl_hv), 10)

# plot variable features with and without labels
plot1_c <- VariableFeaturePlot(ctrl_hv)
plot1_c
plot2_c <- LabelPoints(plot = plot1_c, points = top10_ctrl, repel = TRUE)
plot2_c



stim_hv <- FindVariableFeatures(stim_norm, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10_stim <- head(VariableFeatures(stim_hv), 10)

# plot variable features with and without labels
plot1_st <- VariableFeaturePlot(stim_hv)
plot1_st
plot2_st <- LabelPoints(plot = plot1_st, points = top10_stim, repel = TRUE)
plot2_st


ctrl_hv$stim <- "CTRL"
stim_hv$stim <- "STIM"

#create anchors
immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl_hv, stim_hv), dims = 1:20)
#Integrate data using anchors
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)
#Set integrated assay as default assay
DefaultAssay(immune.combined) <- "integrated"

#Scale Data
immune.combined_scaled <- ScaleData(immune.combined, verbose = FALSE)
#Perform linear dimensionality reduction PCA
immune.combined_pca <- RunPCA(immune.combined_scaled, npcs = 50, verbose = FALSE)

DimPlot(immune.combined_pca, reduction = "pca")

immune.combined_p <- JackStraw(immune.combined_pca, num.replicate = 100)

immune.combined_score <- ScoreJackStraw(immune.combined_p, dims = 1:20)

JackStrawPlot(immune.combined_score, dims = 1:20)

ElbowPlot(immune.combined_score)

#Remove Batch Effect
#Run harmony or RunFastMNN to remove batch effect
immune.combined_b <- RunHarmony(immune.combined_pca, group.by.vars = "stim", assay.use='integrated')

#Perform non-linear dimensionality reduction -  UMAP
immune.combined_umap <- RunUMAP(immune.combined_b, reduction = "pca", dims = 1:15)

#DimPlot(immune.combined_umap, reduction = "pca", group.by = "stim")

#Clustering: Graph based clustering using shared nearest neighbours here, first it find nearest neighbours
immune.combined_neighbour <- FindNeighbors(immune.combined_umap, reduction = "pca")
immune.combined_clus <- FindClusters(immune.combined_neighbour, resolution = 0.5)
#Visualization
p1 <- DimPlot(immune.combined_clus, reduction = "umap", group.by = "stim")
p1

p2 <- DimPlot(immune.combined_clus, reduction = "umap", label = TRUE)
p2

plot_grid(p1, p2)


DimPlot(immune.combined_clus, reduction = "umap", split.by = "stim")


library(metap)
conserved_markers <- FindConservedMarkers(immune.combined_clus, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(conserved_markers)


data.markers <- FindAllMarkers(immune.combined_clus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.1)

Top50_marker_pbmc <- data.markers %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)

write.table(Top50_marker_pbmc,file="Top25_markers_integrated_samples_clustermole.txt", sep='\t',  quote = F,row.names = FALSE)


top_markers <- data.markers %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC)
head(top25_markers,21) 


VlnPlot(immune.combined_clus, features = c( "CTSB", "SELL", "GPR171", "FCGR3A", "IGJ", "GNLY", "ADTRP", "MARCKSL1", "PF4", "SPARC", "TSPAN13", "GZMA" ))

FeaturePlot(immune.combined_clus, features = c("CTSB", "SELL", "GPR171", "FCGR3A", "IGJ", "GNLY", "ADTRP", "MARCKSL1", "PF4", "SPARC", "TSPAN13", "GZMA"), min.cutoff = "q9")


top25_markers

library(clustermole)

Top50_marker_pbmc <- data.markers %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC)


Type <- c()
clusters <- c()
for (i in seq(0,12, by=1)) {
  markers_filt <- as.data.frame(Top50_marker_pbmc %>% filter(cluster == as.character(i)))
  my_overlaps <- clustermole_overlaps(genes = markers_filt$gene, species = "hs")
  type <- as.character(my_overlaps[1,5])
  Type <- c(Type, type)
  clusters <- c(clusters, i)
}

clust_map_pbmc <- as.data.frame(cbind(clusters,Type))
head(clust_map_pbmc,20)



write.table(clust_map,file="clust_map_cell_types_based_on_clustermole.txt", sep='\t',  quote = F,row.names = FALSE)

new.cluster.ids1 <- c("Monocyte", "IMGN_Tgd_vg2+_Sp", "Naive T Cells", "Monocyte", 
                      "B cell", "IMGN T 8Eff- Sp_OT1_d8", "Naive T Cells",
                     "N-Monocyte","Erythroblast","Platelets", 
                     "Megakaryocyte Progenitor","Dendritic Cell", "Paneth cell")

names(new.cluster.ids1) <- levels(immune.combined_clus)

#Rename clusters with cell types
immune.combined_clus_new <- RenameIdents(immune.combined_clus, new.cluster.ids1)


#type based

DimPlot(immune.combined_clus_new, reduction = "umap", label = TRUE, pt.size = 0.4) + NoLegend()


#####
Idents(immune.combined_clus_new) <- factor(Idents(immune.combined_clus_new), levels = c("Monocyte", "IMGN_Tgd_vg2+_Sp", "Naive T Cells", 
                                                                                  "B cell", "IMGN T 8Eff- Sp_OT1_d8", "N-Monocyte","Erythroblast","Platelets", 
                                                                                  "Megakaryocyte Progenitor","Dendritic Cell", "Paneth cell"))
markers.to.plot <- c("CTSB", "SELL", "GPR171", "FCGR3A", "IGJ", "GNLY", "ADTRP", "MARCKSL1", "PF4", "SPARC", "TSPAN13", "GZMA",  "IL3RA")

#Visualization
DotPlot(immune.combined_clus_new, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8,  split.by = "stim") + RotatedAxis()


#Identify conserved cell type markers
###### Conserved markers regardless of condition

DefaultAssay(immune.combined) <- "RNA"
markers_0 <- FindConservedMarkers(immune.combined_clus, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
head(markers_0)
markers_1 <- FindConservedMarkers(immune.combined_clus, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
head(markers_1)
markers_2 <- FindConservedMarkers(immune.combined_clus, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
head(markers_2)
markers_3 <- FindConservedMarkers(immune.combined_clus, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
head(markers_3)
markers_4 <- FindConservedMarkers(immune.combined_clus, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
head(markers_4)
markers_5 <- FindConservedMarkers(immune.combined_clus, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
head(markers_5)
markers_6 <- FindConservedMarkers(immune.combined_clus, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(markers_6)
markers_7 <- FindConservedMarkers(immune.combined_clus, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(markers_7)
markers_8 <- FindConservedMarkers(immune.combined_clus, ident.1 = 8, grouping.var = "stim", verbose = FALSE)
head(markers_8)
markers_9 <- FindConservedMarkers(immune.combined_clus, ident.1 = 9, grouping.var = "stim", verbose = FALSE)
head(markers_9)
markers_10 <- FindConservedMarkers(immune.combined_clus, ident.1 = 10, grouping.var = "stim", verbose = FALSE)
head(markers_10)
markers_11 <- FindConservedMarkers(immune.combined_clus, ident.1 = 11, grouping.var = "stim", verbose = FALSE)
head(markers_11)
markers_12 <- FindConservedMarkers(immune.combined_clus, ident.1 = 12, grouping.var = "stim", verbose = FALSE)
head(markers_12)

markers_12$cluster <- "12"
write.table(Top50_marker_pbmc,file="Top25_markers_integrated_samples_clustermole.txt", sep='\t',  quote = F,row.names = T)
markers_12

Idents(immune.combined_clus_new)

immune.combined_clus_new$celltype.type <- paste(Idents(immune.combined_clus_new), immune.combined_clus_new$stim, sep = "_")
immune.combined_clus_new$celltype <- Idents(immune.combined_clus_new)
Idents(immune.combined_clus_new) <- "celltype.type"

tail(Idents(immune.combined_clus_new))

b_cell_DGE <- FindMarkers(immune.combined_clus_new, ident.1 = "B cell_STIM", ident.2 = "B cell_CTRL", verbose = FALSE)
head(b_cell_DGE, n = 5)

T_cell_DGE <- FindMarkers(immune.combined_clus_new, ident.1 = "Naive T Cells_STIM", ident.2 = "Naive T Cells_CTRL", verbose = FALSE)
head(T_cell_DGE, n = 5)

IMGN_DGE <- FindMarkers(immune.combined_clus_new, ident.1 = "IMGN_Tgd_vg2+_Sp_STIM", ident.2 = "IMGN_Tgd_vg2+_Sp_CTRL", verbose = FALSE)
head(IMGN_DGE, n = 10)


library(ggplot2)
library(patchwork)
plots <- VlnPlot(immune.combined_clus_new, features = c("QKI", "IL1RN", "NADK", "ARMCX2", "BATF3"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE) 
CombinePlots(plots = plots, ncol = 1)
wrap_plots(plots = plots, ncol = 1)

 VlnPlot(immune.combined_clus_new, features = c("CD38", "CXCL11", "SCARB2",  "DAPP1", "IL1RN" ), split.by = "stim", group.by = "celltype", pt.size = 0, combine = F) 

 
 wrap_plots(plots = plots, ncol = 1)
