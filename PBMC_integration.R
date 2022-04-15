library(Seurat)
library(cowplot)

#Setup the Seurat objects


ctrl.data <- read.table(file = "immune_control_expression_matrix.txt", sep = "\t")
stim.data <- read.table(file = "immune_stimulated_expression_matrix.txt", sep = "\t")

# Set up control object
ctrl <- CreateSeuratObject(counts = ctrl.data, project = "IMMUNE_CTRL", min.cells = 5)
stim <- CreateSeuratObject(counts = stim.data, project = "IMMUNE_STIM", min.cells = 5)


# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
ctrl [["percent.mt"]] <- PercentageFeatureSet(ctrl , pattern = "^MT-")

VlnPlot(ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
stim[["percent.mt"]] <- PercentageFeatureSet(stim , pattern = "^MT-")

VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


ctrl$type <- "CTRL"
stim$type <- "STIM"

ctrl <- subset(ctrl, subset = nFeature_RNA > 250)
stim <- subset(stim, subset = nFeature_RNA > 250)

ctrl_norm <- NormalizeData(ctrl , verbose = FALSE)
stim_norm <- NormalizeData(stim , verbose = FALSE)


ctrl_hvf <- FindVariableFeatures(ctrl_norm, selection.method = "vst", nfeatures = 2000)

stim_hvf <- FindVariableFeatures(stim_norm, selection.method = "vst", nfeatures = 2000)

#Perform integration

#We then identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input, and use these anchors to integrate the two datasets together with IntegrateData.

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl_hvf, stim_hvf), dims = 1:50)
immune.anchors


immune.combined <- IntegrateData(anchorset = immune.anchors)


DefaultAssay(immune.combined) <- "integrated"


# Run the standard workflow for visualization and clustering
immune.combined_scaled <- ScaleData(immune.combined, verbose = FALSE)
immune.combined_pca <- RunPCA(immune.combined_scaled, npcs = 30, verbose = FALSE)


#Run harmony to remove batch effect
immune.combined_b <- RunHarmony(immune.combined_pca, group.by.vars = "stim", assay.use='integrated')

# Run UMAP 
immune.combined_umap <- RunUMAP(immune.combined_b, reduction = "pca", dims = 1:20)


# Run Clustering
immune.combined_neighbour <- FindNeighbors(immune.combined_umap, reduction = "pca", dims = 1:20)
immune.combined_clus <- FindClusters(immune.combined_neighbour, resolution = 0.5)


# Visualization
p1 <- DimPlot(immune.combined_clus, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined_clus, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


DimPlot(immune.combined_clus, reduction = "umap", split.by = "type")

DefaultAssay(immune.combined_clus) <- "RNA"


FeaturePlot(immune.combined_clus, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A", 
                                               "CCL2", "PPBP"), min.cutoff = "q9")



immune.combined_clus1 <- RenameIdents(immune.combined_clus, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T", 
                                      `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "T activated", `7` = "NK", `8` = "DC", `9` = "B Activated", 
                                      `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets")

DimPlot(immune.combined_clus1, label = TRUE)






library(metap)
nk.markers <- FindConservedMarkers(immune.combined_clus, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
head(nk.markers)
