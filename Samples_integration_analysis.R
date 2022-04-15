#
install.packages('cowplot')

#install.packages("harmony")


# Install the remotes package
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)

remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

install.packages("patchwork")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multtest", force = TRUE)

install.packages("metap")


library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

library(multtest)
library(metap)
library(harmony)


data1 <- Read10X(data.dir = "E1_S1/filtered_feature_bc_matrix/")

seurat_ob_S1 <- CreateSeuratObject(counts = data1, project = "E1_S1", min.cells = 3, min.features = 200)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S1[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S1, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob_S1@meta.data, 5)

jpeg('QC_sample1.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

data4 <- Read10X(data.dir = "E4_S4/filtered_feature_bc_matrix/")

seurat_ob_S4 <- CreateSeuratObject(counts = data4, project = "E4_S4", min.cells = 3, min.features = 200)


head(seurat_ob_S4@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S4[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S4, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob_S4@meta.data, 5)

jpeg('QC_Encap.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

seurat_ob_S1$type <- "EB_Ctrl"
seurat_ob_S4$type <- "Encap_Ctrl"


# filter genes with low values (Here, we are removing empty cells/droplet)
EB <- subset(seurat_ob_S1, subset = nFeature_RNA > 250)
Encap <- subset(seurat_ob_S4, subset = nFeature_RNA > 250)


#normalize data

EB_norm<- NormalizeData(EB, verbose = FALSE)
Encap_norm<- NormalizeData(Encap, verbose = FALSE)
# Show QC metrics for the first 5 cells
head(EB_norm@meta.data, 5)

#select highly variable features from each data
EB_ftr<- FindVariableFeatures(EB_norm, selection.method = "vst", nfeatures = 5000)
Encap_ftr<- FindVariableFeatures(Encap_norm, selection.method = "vst", nfeatures = 5000)


# Show QC metrics for the first 5 cells
head(EB_ftr@meta.data, 5)


## integrate samples
integ_anchors <- FindIntegrationAnchors(object.list = list(EB_ftr, Encap_ftr), dims = 1:20)
integ_anchors

Integ_data <- IntegrateData(anchorset = integ_anchors, dims = 1:20)

DefaultAssay(Integ_data) <- "integrated"


# Run the standard workflow for visualization and clustering
Integ_data_scaled <- ScaleData(Integ_data, verbose = FALSE)
Integ_data_pca <- RunPCA(Integ_data_scaled, npcs = 30, verbose = FALSE)

#run harmony or RunFastMNN to remove batch effect

Integ_data_pca_b <- RunHarmony(Integ_data_pca, group.by.vars = "type", assay.use='integrated')

#Integ_data_pca_b <- RunFastMNN(object.list = SplitObject(Integ_data_pca, split.by = "type"))
head(Integ_data_pca_b@meta.data, 5)
Integ_data_pca_b

# t-SNE/UMAPP and Clustering
Integ_data_umap <- RunUMAP(Integ_data_pca_b, reduction = "pca", dims = 1:20)


#Graph based clustering using shared nearest neighbours 
#here, first it find nearest neighbours
Integ_data_neighbour <- FindNeighbors(Integ_data_umap, reduction = "pca", dims = 1:20)

#clustering
Integ_data_clus <- FindClusters(Integ_data_neighbour, resolution = 0.5)


# Visualization
p1 <- DimPlot(Integ_data_clus, reduction = "umap", group.by = "type")
p2 <- DimPlot(Integ_data_clus, reduction = "umap", label = TRUE)


jpeg('UMAP_plot_comb1.jpg', units="in", width=10, height=8, res=300)
plot_grid(p1, p2)

dev.off()

#visualize both samples side by side

jpeg('UMAP_plot_comb2_side_by_side.jpg', units="in", width=10, height=8, res=300)
DimPlot(Integ_data_clus, reduction = "umap", split.by = "type")
dev.off()

DefaultAssay(Integ_data_clus) <- "RNA"

#find conserved markers: calculated the genes that are conserved markers irrespective of stimulation condition in cluster 0 (large cluster).
clus0_markers <- FindConservedMarkers(Integ_data_clus, ident.1 = 0, grouping.var = "type", verbose = FALSE)
head(clus0_markers, 15)

dim(clus0_markers)

clus0_markers_df <- as.data.frame(clus0_markers)
#Write results into a file
write.table(clus0_markers_df,file="cluster0_markers.txt", sep='\t',  quote = F,row.names = T)


######### 

jpeg('Marker_feature_plot_ED.jpg', units="in", width=10, height=8, res=300)

#We can explore marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(Integ_data_clus, features = c("Hspd1", "Ncl", "Supt16", "Set", "Nasp","Hnrnpf", "mt-Co1", "Gm49359", "mt-Cytb",  "Gm42418", "Gm47283", "Gm26917"), min.cutoff = "q9")
dev.off()

Integ_data_clus1 <- RenameIdents(Integ_data_clus, `0` = "cell_type1", `1` = "cell_type2", `2` = "cell_type3", 
                                 `3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9")

jpeg('UMAP_plot_with_diff_cell_types.jpg', units="in", width=10, height=8, res=300)

DimPlot(Integ_data_clus1, label = TRUE)

dev.off()

jpeg('UMAP_plot_both_samples.jpg', units="in", width=10, height=8, res=300)

DimPlot(Integ_data_clus1, label = TRUE, group.by = "type")
dev.off()

# Visualization
pp1 <- DimPlot(Integ_data_clus1, reduction = "umap", group.by = "type")
pp2 <- DimPlot(Integ_data_clus1, reduction = "umap", label = TRUE)

jpeg('UMAP_plot_both_samples_side_by_side.jpg', units="in", width=10, height=8, res=300)
plot_grid(pp1, pp2)

dev.off()


####### Based on ED ##########


#We can explore marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(Integ_data_clus, features = c("Oct4", "Foxa2", "Sox17", "Nanog", "Bmp2","TH", "Nirr2", "Otx2", "Pax6",  "SMest", "GAD65", "Gata3", "Sox21", "Ng2", "GFAP", "Krt3", "Sox10", "Nestin"))


#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
#FeaturePlot(Integ_data_clus, features = c("mt-Nd3", "Rpl41", "Rps20", "Rpl37", "Hsp90aa1", "mt-Nd1",  "Sox10"), min.cutoff = "q9")


Integ_data_clus1 <- RenameIdents(Integ_data_clus, `0` = "cell_type1", `1` = "cell_type2", `2` = "cell_type3", 
`3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9")

DimPlot(Integ_data_clus1, label = TRUE)
DimPlot(Integ_data_clus1, label = TRUE, group.by = "type")

# Visualization
pp1 <- DimPlot(Integ_data_clus1, reduction = "umap", group.by = "type")
pp2 <- DimPlot(Integ_data_clus1, reduction = "umap", label = TRUE)
plot_grid(pp1, pp2)


#assign clusters with cell type
Integ_data_clus1 <- RenameIdents(Integ_data_clus, `0` = "cell_type1", `1` = "cell_type2", `2` = "cell_type3", 
                                 `3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9")


Idents(Integ_data_clus1) <- factor(Idents(Integ_data_clus1), levels = c("cell_type1", "cell_type2", 
                                                                        "cell_type3", "cell_type4", "cell_type5", "cell_type6", "cell_type7", "cell_type8", "cell_type9"))
markers.to.plot <- c("Oct4", "Foxa2", "Sox17", "Nanog", "Bmp2","TH", "Nirr2", "Otx2", "Pax6",  "SMest", "GAD65", "Gata3", "Sox21", "Ng2", "GFAP", "Krt3", "Sox10", "Nestin")

jpeg('Dotplot_genes_cells_ED.jpg', units="in", width=10, height=8, res=300)

DotPlot(Integ_data_clus1, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "type") + RotatedAxis()
dev.off()

################ based on markers avialble in CellMarkers ######


#We can explore marker genes for each cluster and use them to annotate our clusters as specific cell types.

PDF('plot_genes_marker_genes.pdf', units="in", width=20, height=16, res=300)

FeaturePlot(Integ_data_clus, features = c(
  "Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4", ##ESC marker
  "CD146", "CD45", "PDGFR-beta", ##stem markers
  "Ng2", "GFAP", "S100", "Iba1", "Aif1", "CD68", "Cx3cr1", "Itgam", "Aif1", "CD68", "Cx3cr1", "Itgam", "Dock2", "BLBP", "ACSA-2", "CD11b", ##Glial markers
  "Sox10", "NeuN", "DCX", "GFAP", "Iba1", "MAP2", "NF", "MAP2", "Nestin", "Meg3", "PEP-19", ##Neuron markers
  "Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1", ##endodem
  "Mest"," Bmp2"," Eomes"," Hand1"," Isl1"," Kdr"," Mesdc1"," Mesdc2"," Myf5"," Myod1"," Nkx2-5"," T"," Tbx2"," Prdm1"," Tbx1", ##mesoderm
  "Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6", ##ectpderm
  "CD24", "Sca-1", "CD119", "CD13", "CD1d", "CD282", "CD44", "CD49b" ##epithelial
   ))
dev.off()

#We can explore these marker genes for each cluster and use them to annotate our clusters as specific cell types.
#FeaturePlot(Integ_data_clus, features = c("mt-Nd3", "Rpl41", "Rps20", "Rpl37", "Hsp90aa1", "mt-Nd1",  "Sox10"), min.cutoff = "q9")


##Asign cell types ###
Integ_data_clus1 <- RenameIdents(Integ_data_clus, `0` = "cell_type1", `1` = "cell_type2", `2` = "cell_type3", 
                                 `3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9")

DimPlot(Integ_data_clus1, label = TRUE)
DimPlot(Integ_data_clus1, label = TRUE, group.by = "type")

# Visualization
pp1 <- DimPlot(Integ_data_clus1, reduction = "umap", group.by = "type")
pp2 <- DimPlot(Integ_data_clus1, reduction = "umap", label = TRUE)
plot_grid(pp1, pp2)

# Dotplot 

Idents(Integ_data_clus1) <- factor(Idents(Integ_data_clus1), levels = c("cell_type1", "cell_type2", 
                                                                        "cell_type3", "cell_type4", "cell_type5", "cell_type6", "cell_type7", "cell_type8", "cell_type9"))
#markers.to.plot <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4") ##Esc
#markers.to.plot <- c("CD146", "CD45", "PDGFR-beta") ##stem markers
#markers.to.plot <- c("Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6") ##Ectoderm
#markers.to.plot <- c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1") ##Endoderm
#markers.to.plot <- c("Mest"," Bmp2"," Eomes"," Hand1"," Isl1"," Kdr"," Mesdc1"," Mesdc2"," Myf5"," Myod1"," Nkx2-5"," T"," Tbx2"," Prdm1"," Tbx1") ##mesoderm

#markers.to.plot <- c("CD24", "Sca-1", "CD119", "CD13", "CD1d", "CD282", "CD44", "CD49b") ##epithelial

#markers.to.plot <- c("Sox10", "NeuN", "DCX", "GFAP", "Iba1", "MAP2", "NF", "MAP2", "Nestin", "Meg3", "PEP-19") ##Neuron markers
#markers.to.plot <- c("Ng2", "GFAP", "S100", "Iba1", "Aif1", "CD68", "Cx3cr1", "Itgam", "Aif1", "CD68", "Cx3cr1", "Itgam", "Dock2", "BLBP", "ACSA-2", "CD11b") ##Glial markers

#markers.to.plot <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4", "Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6", "Mest", "Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1", "Sox10")

#markers.to.plot <- c("Sox10","NeuN","DCX","GFAP","Iba1","MAP2","NF","Nestin"," Meg3"," PEP-19")

#markers.to.plot <- c("Adcy3"," Ano2"," Cnga2"," Cnga4"," Cngb1"," Gnal"," Ncam1"," Neurog1"," Omp")

jpeg('Dotplot_ESC_Ecto_mes_endo_neuron_markers.jpg', units="in", width=10, height=8, res=300)

DotPlot(Integ_data_clus1, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "type") + RotatedAxis()

dev.off()




###########################################################################

#find  markers for a cluster

clus8_markers <- FindMarkers(Integ_data_clus, ident.1 = 8, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
clus8_markers_df <- as.data.frame(clus8_markers)
#Write results into a file
write.table(clus8_markers_df,file="cluster8_markers.txt", sep='\t',  quote = F,row.names = FALSE)
head(clus8_markers, 15)

Idents(Integ_data_clus1) <- factor(Idents(Integ_data_clus1), levels = c("cell_type1", "cell_type2", 
"cell_type3", "cell_type4", "cell_type5", "cell_type6", "cell_type7", "cell_type8", "cell_type9", "cell_type10"))
markers.to.plot <- c("Rpl13", "Rpl8", "Rps18", "Eef1a1", "Rpl12", "Rpsa", "Tpt1", "Rpl18a", "Rps12", "Rpl32",
                     "Rps19", "Rps23", "Rplp0", "Rplp1", "Rps5")
DotPlot(Integ_data_clus1, features = rev(markers.to.plot), cols = c("blue", "red"), dot.scale = 8, 
        split.by = "type") + RotatedAxis()


#########

Endodem <- subset(Integ_data_clus1, idents = "cell_type9")
Idents(Endodem) <- "type"
avg_cells_9 <- log1p(AverageExpression(Endodem, verbose = FALSE)$RNA)
avg_cells_9$gene <- rownames(avg_cells_9 )

cells_10 <- subset(Integ_data_clus1, idents = "cell_type10")
Idents(cells_10) <- "type"
avg_cells_10 <- log1p(AverageExpression(cells_10, verbose = FALSE)$RNA)
avg_cells_10$gene <- rownames(avg_cells_10)


library(ggplot2)

#genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
#p_11 <- ggplot(avg_cells_9, aes(EB_Ctrl, Encap_Ctrl)) + geom_point() + ggtitle("CD4 Naive T Cells")
#p11 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
#p_22 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
#p22 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
#plot_grid(p_11, p_22)


Integ_data_clus1$celltype.type <- paste(Idents(Integ_data_clus1), Integ_data_clus1$type, sep = "_")
Integ_data_clus1$celltype <- Idents(Integ_data_clus1)
Idents(Integ_data_clus1) <- "celltype.type"
tail(Integ_data_clus1@meta.data, 5)

markers_endo <- FindMarkers(Integ_data_clus1, ident.1 = "cell_type9_EB_Ctrl", ident.2 = "cell_type9_Encap_Ctrl", verbose = FALSE)
head(markers_endo, n = 10)
dim(markers_endo)

write.table(markers_endo,file="cell_type1_DGE.txt", sep='\t',  quote = F,row.names = T)
head(clus8_markers, 15)


FeaturePlot(Integ_data_clus1, features = c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1", "Cib1", "Vcp", "Tmem59", "Idh2"), split.by = "type", max.cutoff = 3, 
            cols = c("grey", "red"))



jpeg('violin_plot_sig_genes_cell_type9.jpg', units="in", width=10, height=12, res=300)

plots <- VlnPlot(Integ_data_clus1, features = c( "Cib1", "Vcp", "Tmem59", "Idh2"), split.by = "type", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

dev.off()


jpeg('violin_plot_common_markers_cell_type9.jpg', units="in", width=10, height=12, res=300)

plots <- VlnPlot(Integ_data_clus1, features = c( "Foxa2", "Sox17",  "Gata4", "Gata6"), split.by = "type", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)

dev.off()