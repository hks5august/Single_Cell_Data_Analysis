library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

library(multtest)
library(metap)
library(harmony)
library(SeuratWrappers)




data1 <- Read10X(data.dir = "E1_S1/filtered_feature_bc_matrix/")

seurat_ob_S1 <- CreateSeuratObject(counts = data1, project = "E1_S1", min.cells = 3, min.features = 200)

seurat_ob_S1

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S1[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S1, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob_S1@meta.data, 5)

#jpeg('QC_sample1.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

data2 <- Read10X(data.dir = "E2_S2/filtered_feature_bc_matrix/")

seurat_ob_S2 <- CreateSeuratObject(counts = data2, project = "E2_S2", min.cells = 3, min.features = 200)


head(seurat_ob_S2@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S2[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S2, pattern = "^MT-")

seurat_ob_S2

# Show QC metrics for the first 5 cells
head(seurat_ob_S2@meta.data, 5)

#jpeg('QC_EB_2D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()


data3 <- Read10X(data.dir = "E3_S3/filtered_feature_bc_matrix/")

seurat_ob_S3 <- CreateSeuratObject(counts = data3, project = "E3_S3", min.cells = 3, min.features = 200)


head(seurat_ob_S3@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S3[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S3, pattern = "^MT-")

seurat_ob_S3

# Show QC metrics for the first 5 cells
head(seurat_ob_S3@meta.data, 5)

#jpeg('QC_EB_4D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()



seurat_ob_S1$type <- "EB_Ctrl"
seurat_ob_S2$type <- "EB_2D"
seurat_ob_S3$type <- "EB_4D"


# filter genes with low values (Here, we are removing empty cells/droplet)
EB <- subset(seurat_ob_S1, subset = nFeature_RNA > 250)
EB_2D <- subset(seurat_ob_S2, subset = nFeature_RNA > 250)
EB_4D <- subset(seurat_ob_S3, subset = nFeature_RNA > 250)


#normalize data
EB_norm<- NormalizeData(EB, verbose = FALSE)
EB_2D_norm<- NormalizeData(EB_2D, verbose = FALSE)
EB_4D_norm<- NormalizeData(EB_4D, verbose = FALSE)

# Show QC metrics for the first 5 cells
head(EB_2D_norm@meta.data, 5)

#select highly variable features from each data
EB_ftr<- FindVariableFeatures(EB_norm, selection.method = "vst", nfeatures = 5000)
EB_2D_ftr<- FindVariableFeatures(EB_2D_norm, selection.method = "vst", nfeatures = 5000)
EB_4D_ftr<- FindVariableFeatures(EB_4D_norm, selection.method = "vst", nfeatures = 5000)


# Show QC metrics for the first 5 cells
head(EB_ftr@meta.data, 5)




EB_ctrl_cds <- as.cell_data_set(EB_ftr)
EB_cds_p<- preprocess_cds(EB_ctrl_cds, num_dim = 100)

plot_pc_variance_explained(EB_cds_p)


EB_cds_umap <- reduce_dimension(EB_cds_p)

plot_cells(EB_cds_umap)



EB_cds_clus = cluster_cells(EB_cds_umap, resolution=1e-5)


plot_cells(EB_cds_clus )

plot_cells(EB_cds_clus, color_cells_by="cluster")
plot_cells(EB_cds_clus, color_cells_by="partition")

#plot_cells(EB_cds_umap, color_cells_by="cao_cell_type")





EB_2D_cds <- as.cell_data_set(EB_2D_ftr)
EB_2D_cds_p<- preprocess_cds(EB_2D_cds, num_dim = 100)

plot_pc_variance_explained(EB_2D_cds_p)


EB_2D_cds_umap <- reduce_dimension(EB_2D_cds_p)

plot_cells(EB_2D_cds_umap)

EB_2D_cds_clus = cluster_cells(EB_2D_cds_umap, resolution=1e-5)

plot_cells(EB_2D_cds_clus, color_cells_by="cluster")





EB_4D_cds <- as.cell_data_set(EB_4D_ftr)
EB_4D_cds_p<- preprocess_cds(EB_4D_cds, num_dim = 100)

plot_pc_variance_explained(EB_4D_cds_p)


EB_4D_cds_umap <- reduce_dimension(EB_4D_cds_p)

plot_cells(EB_4D_cds_umap)

EB_4D_cds_clus = cluster_cells(EB_4D_cds_umap, resolution=1e-5)

plot_cells(EB_4D_cds_clus, color_cells_by="cluster")




## integrate samples
integ_anchors <- FindIntegrationAnchors(object.list = list(EB_ftr, EB_2D_ftr, EB_4D_ftr), dims = 1:20)
integ_anchors

Integ_data <- IntegrateData(anchorset = integ_anchors, dims = 1:100)

DefaultAssay(Integ_data) <- "integrated"


#Integ_data_cds <- as(as.matrix(GetAssayData(Integ_data, assay = "integrated", slot = "scale.data")), 'sparseMatrix')

Integ_data_cds <- as(as.matrix(GetAssayData(Integ_data, assay = "integrated")),  'sparseMatrix')
pd <- data.frame(Integ_data@meta.data)
#keep only the columns that are relevant
pData <- pd %>% select(orig.ident, nCount_RNA, nFeature_RNA, type)
fData <- data.frame(gene_short_name = row.names(Integ_data_cds), row.names = row.names(Integ_data_cds))


#Construct monocle cds
monocle.object <- new_cell_data_set(expression_data = Integ_data_cds, cell_metadata = pData, gene_metadata = fData)


#preprocess
monocle_object_p = preprocess_cds(monocle.object, num_dim = 25)
monocle_object_d = reduce_dimension(monocle_object_p)


#Cluster your cells: each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
monocle_cluster <- cluster_cells(monocle_object_d, reduction_method = "UMAP" )

plot_cells(monocle_cluster, color_cells_by = "partition")

plot_cells(monocle_cluster, color_cells_by = "type")


#map pseudotime
monocle.object = order_cells(monocle_cluster)
monocle.object = learn_graph(monocle.object)

#plot trajectories
plot_cells(monocle.object, color_cells_by = "pseudotime")


