if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")


BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))


install.packages("devtools")


devtools::install_github('cole-trapnell-lab/leidenbase')

devtools::install_github('cole-trapnell-lab/monocle3')


devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")



#install.packages('remotes')
# Replace '2.3.0' with your desired version
#remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
#library(Seurat)


library(monocle3)
library(dplyr)

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

#Clustering and classifying your cells

#https://cole-trapnell-lab.github.io/monocle3/docs/clustering/


#load data from 10X cell ranger
cds <- load_mm_data(mat_path = "E1_S1/filtered_feature_bc_matrix/matrix.mtx.gz", 
                    feature_anno_path = "E1_S1/filtered_feature_bc_matrix/features.tsv.gz", 
                    cell_anno_path = "E1_S1/filtered_feature_bc_matrix/barcodes.tsv.gz")




library(Seurat)
#library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)


#Load Data
data4 <- Read10X(data.dir = "E4_S4/filtered_feature_bc_matrix/")


seurat_ob_S4 <- CreateSeuratObject(counts = data4, project = "E4_S4", min.cells = 3, min.features = 200)

seurat_ob_S4

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S4[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S4, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob_S4@meta.data, 5)

#jpeg('QC_sample1.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()


#load data
data5 <- Read10X(data.dir = "E5_S5/filtered_feature_bc_matrix/")

#Create Seurat object
seurat_ob_S5 <- CreateSeuratObject(counts = data5, project = "E5_S5", min.cells = 3, min.features = 200)


# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S5[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S5, pattern = "^MT-")

seurat_ob_S5

# Show QC metrics for the first 5 cells
head(seurat_ob_S5@meta.data, 5)

#jpeg('QC_Encap_Encap_EB_2D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()


#Load data
data6 <- Read10X(data.dir = "E6_S6/filtered_feature_bc_matrix/")

#Create Seurat object
seurat_ob_S6 <- CreateSeuratObject(counts = data6, project = "E6_S6", min.cells = 3, min.features = 200)


head(seurat_ob_S6@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S6[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S6, pattern = "^MT-")

seurat_ob_S6

# Show QC metrics for the first 5 cells
head(seurat_ob_S6@meta.data, 5)

#jpeg('QC_Encap_EB_4D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()



seurat_ob_S4$type <- "Encap_EB_Ctrl"
seurat_ob_S5$type <- "Encap_EB_2D"
seurat_ob_S6$type <- "Encap_EB_4D"

seurat_ob_S4$treat <- "Control"
seurat_ob_S5$treat <- "RA_treat_2D"
seurat_ob_S6$treat<- "RA_treat_4D"

seurat_ob_S4$treat.time <- "0h"
seurat_ob_S5$treat.time <- "48h"
seurat_ob_S6$treat.time<- "96h"

seurat_ob_S4$treat.time1 <- 0
seurat_ob_S5$treat.time1 <- 48
seurat_ob_S6$treat.time1<- 96

# filter genes with low values (Here, we are removing empty cells/droplet)
Encap_EB <- subset(seurat_ob_S4, subset = nFeature_RNA > 250)
Encap_EB_2D <- subset(seurat_ob_S5, subset = nFeature_RNA > 250)
Encap_EB_4D <- subset(seurat_ob_S6, subset = nFeature_RNA > 250)

#normalize data
Encap_EB_norm<- NormalizeData(Encap_EB, verbose = FALSE)
Encap_EB_2D_norm<- NormalizeData(Encap_EB_2D, verbose = FALSE)
Encap_EB_4D_norm<- NormalizeData(Encap_EB_4D, verbose = FALSE)

# Show QC metrics for the first 5 cells
head(Encap_EB_2D_norm@meta.data, 5)

#select highly variable features from each data
Encap_EB_ftr<- FindVariableFeatures(Encap_EB_norm, selection.method = "vst", nfeatures = 5000)
Encap_EB_2D_ftr<- FindVariableFeatures(Encap_EB_2D_norm, selection.method = "vst", nfeatures = 5000)
Encap_EB_4D_ftr<- FindVariableFeatures(Encap_EB_4D_norm, selection.method = "vst", nfeatures = 5000)


# Show QC metrics for the first 5 cells
head(Encap_EB_ftr@meta.data, 5)



Encap_EB_ctrl_cds <- as.cell_data_set(Encap_EB_ftr)
Encap_EB_cds_p<- preprocess_cds(Encap_EB_ctrl_cds, num_dim = 100)

plot_pc_variance_explained(Encap_EB_cds_p)


Encap_EB_cds_umap <- reduce_dimension(Encap_EB_cds_p)

Encap_EB_cds_umap

plot_cells(Encap_EB_cds_umap)



Encap_EB_cds_clus = cluster_cells(Encap_EB_cds_umap, resolution=1e-5)


plot_cells(Encap_EB_cds_clus )

jpeg('Encap_EB_ctrl_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(Encap_EB_cds_clus, color_cells_by="cluster")
dev.off()

plot_cells(Encap_EB_cds_clus, color_cells_by="partition")



#Find marker genes expressed by each cluster
#, start by calling the top_markers() function:

marker_test_Encap_EB<- top_markers(Encap_EB_cds_clus, group_cells_by="cluster", reference_cells=1000, cores=8)


top_specific_markers<- marker_test_Encap_EB %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(20, pseudo_R2)

top_specific_markers

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

top_specific_marker_ids

marker_df <- as.data.frame(top_specific_marker_ids)

#Write markers into a file
write.table(marker_df,file="marker_list_Encap_EB_ctrl.txt", sep='\t',  quote = F,row.names = FALSE)



#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:

#since gene_short name is required to plot
rowData(Encap_EB_cds_clus)$gene_short_name <- row.names(rowData(Encap_EB_cds_clus))


plot_genes_by_group(Encap_EB_cds_clus, top_specific_marker_ids,
                    group_cells_by="cluster", ordering_type="maximal_on_diag",max.size=2)




Encap_EB_2D_cds <- as.cell_data_set(Encap_EB_2D_ftr)
Encap_EB_2D_cds_p<- preprocess_cds(Encap_EB_2D_cds, num_dim = 100)

plot_pc_variance_explained(Encap_EB_2D_cds_p)


Encap_EB_2D_cds_umap <- reduce_dimension(Encap_EB_2D_cds_p)

plot_cells(Encap_EB_2D_cds_umap)

Encap_EB_2D_cds_clus = cluster_cells(Encap_EB_2D_cds_umap, resolution=1e-5)


jpeg('Encap_EB_2D_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(Encap_EB_2D_cds_clus, color_cells_by="cluster")
plot_cells(Encap_EB_2D_cds_clus, color_cells_by="partition")
dev.off()
#Find marker genes expressed by each cluster
#, start by calling the top_markers() function:

marker_test_2D <- top_markers(Encap_EB_2D_cds_clus, group_cells_by="partition", reference_cells=1000, cores=8)


top_specific_markers_2d <- marker_test_2D %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(10, pseudo_R2)

top_specific_markers_2d

top_specific_marker_ids_2d <- unique(top_specific_markers_2d %>% pull(gene_id))

top_specific_marker_ids_2d

marker_df_2d <- as.data.frame(top_specific_marker_ids_2d)

#Write markers into a file
write.table(marker_df_2d,file="marker_list_2D.txt", sep='\t',  quote = F,row.names = FALSE)



#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:

#since gene_short name is required to plot
rowData(Encap_EB_2D_cds_clus)$gene_short_name <- row.names(rowData(Encap_EB_2D_cds_clus))


plot_genes_by_group(Encap_EB_2D_cds_clus, top_specific_marker_ids_2d,
                    group_cells_by="cluster", ordering_type="maximal_on_diag",max.size=2)








Encap_EB_4D_cds <- as.cell_data_set(Encap_EB_4D_ftr)
Encap_EB_4D_cds_p<- preprocess_cds(Encap_EB_4D_cds, num_dim = 100)

plot_pc_variance_explained(Encap_EB_4D_cds_p)


Encap_EB_4D_cds_umap <- reduce_dimension(Encap_EB_4D_cds_p)

plot_cells(Encap_EB_4D_cds_umap)

Encap_EB_4D_cds_clus = cluster_cells(Encap_EB_4D_cds_umap, resolution=1e-5)


jpeg('Encap_EB_4D_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(Encap_EB_4D_cds_clus, color_cells_by="cluster")
dev.off()



#Find marker genes expressed by each cluster
#, start by calling the top_markers() function:

marker_test_4D <- top_markers(Encap_EB_4D_cds_clus, group_cells_by="cluster", reference_cells=1000, cores=8)


top_specific_markers_4d <- marker_test_4D %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(10, pseudo_R2)

top_specific_markers_4d

top_specific_marker_ids_4d <- unique(top_specific_markers_4d %>% pull(gene_id))

top_specific_marker_ids_4d

marker_df_4d <- as.data.frame(top_specific_marker_ids_4d)

#Write markers into a file
write.table(marker_df_4d,file="marker_list_4D.txt", sep='\t',  quote = F,row.names = FALSE)



#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:

#since gene_short name is required to plot
rowData(Encap_EB_4D_cds_clus)$gene_short_name <- row.names(rowData(Encap_EB_4D_cds_clus))


plot_genes_by_group(Encap_EB_4D_cds_clus, top_specific_marker_ids_4d,
                    group_cells_by="cluster", ordering_type="maximal_on_diag",max.size=2)



################################################################################
## integrate samples

integ_anchors_encap <- FindIntegrationAnchors(object.list = list(Encap_EB_ftr, Encap_EB_2D_ftr, Encap_EB_4D_ftr), dims = 1:25)
Integ_data_encap <- IntegrateData(anchorset = integ_anchors_encap, dims = 1:100)
DefaultAssay(Integ_data_encap) <- "integrated"


Integ_data_cds_encap <- as(as.matrix(GetAssayData(Integ_data_encap, assay = "integrated")),  'sparseMatrix')
pd_encap <- data.frame(Integ_data_encap@meta.data)
head(pd_encap)

#keep only the columns that are relevant
pData_encap <- pd_encap %>% select(orig.ident, nCount_RNA, nFeature_RNA, type, treat, treat.time1)
fData_encap <- data.frame(gene_short_name = row.names(Integ_data_cds_encap), row.names = row.names(Integ_data_cds_encap))
#Construct monocle cds
monocle.object_encap <- new_cell_data_set(expression_data = Integ_data_cds_encap, cell_metadata = pData_encap, gene_metadata = fData_encap)

#preprocess
monocle_object_p_encap = preprocess_cds(monocle.object_encap, num_dim = 50)

#PCA plot
plot_pc_variance_explained(monocle_object_p_encap)
#monocle_umap_encap = reduce_dimension(monocle_object_p_encap)

#batch correction: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
#cds_b <- align_cds(cds_p, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
monocle_object_p1_encap <- align_cds(monocle_object_p_encap, alignment_group = "type")

#UMAP
monocle_umap_encap1 = reduce_dimension(monocle_object_p1_encap)
monocle_umap_encap = reduce_dimension(monocle_object_p_encap)


# clustering
monocle_cluster_encap1 <- cluster_cells(monocle_umap_encap1)
monocle_cluster_encap <- cluster_cells(monocle_umap_encap)

#UMAP plot based on sample
plot_cells(monocle_cluster_encap , label_groups_by_cluster=FALSE,  color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


plot_cells(monocle_cluster_encap1 , label_groups_by_cluster=FALSE,  color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


#UMAP plot based on cluster
plot_cells(monocle_cluster_encap1, label_groups_by_cluster=FALSE,  color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(monocle_cluster_encap, label_groups_by_cluster=FALSE,  color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


##########################################
#Find marker genes for each cluster
marker_clusters_Encap <- top_markers(monocle_cluster_encap, group_cells_by="cluster",  reference_cells=5000, cores=8)
write.table(marker_clusters_Encap,file="markers_Encap_9_cluster.txt", sep='\t',  quote = F,row.names = FALSE)

#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers_Encap <- marker_clusters_Encap %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(25, pseudo_R2)
top_specific_markers_Encap

#plot markers for all clusters together in a dotplot
top_specific_marker_ids <- unique(top_specific_markers_encap %>% pull(gene_id))

#Dotplot
plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=4)) 

######## Markers For each cluster one by one ####
#Cluster1 markers
top_specific_markers_Encap <- marker_clusters_Encap %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(10, pseudo_R2)

dim(top_specific_markers_Encap)
top_specific_markers1 <- subset(top_specific_markers_Encap, cell_group==1)

top_specific_marker_ids1 <- unique(top_specific_markers1 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids1,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster2 markers
top_specific_markers2 <- subset(top_specific_markers_Encap, cell_group==2)

top_specific_marker_ids2 <- unique(top_specific_markers2 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster3 markers
top_specific_markers3 <- subset(top_specific_markers_Encap, cell_group==3)

top_specific_marker_ids3 <- unique(top_specific_markers3 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids3,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster4 markers
top_specific_markers4 <- subset(top_specific_markers_Encap, cell_group==4)

top_specific_marker_ids4 <- unique(top_specific_markers4 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids4,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 

#Cluster4 markers
top_specific_markers5 <- subset(top_specific_markers_Encap, cell_group==5)

top_specific_marker_ids5 <- unique(top_specific_markers5 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids5,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster6 markers
top_specific_markers6 <- subset(top_specific_markers_Encap, cell_group==6)

top_specific_marker_ids6 <- unique(top_specific_markers6 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids6,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster7 markers
top_specific_markers7 <- subset(top_specific_markers_Encap, cell_group==4)

top_specific_marker_ids7 <- unique(top_specific_markers7 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids7,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster8 markers
top_specific_markers8 <- subset(top_specific_markers_Encap, cell_group==8)

top_specific_marker_ids8 <- unique(top_specific_markers8 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids8,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 

#Cluster9 markers
top_specific_markers9 <- subset(top_specific_markers_Encap, cell_group==9)

top_specific_marker_ids9 <- unique(top_specific_markers9 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster_encap,
                    top_specific_marker_ids9,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


##################################################

#############################################
###Assign cell types based on clustermole
library(clustermole)

Type1 <- c()
clusters1 <- c()
for (i in seq(1,9, by=1)) {
  markers_Encap.filt <- as.data.frame(top_specific_markers_Encap %>% filter(cell_group == as.character(i)))
  my_overlaps <- clustermole_overlaps(genes = markers_Encap.filt$gene_id, species = "mm")
  type <- as.character(my_overlaps[1,5])
  Type1 <- c(Type1, type)
  clusters1 <- c(clusters1, i)
}


clust_map1 <- as.data.frame(cbind(clusters1,Type1))
head(clust_map1,10)



write.table(clust_map1,file="clust_map_cell_types_based_on_clustermole1.txt", sep='\t',  quote = F,row.names = FALSE)

colData(monocle_cluster_encap)$assigned_cell_type <- as.character(clusters(monocle_cluster_encap))
colData(monocle_cluster_encap)$assigned_cell_type = dplyr::recode(colData(monocle_cluster_encap)$assigned_cell_type,
                                                               "1"="embryonic_stem_line",
                                                               "2"="Hu Fetal Retina fibroblast",
                                                               "3"="OVARY (BULK TISSUE)",
                                                               "4"="fetal kidney 2 stromal cells",
                                                               "5"="Kidney C8 thin ascending limb",
                                                               "6"="Heart vascular Endothelial",
                                                               "7"="Kidney C13 thick ascending limb",
                                                               "8"="Schwalie et al.Nature.G4",
                                                               "9"="Aizarani liver C10 MVECS1")


#type based
plot_cells(monocle_cluster_encap, color_cells_by = "assigned_cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


#Assign cell types based on singleR
library(SingleR)
library(celldex)
library(scuttle)
ref <- MouseRNAseqData()
ref


colData(monocle_cluster_encap)$cluster <- as.character(clusters(monocle_cluster_encap))

pred_encap <- SingleR(test=monocle_cluster_encap, ref=ref, assay.type.test=1, labels=ref$label.fine)

pred_encap

table(pred_encap$labels)
table(pred_encap$first.labels)
table(pred_encap$pruned.labels)
plotScoreHeatmap(pred_encap)

plotDeltaDistribution(pred_encap, ncol = 4)

tab1_encap_cell_types <- table(Assigned=pred_encap$labels, Cluster=clusters(monocle_cluster_encap))
tab1_encap_cell_types

tab1_encap_cell_types_df <- as.data.frame(tab1_encap_cell_types)
write.table(tab1_encap_cell_types_df,file="tab1_encap_cell_types_single_R.txt", sep='\t',  quote = F,row.names = FALSE)


# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab1_encap_cell_types+10), color=colorRampPalette(c("white", "blue"))(50))

colData(monocle_cluster_encap)$cell_type <- as.character(pred_encap$labels)

plot_cells(monocle_cluster_encap, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


################################################

#Differential gene expression analysis between partitions
all_genes_d_encap <- unique(rowData(monocle_cluster_encap)$gene_short_name)
all_genes_d_encap 

all_gene_subset_d_encap <- monocle_cluster_encap[rowData(monocle_cluster_encap)$gene_short_name %in% all_genes_d_encap,]


all_gene_fits_cls_d_encap <- fit_models(all_gene_subset_d_encap, model_formula_str = "~cluster")
head(all_gene_fits_cls_d_encap)

#all_gene_fits <- fit_models(all_gene_subset, model_formula_str = "~treat.time1")

all_gene_fit_coefs_d_encap <- coefficient_table(all_gene_fits_cls_d_encap)


head(all_gene_fit_coefs_d_encap)
#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:

all_gene_time_terms_d_encap <- all_gene_fit_coefs_d_encap %>% filter(term != "(Intercept)")
all_gene_time_terms1_d_encap <- all_gene_time_terms_d_encap %>% filter(status == "OK")

#Extract significant genes based on q-values
sig_genes_d_encap <- all_gene_time_terms1_d_encap %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
sig_genes_df_d_encap  <- as.data.frame(sig_genes_d_encap)

write.table(sig_genes_df_d_encap,file="all_sig_genes_between_clusters_encap.txt", sep='\t',  quote = F,row.names = FALSE)
###############################
sig_genes_d_encap


#Visualize top genes
sig_1<- c("Rhox5", "Sparc", "Col4a1", "Ctsl", "Dab2", "Lama1", "Htra1", "Gata6", "Foxq1", "Sox17")

sig_2 <- c( "H3f3b", "Krt18",   "Fau", "Pga5", "Glipr2", "Gpc3", "Lgmn", "Actb", "Cubn", "Peg3")




#colData(monocle_cluster)$cluster <- as.character(partitions(monocle_cluster))

top_sig_set_1 <- monocle_cluster_encap[rowData(monocle_cluster_encap)$gene_short_name %in% sig_1,]
top_sig_set_2 <- monocle_cluster_encap[rowData(monocle_cluster_encap)$gene_short_name %in% sig_2,]


#top_sig_set <- monocle_cluster[rowData(monocle_cluster)$gene_short_name %in% sig,]

plot_genes_violin(top_sig_set_1, group_cells_by="cluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

plot_genes_violin(top_sig_set_2, group_cells_by="cluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


#label_cell_groups

markers_ESC <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4") ##Esc
markers_stem <- c("CD146", "CD45", "PDGFR-beta") ##stem markers
markers_ectodem <- c("Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6") ##Ectoderm
markers_endoderm <- c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1") ##Endoderm
markers_mesoderm <- c("Mest" , "Bmp2","Eomes", "Hand1" , "Isl1", "Kdr", "Mesdc1" , "Mesdc2", "Myf5", "Myod1", "Nkx2-5", "T", "Tbx2", "Prdm1", "Tbx1") ##mesoderm
markers_neural_stem <- c("Apoe", "Blbp", "Gfap", "Slc1a3", "Sox2", "Sox9", "Thrsp") ## neural stem cells
markers_microglial <- c("1700112E06Rik", "4632428N05Rik", "Abca9", "Abi3", "Adap2", "AF251705", "Aif1", "Apbb1ip", "Arsb", "Bmp2k", "C1qa", "C1qb", "C1qc", "Ccr5", "Cd14", "Cd37", "Cd53", "Cd68", "Csf1r", "Ctsb", "Ctss", "Cx3cr1", "Cyth4", "Dock2", "Dock8", "Emr1", "Entpd1", "Fcer1g", "Fcgr3", "Fcrls", "Fli1", "Fyb", "Gpr34", "Hexb", "Hmha1", "Hpgd", "Hpgds", "Ikzf1", "Il6ra", "Inpp5d", "Irf8", "Itgam", "Itgb5", "Lair1", "Laptm5", "Lcp2", "Lgmn", "Lpcat2", "Ltc4s", "Ly86", "Lyn", "Mafb", "Mertk", "Mpeg1", "Myo1f", "Ncf1", "Nckap1l", "Olfml3", "P2ry12", "P2ry13", "Pik3ap1", "Pld4", "Pros1", "Ptgs1", "Ptprc", "Rasgrp3", "Rnase4", "RrEncap_EB1", "Runx1", "Selplg", "Serinc3", "Siglech", "Sirpa", "Skap2", "Slco2b1", "Tbxas1", "Tgfbr1", "Tgfbr2", "Tmem119", "Trem2", "Tyrobp", "Unc93b1", "Vav1", "Zfp710", "0610040J01Rik")
markers_Early_intermediate_precursor <- c("Sox11", "Tbr2") #Early intermediate precursor cell
markers_Early_neuroblast  <- c("Acvr2a", "Fgf10", "Fgf3", "Fgf8", "Hes5", "Neurog1", "Notch1", "Six1") #Early_neuroblast
markers_Interneuron_selective_cell <- c("Bcar3", "Celsr3", "Chrna4", "Cpox", "Dlgap3", "Dok7", "Egfr", "Kcnh3", "Lrrc3b", "Miat", "Npy2r", "Nr2e1", "Ostf1", "Ptms", "Rasd2", "Rgs16", "Rprm", "Sema5b", "Slc12a4", "St3gal1", "Tpbg", "Unc5a", "Wscd2", "2010300C02Rik")#Interneuron-selective cell
markers_Long_projecting_GABAergic <- c("1700086L19Rik" , "A930011G23Rik" , "Adamts16" , "Arhgdig" , "Arhgef15" , "Arl4a" , "Ass1" , "B2m" , "Bace2" , "Carhsp1" , "Ccdc109b" , "Chodl" , "Clic1" , "Cnih3" , "Col11a1" , "Cort" , "Ctxn2" , "D430019H16Rik" , "Dbpht2" , "Dpy19l1" , "Egfl7" , "Fam46a" , "Fndc1" , "Gabra2" , "Gabrg1" , "Gm11744" , "Gpr126" , "Gpr88" , "Hey1" , "Htr1a" , "Kcnmb4" , "N28178" , "Ndrg1" , "Ndst4" , "Nnat" , "Nos1" , "Nptx1" , "Ntn1" , "Opn3" , "Oxtr" , "Patl2" , "Pnma3" , "Pou3f2" , "Prex1" , "Ptn" , "Ptpru" , "Rasgef1b" , "Rbp1" , "S100a10" , "Sc5d" , "Sfrp1" , "Slc5a5" , "Slc5a7" , "Slc7a3" , "Sst" , "St6galnac2" , "Sv2b" , "Tppp3" , "Tspan17" , "Wnt2" , "Wwc1" , "1700058G18Rik") #Long-projecting GABAergic cell
markers_endothelil <- c("Abcb1a", "Abcg2", "Acvrl1", "Apcdd1", "Car4", "Ccdc141", "CD34", "Cd93", "Cdh5", "Cgnl1", "Cldn5", "Clic5", "Ctla2a", "Cxcl12", "Cyyr1", "Egfl7", "Eltd1", "Emcn", "Eng", "Epas1", "Erg", "Esam", "Ets1", "Fli1", "Flt1", "Fn1", "Foxq1", "Fzd6", "Gpr116", "Hmcn1", "Hspb1", "Hspg2", "Ifitm3", "Ifltd1", "Igfbp7", "Itga1", "Itga4", "Itm2a", "Kank3", "Kdr", "Kitl", "Klf2", "Lama4", "Lsr", "Ly6a", "Ly6c1", "Ly6c2", "Ly75", "Mecom", "Mfsd2a", "Nos3", "Nostrin", "Palmd", "Paqr5", "Pcp4l1", "Pecam1", "Pglyrp1", "Pltp", "Podxl", "Ptprb", "Ramp2", "Rasgrp3", "Rassf9", "Rbpms", "Rgs5", "Rhoj", "Sdpr", "Sema3c", "Sgms1", "Slc16a1", "Slc22a8", "Slc2a1", "Slc40a1", "Slc7a1", "Slc7a5", "Slc9a3r2", "Slco1a4", "Slco1c1", "Slco2b1", "Sox17", "Sparc", "Srgn", "St3gal6", "Tek", "Tm4sf1", "Vwf", "Wfdc1", "Wwtr1", "Zfp366", "9430020K01Rik") ###endothelial cells
markers_neurons <- c("6330403K07Rik", "A030009H04Rik", "Agap2", "Ahi1", "Arf3", "Atp1a3", "Atp1b1", "Bex2", "Camk2a", "Camk2b", "Camk2n1", "Celf4", "Chd3", "Cnih2", "Cplx2", "Ctxn1", "Dynll2", "Fbxl16", "Gap43", "Gng13", "Hap1", "Hpcal4", "Kif5a", "Kif5c", "Ly6h", "Map1a", "Map1b", "Map2", "Myt1l", "Ncdn", "Nefl", "Nrgn", "Nrsn1", "Pcp4", "Pcsk1n", "Rab3c", "Rasgef1a", "Rph3a", "Rtn1", "S100b", "Scg2", "Sepw1", "Sez6", "Snap25", "Snrpn", "Stmn2", "Stmn3", "Syt1", "Ttc3", "Unc13c") ##neurons
markers_neuroblast <- c("Ccnd2", "Crmp1", "Dbn1", "Dchs1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Evf1", "Foxc1", "Fxyd6", "Klf5a", "Olfm1", "Pfn2", "Pou3f3", "Sox5", "Trn1", "Btg1", "Celf4", "Dcx", "Meis2", "Stmn2") #neuroblast
markers_Late_activated_neural_stem <- c("2810417H13Rik", "Cdk1", "Cenpf", "Hist1h2ap", "Hmgb2", "Top2a")##Late activated neural stem cell
markers_Late_neuroblast <- c("Eya2", "Fzd10", "Isl1", "Jag2", "Neurod1", "Neurog2") #Late neuroblast
markers_Activated_neural_stem <- c("Csf3", "Fbxo2", "Fxyd1", "Gfap", "Id3", "Thbs4") #Activated neural stem cell
markers_precursor <- c("Neurod1", "Neurog1") #precusor cells
markers_mature_olf_sensory_neuron <- c("Omp") #Mature olfactory sensory neuron
markers_Progenitor <- c("Ascl1") #Progenitor cell
markers_immature_olf_sensory_neuron <- c("Gap43", "Gng8") #immature olfactory sensory neuron
markers_ependymal <- c( "3300002A11Rik" , "4930529M08Rik", "9330101J02Rik", "Acta2", "Adamts20", "Ak7", "Ak9", "Armc3", "Ccdc108", "Ccdc114", "Ccdc146", "Ccdc153", "Ccdc162", "Ccdc170", "Ccdc60", "Clu", "Daw1", "Dnah10", "Dnah11", "Dnah12", "Dnah3", "Dnah5", "Dnah6", "E230008N13Rik", "Efhb", "Enkur", "Fhad1", "Foxj1", "Gja1", "Gm973", "Hydin", "Ifltd1", "Iqck", "Kif6", "Kif9", "Lrguk", "Lrrc48", "Lrriq1", "Mia", "Nek11", "Rarres2", "Rgs22", "S100b", "Sox9", "Spag16", "Spag17", "Spata17", "Spef2", "Tm4sf1", "Tmem212", "Tmem232", "Ttc21a", "Ttc29", "Ttll8", "Vwa3a", "Wdr52", "Wdr63", "Wdr96", "Zbbx", "1700007G11Rik") ##Ependymal cells
markers_type1_ganglion <-c("Vars", "Vdac2", "Vegfa", "Vkorc1l1", "Vps36", "Vps52", "Vsnl1", "Vstm5", "Vti1a", "Vwa7", "Wbscr17", "Wdr26", "Wdr75", "Wee1", "Wfs1", "Wipf3", "Wwp1", "Xylt2", "Yars2", "Yipf1", "Ypel4", "Ythdc2", "Ythdf1", "Ythdf3", "Zbtb17", "Zcchc2", "Zcrb1", "Zdhhc1", "Zdhhc5", "Zfp157", "Zfp2", "Zfp24", "Zfp30", "Zfp322a", "Zfp354c", "Zfp385a", "Zfp397", "Zfp422", "Zfp426", "Zfp511", "Zfp521", "Zfp62", "Zfp672", "Zfp687", "Zfp770", "Zfp804a", "Zfp930", "Zfp935", "Zfp958", "Zfyve28", "Zrsr2", "Zscan21", "Zyx", "Kcnc3", "Nefl", "Scn4b", "Tubb3")
markers_typeII_ganglion <-c("Zfp142", "Zfp185", "Zfp319", "Zfp57", "Zfp646", "Zfyve1", "Zmynd11", "Gata3", "Mafb", "Prph", "Th")
markers_oligodendrocytes <- c("Aspa", "Ccp110", "Cdn11", "Cldn11", "Cntn2", "Efnb3", "Enpp2", "Ermn", "Evi2a", "Fa2h", "Galnt6", "Gjb1", "Gjc2", "Gjc3", "Gm21984", "Gpr37", "Grb14", "Gsn", "Hapln2", "Josd2", "Lgi3", "Lpar1", "Mag", "Mal", "Mbp", "Mog", "Ndrg1", "Olig1", "Opalin", "Pdlim2", "Pex5l", "Plekhh1", "Pllp", "Plp1", "Pls1", "Ppp1r14a", "Ptgds", "S1pr5", "Sept4", "Slc12a2", "Smco3", "Sox10", "Stmn4", "Tmeff2", "Tmem125", "Tmem151a", "Tmem88b", "Tnfaip6", "Trf", "Trim59", "Tspan2", "Ttll7", "Tubb4a", "Ugt8a", "Apod")



#Let's look at some genes with interesting patterns of expression in spefic cell type

# change marker genes
marker_genes <- markers_ependymal  #### change marker genes according to cells we are looking for
marker_genes <- markers_oligodendrocytes 

plot_cells(monocle_cluster_encap,  
           min_expr= 1.0, 
           norm_method = "log",
           #norm_method = "size_only",
           label_groups_by_cluster=TRUE, 
           genes=marker_genes, rasterize = TRUE,
           label_cell_groups=FALSE, 
           show_trajectory_graph=FALSE)


plot_cells(monocle_cluster, color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")



plot_cells(monocle_cluster, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

jpeg('Integrated_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_cluster, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

jpeg('Integrated_type.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_cluster, color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()



############# Count cells ######################


# cluster distribution in each group
#table(Idents(monocle_cluster), monocle_cluster$type)
colData(monocle_cluster_encap)$cluster <- as.character(clusters( monocle_cluster_encap))

group_dist_encap  <- table(clusters(monocle_cluster_encap), monocle_cluster_encap$type)

group_dist_encap

#Write results into a file
write.table(group_dist_encap,file="cell_distribution_in_9_clusters.txt", sep='\t',  quote = F,row.names = FALSE)


#Find marker genes for each cluster

marker_clusters_encap <- top_markers(monocle_cluster_encap, group_cells_by="cluster",  reference_cells=1000, cores=8)
write.table(marker_clusters_encap,file="markers_integrated_data_cluster.txt", sep='\t',  quote = F,row.names = FALSE)



#The data frame marker_clusters contains a number of metrics for how specifically expressed each gene is in each partition


#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers_encap <- marker_clusters_encap %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)


top_specific_markers_clus1 <- subset(top_specific_markers_encap, cell_group == 1)
top_specific_markers_clus2 <- subset(top_specific_markers_encap, cell_group == 2)
top_specific_markers_clus3 <- subset(top_specific_markers_encap, cell_group == 3)
top_specific_markers_clus4 <- subset(top_specific_markers_encap, cell_group == 4)
top_specific_markers_clus5 <- subset(top_specific_markers_encap, cell_group == 5)
top_specific_markers_clus6 <- subset(top_specific_markers_encap, cell_group == 6)
top_specific_markers_clus7 <- subset(top_specific_markers_encap, cell_group == 7)
top_specific_markers_clus8 <- subset(top_specific_markers_encap, cell_group == 8)



top_specific_marker_ids1 <- unique(top_specific_markers_clus1 %>% pull(gene_id))

top_specific_marker_ids2 <- unique(top_specific_markers_clus2 %>% pull(gene_id))

top_specific_marker_ids3 <- unique(top_specific_markers_clus3 %>% pull(gene_id))

top_specific_marker_ids4 <- unique(top_specific_markers_clus4 %>% pull(gene_id))
top_specific_marker_ids5 <- unique(top_specific_markers_clus5 %>% pull(gene_id))
top_specific_marker_ids6 <- unique(top_specific_markers_clus6 %>% pull(gene_id))
top_specific_marker_ids7 <- unique(top_specific_markers_clus7 %>% pull(gene_id))
top_specific_marker_ids8 <- unique(top_specific_markers_clus8 %>% pull(gene_id))


plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids1,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids2,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids3,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 



plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids4,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids5,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids6,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids7,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 



plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids8,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




top_specific_marker_ids 
top_specific_markers 

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=4)) 


#Learn the trajectory graph: we will fit a principal graph within each parition using the learn_graph() function:
monocle_traj_encap <- learn_graph(monocle_cluster_encap,  use_partition = FALSE)

jpeg('Integrated_trajectory_type.jpg', units="in", width=10, height=8, res=300)

plot_cells(monocle_traj_encap, color_cells_by = "type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

dev.off()


jpeg('Integrated_trajectory_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_traj_encap, color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

dev.off()



jpeg('Integrated_trajectory_time.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_traj, color_cells_by = "treat.time1",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()


#oredr cells

monocle_traj_order_encap <- order_cells(monocle_traj_encap)

plot_cells(monocle_traj_order_encap , color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)

cluster_table <- table(round(pseudotime(monocle_traj_order_encap),2), monocle_traj_order_encap$type)
write.table(cluster_table,file="cluster_1_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_table

cluster_table_c <- table(round(pseudotime(monocle_traj_order_encap),2), monocle_traj_order_encap$cluster)
write.table(cluster_table_c,file="cluster_c_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_table



#order_cells()needs you to specify the root nodes of the trajectory graph. 
monocle_traj_order1 <- order_cells(monocle_traj)

plot_cells(monocle_traj_order1, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)

cluster_1_table <- table(round(pseudotime(monocle_traj_order1),2), monocle_traj_order1$type)
write.table(cluster_1_table,file="cluster_1_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_1_table



monocle_traj_order2 <- order_cells(monocle_traj)

plot_cells(monocle_traj_order2, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


cluster_2_table <- table(round(pseudotime(monocle_traj_order2),2), monocle_traj_order2$type)
write.table(cluster_2_table,file="cluster_2_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_2_table



monocle_traj_order3 <- order_cells(monocle_traj)

plot_cells(monocle_traj_order3, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


cluster_3_table <- table(round(pseudotime(monocle_traj_order3),2), monocle_traj_order3$type)
write.table(cluster_3_table,file="cluster_3_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_3_table



monocle_traj_order4 <- order_cells(monocle_traj)

plot_cells(monocle_traj_order4, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


cluster_4_table <- table(round(pseudotime(monocle_traj_order4),2), monocle_traj_order4$type)
write.table(cluster_4_table,file="cluster_4_table_pseudotime.txt", sep='\t',  quote = F,row.names = T)
cluster_4_table

#######################################################
monocle_traj_order_all <- order_cells(monocle_traj)

plot_cells(monocle_traj_order_all, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


############ Choose nodes with function #########


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(monocle_traj, treat.time1="0"){
  cell_ids1 <- which(colData(monocle_traj)[, "treat.time1"] == treat.time1)
  
  closest_vertex <- monocle_traj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(monocle_traj), ])
  root_pr_nodes <-igraph::V(principal_graph(monocle_traj)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids1,]))))]
  root_pr_nodes
}
cds11 <- order_cells(monocle_traj, root_pr_nodes=get_earliest_principal_node(monocle_traj))

plot_cells(cds11,
           color_cells_by = "pseudotime",
           label_cell_groups=TRUE,
           label_leaves=F,
           label_branch_points=TRUE,
           graph_label_size=1.5)



################# View trajectory for specific Gene type ########

################# View trajectory for specific Gene type ########
#First let's see ESC markers
ESC_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_ESC,]
#ESC_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_ESC,]
#ESC_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_ESC,]
#ESC_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_ESC,]



#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
plot_genes_in_pseudotime(ESC_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(ESC_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(ESC_cds1, color_cells_by="cluster", min_expr=2, cell_size = 1)

#plot_genes_in_pseudotime(ESC_cds2, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(ESC_cds3, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(ESC_cds4, color_cells_by="type", min_expr=0.5)



#Neural stem cells markers

neural_stem_lineage_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neural_stem,]
#neural_stem_lineage_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_neural_stem,]
#neural_stem_lineage_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_neural_stem,]
#neural_stem_lineage_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_neural_stem,]

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

###type ###
plot_genes_in_pseudotime(neural_stem_lineage_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neural_stem_lineage_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

#plot_genes_in_pseudotime(neural_stem_lineage_cds2, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(neural_stem_lineage_cds3, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(neural_stem_lineage_cds4, color_cells_by="type", min_expr=0.5)

#### Neuroblast
markers_neuroblast1 <- c("Ccnd2", "Crmp1", "Dbn1", "Dchs1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Evf1", "Foxc1") #neuroblast
markers_neuroblast2 <- c( "Fxyd6", "Klf5a", "Olfm1", "Pfn2", "Pou3f3", "Sox5", "Trn1", "Btg1", "Celf4", "Dcx", "Meis2", "Stmn2") #neuroblast


neuroblast1_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neuroblast1,]
neuroblast2_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neuroblast2,]

###type ###
plot_genes_in_pseudotime(neuroblast1_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neuroblast2_cds1, color_cells_by="type", min_expr=0.5)

plot_genes_in_pseudotime(neuroblast1_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(neuroblast2_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)



### Endoderm marker
endoderm_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endoderm,]
#endoderm_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_endoderm,]
#endoderm_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_endoderm,]
#endoderm_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_endoderm,]


###type ###
plot_genes_in_pseudotime(endoderm_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(endoderm_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

#plot_genes_in_pseudotime(endoderm_cds2, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(endoderm_cds3, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(endoderm_cds4, color_cells_by="type", min_expr=0.5)

### Ectoderm marker
ectoderm_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_ectodem,]
#ectoderm_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_ectodem,]
#ectoderm_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_ectodem,]
#ectoderm_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_ectodem,]


###type ###
plot_genes_in_pseudotime(ectoderm_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(ectoderm_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

#plot_genes_in_pseudotime(ectoderm_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(ectoderm_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(ectoderm_cds4, color_cells_by="type", min_expr=0.5)


### mesoderm marker
mesoderm_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_mesoderm,]
#mesoderm_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_mesoderm,]
#mesoderm_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_mesoderm,]
#mesoderm_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_mesoderm,]


###type ###
plot_genes_in_pseudotime(mesoderm_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(mesoderm_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
#plot_genes_in_pseudotime(mesoderm_cds2, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(mesoderm_cds3, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(mesoderm_cds4, color_cells_by="type", min_expr=0.5)



#Neurons markers

markers_neurons1 <- c("6330403K07Rik", "A030009H04Rik", "Agap2", "Ahi1", "Arf3", "Atp1a3", "Atp1b1", "Bex2", "Camk2a", "Camk2b", "Camk2n1", "Celf4", "Chd3", "Cnih2", "Cplx2", "Ctxn1") ##neurons

neurons1_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neurons1]
#neurons_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_neurons1]
#neurons_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_neurons1]
#neurons_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_neurons1]

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
## type
plot_genes_in_pseudotime(neurons1_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons1_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

#plot_genes_in_pseudotime(neurons_cds2, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(neurons_cds3, color_cells_by="type", min_expr=0.5)
#plot_genes_in_pseudotime(neurons_cds4, color_cells_by="type", min_expr=0.5)


markers_neurons2 <- c( "Dynll2", "Fbxl16", "Gap43", "Gng13", "Hap1", "Hpcal4", "Kif5a", "Kif5c", "Ly6h", "Map1a", "Map1b", "Map2", "Myt1l")

neurons2_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neurons2]
neurons2_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_neurons2]
neurons2_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_neurons2]
neurons2_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_neurons2]

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type
plot_genes_in_pseudotime(neurons2_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons2_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

plot_genes_in_pseudotime(neurons2_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons2_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons2_cds4, color_cells_by="type", min_expr=0.5)


#neurons marker set3
markers_neurons3 <- c("Ncdn", "Nefl", "Nrgn", "Nrsn1", "Pcp4", "Pcsk1n", "Rab3c", "Rasgef1a", "Rph3a", "Rtn1", "S100b", "Scg2", "Sepw1", "Sez6", "Snap25", "Snrpn", "Stmn2", "Stmn3", "Syt1", "Ttc3", "Unc13c")

neurons3_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_neurons3]
#neurons3_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% markers_neurons3]
#neurons3_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% markers_neurons3]
#neurons3_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% markers_neurons3]

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type
plot_genes_in_pseudotime(neurons3_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons3_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

plot_genes_in_pseudotime(neurons3_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons3_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(neurons3_cds4, color_cells_by="type", min_expr=0.5)


#First let's see Long projecting Gabergic neurons markers

gab1 <- c("1700086L19Rik" , "A930011G23Rik" , "Adamts16" , "Arhgdig" , "Arhgef15" , "Arl4a" , "Ass1" , "B2m" , "Bace2" , "Carhsp1" , "Ccdc109b" , "Chodl" , "Clic1" , "Cnih3" , "Col11a1" , "Cort" , "Ctxn2" , "D430019H16Rik" , "Dbpht2" , "Dpy19l1" , "Egfl7" , "Fam46a" , "Fndc1" , "Gabra2" , "Gabrg1" , "Gm11744" , "Gpr126" , "Gpr88" , "Hey1" , "Htr1a" , "Kcnmb4" , "N28178" , "Ndrg1" , "Ndst4" , "Nnat" , "Nos1" , "Nptx1" , "Ntn1" , "Opn3" , "Oxtr" , "Patl2" , "Pnma3" , "Pou3f2" , "Prex1" , "Ptn" , "Ptpru" , "Rasgef1b" , "Rbp1") #Long-projecting GABAergic cell
gab2 <- c( "S100a10" , "Sc5d" , "Sfrp1" , "Slc5a5" , "Slc5a7" , "Slc7a3" , "Sst" , "St6galnac2" , "Sv2b" , "Tppp3" , "Tspan17" , "Wnt2" , "Wwc1" , "1700058G18Rik") #Long-projecting GABAergic cell

Gab1_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% gab1,]
Gab1_cds2 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% gab1,]
Gab1_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% gab1,]
Gab1_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% gab1,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

#type
plot_genes_in_pseudotime(Gab1_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab1_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)

plot_genes_in_pseudotime(Gab1_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab1_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab1_cds4, color_cells_by="type", min_expr=0.5)



##### se2
gab2 <- c( "S100a10" , "Sc5d" , "Sfrp1" , "Slc5a5" , "Slc5a7" , "Slc7a3" , "Sst" , "St6galnac2" , "Sv2b" , "Tppp3" , "Tspan17" , "Wnt2" , "Wwc1" , "1700058G18Rik") #Long-projecting GABAergic cell

Gab2_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% gab2,]
Gab2_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% gab2,]
Gab2_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% gab2,]
Gab2_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% gab2,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type
plot_genes_in_pseudotime(Gab2_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab2_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)


plot_genes_in_pseudotime(Gab2_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab2_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gab2_cds4, color_cells_by="type", min_expr=0.5)


####

gang1  <-c("Vars", "Vdac2", "Vegfa", "Vkorc1l1", "Vps36", "Vps52", "Vsnl1", "Vstm5", "Vti1a", "Vwa7", "Wbscr17", "Wdr26", "Wdr75")


Gang1_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% gang1,]
Gang1_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% gang1,]
Gang1_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% gang1,]
Gang1_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% gang1,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type
plot_genes_in_pseudotime(Gang1_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang1_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)


plot_genes_in_pseudotime(Gang1_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang1_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang1_cds4, color_cells_by="type", min_expr=0.5)

### ganglion type 1 set2
gang2 <- c( "Wee1", "Wfs1", "Wipf3", "Wwp1", "Xylt2", "Yars2", "Yipf1", "Ypel4", "Ythdc2", "Ythdf1", "Ythdf3", "Zbtb17", "Zcchc2", "Zcrb1", "Zdhhc1", "Zdhhc5", "Zfp157", "Zfp2", "Zfp24", "Zfp30", "Zfp322a", "Zfp354c", "Zfp385a", "Zfp397", "Zfp422", "Zfp426", "Zfp511", "Zfp521", "Zfp62", "Zfp672", "Zfp687", "Zfp770", "Zfp804a", "Zfp930", "Zfp935", "Zfp958", "Zfyve28", "Zrsr2", "Zscan21", "Zyx", "Kcnc3", "Nefl", "Scn4b", "Tubb3")


Gang2_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% gang2,]
Gang2_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% gang2,]
Gang2_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% gang2,]
Gang2_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% gang2,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type- gang_type1_set2_clus
plot_genes_in_pseudotime(Gang2_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang2_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)



plot_genes_in_pseudotime(Gang2_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang2_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(Gang2_cds4, color_cells_by="type", min_expr=0.5)


######### Type 2 ganglion #####
typeII_ganglion <-c("Zfp142", "Zfp185", "Zfp319", "Zfp57", "Zfp646", "Zfyve1", "Zmynd11", "Gata3", "Mafb", "Prph", "Th")

type2_Gang_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% typeII_ganglion,]
type2_Gang_cds2 <- monocle_traj_order2[rowData(monocle_traj_order2)$gene_short_name %in% typeII_ganglion,]
type2_Gang_cds3 <- monocle_traj_order3[rowData(monocle_traj_order3)$gene_short_name %in% typeII_ganglion,]
type2_Gang_cds4 <- monocle_traj_order4[rowData(monocle_traj_order4)$gene_short_name %in% typeII_ganglion,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type gang_type2_clus
plot_genes_in_pseudotime(type2_Gang_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(type2_Gang_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)



plot_genes_in_pseudotime(type2_Gang_cds2, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(type2_Gang_cds3, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(type2_Gang_cds4, color_cells_by="type", min_expr=0.5)


##Oligodendrocytes
markers_oligodendrocytes <- c("Aspa", "Ccp110", "Cdn11", "Cldn11", "Cntn2", "Efnb3", "Enpp2", "Ermn", "Evi2a", "Fa2h", "Galnt6", "Gjb1", "Gjc2", "Gjc3", "Gm21984", "Gpr37", "Grb14", "Gsn", "Hapln2", "Josd2", "Lgi3", "Lpar1", "Mag", "Mal", "Mbp", "Mog", "Ndrg1", "Olig1", "Opalin", "Pdlim2", "Pex5l", "Plekhh1", "Pllp", "Plp1", "Pls1", "Ppp1r14a", "Ptgds", "S1pr5", "Sept4", "Slc12a2", "Smco3", "Sox10", "Stmn4", "Tmeff2", "Tmem125", "Tmem151a", "Tmem88b", "Tnfaip6", "Trf", "Trim59", "Tspan2", "Ttll7", "Tubb4a", "Ugt8a", "Apod")


oligo_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_oligodendrocytes,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type gang_type2_clus
plot_genes_in_pseudotime(oligo_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(oligo_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)


### Endothelial
markers_endothelil <- c("Abcb1a", "Abcg2", "Acvrl1", "Apcdd1", "Car4", "Ccdc141", "CD34", "Cd93", "Cdh5", "Cgnl1", "Cldn5", "Clic5", "Ctla2a", "Cxcl12", "Cyyr1", "Egfl7", "Eltd1", "Emcn", "Eng", "Epas1", "Erg", "Esam", "Ets1", "Fli1", "Flt1", "Fn1", "Foxq1", "Fzd6", "Gpr116", "Hmcn1", "Hspb1", "Hspg2", "Ifitm3", "Ifltd1", "Igfbp7", "Itga1", "Itga4", "Itm2a", "Kank3", "Kdr", "Kitl", "Klf2", "Lama4", "Lsr", "Ly6a", "Ly6c1", "Ly6c2", "Ly75", "Mecom", "Mfsd2a", "Nos3", "Nostrin", "Palmd", "Paqr5", "Pcp4l1", "Pecam1", "Pglyrp1", "Pltp", "Podxl", "Ptprb", "Ramp2", "Rasgrp3", "Rassf9", "Rbpms", "Rgs5", "Rhoj", "Sdpr", "Sema3c", "Sgms1", "Slc16a1", "Slc22a8", "Slc2a1", "Slc40a1", "Slc7a1", "Slc7a5", "Slc9a3r2", "Slco1a4", "Slco1c1", "Slco2b1", "Sox17", "Sparc", "Srgn", "St3gal6", "Tek", "Tm4sf1", "Vwf", "Wfdc1", "Wwtr1", "Zfp366", "9430020K01Rik") ###endothelial cells

markers_endothelil1 <- c("Abcb1a", "Abcg2", "Acvrl1", "Apcdd1", "Car4", "Ccdc141", "CD34", "Cd93", "Cdh5", "Cgnl1", "Cldn5", "Clic5", "Ctla2a", "Cxcl12", "Cyyr1", "Egfl7", "Eltd1", "Emcn", "Eng", "Epas1", "Erg", "Esam") ###endothelial cells
markers_endothelil2 <- c( "Ets1", "Fli1", "Flt1", "Fn1", "Foxq1", "Fzd6", "Gpr116", "Hmcn1", "Hspb1", "Hspg2") ###endothelial cells
markers_endothelil3 <- c( "Ifitm3", "Ifltd1", "Igfbp7", "Itga1", "Itga4", "Itm2a", "Kank3", "Kdr", "Kitl", "Klf2", "Lama4", "Lsr") ###endothelial cells
markers_endothelil4 <- c("Ly6a", "Ly6c1", "Ly6c2", "Ly75", "Mecom", "Mfsd2a", "Nos3", "Nostrin", "Palmd", "Paqr5", "Pcp4l1", "Pecam1", "Pglyrp1", "Pltp", "Podxl", "Ptprb", "Ramp2", "Rasgrp3", "Rassf9", "Rbpms", "Rgs5", "Rhoj", "Sdpr", "Sema3c", "Sgms1")
markers_endothelil5 <- c( "Slc16a1", "Slc22a8", "Slc2a1", "Slc40a1", "Slc7a1", "Slc7a5", "Slc9a3r2", "Slco1a4", "Slco1c1", "Slco2b1", "Sox17", "Sparc", "Srgn", "St3gal6", "Tek", "Tm4sf1", "Vwf", "Wfdc1", "Wwtr1", "Zfp366", "9430020K01Rik")


endothelial1_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endothelil1,]
endothelial2_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endothelil2,]
endothelial3_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endothelil3,]
endothelial4_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endothelil4,]
endothelial5_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_endothelil5,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type gang_type2_clus
plot_genes_in_pseudotime(endothelial1_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(endothelial1_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(endothelial2_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(endothelial3_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(endothelial4_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)
plot_genes_in_pseudotime(endothelial5_cds1, color_cells_by="cluster", min_expr=0.5, cell_size = 1)


plot_genes_in_pseudotime(endothelial5_cds1, color_cells_by="cluster", min_expr=3, cell_size = 1)



### Ependymal
ependymal_cds1 <- monocle_traj_order_encap[rowData(monocle_traj_order_encap)$gene_short_name %in% markers_ependymal,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type gang_type2_clus
plot_genes_in_pseudotime(ependymal_cds1 , color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(ependymal_cds1 , color_cells_by="cluster", min_expr=0.5, cell_size = 1)




### precursors
precursor_cds1 <- monocle_traj_order1[rowData(monocle_traj_order1)$gene_short_name %in% markers_precursor,]


#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:
#type gang_type2_clus
plot_genes_in_pseudotime(precursor_cds1, color_cells_by="type", min_expr=0.5)
plot_genes_in_pseudotime(precursor_cds1, color_cells_by="cluster", min_expr=2, cell_size = 1)




#Find marker genes for each cluster

#Find marker genes for each partition
marker_clusters <- top_markers(monocle_cluster, group_cells_by="partition",  reference_cells=1000, cores=8)
write.table(marker_clusters,file="markers_integrated_data_partition.txt", sep='\t',  quote = F,row.names = FALSE)

#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers <- marker_clusters %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)


#plot markers for all clusters together in a dotplot
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

#Dotplot
plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=4)) 

######## Markers For each cluster one by one ####
#Cluster1 markers
top_specific_markers1 <- subset(top_specific_markers, cell_group==1)

top_specific_marker_ids1 <- unique(top_specific_markers1 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids1,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster2 markers
top_specific_markers2 <- subset(top_specific_markers, cell_group==2)

top_specific_marker_ids2 <- unique(top_specific_markers2 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids2,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster3 markers
top_specific_markers3 <- subset(top_specific_markers, cell_group==3)

top_specific_marker_ids3 <- unique(top_specific_markers3 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids3,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 


#Cluster4 markers
top_specific_markers4 <- subset(top_specific_markers, cell_group==4)

top_specific_marker_ids4 <- unique(top_specific_markers4 %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids4,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=4) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=8), axis.text.y=element_text(angle=0, hjust=1, size=7)) 




#The data frame marker_clusters contains a number of metrics for how specifically expressed each gene is in each partition


#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers1 <- marker_clusters1 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids1 <- unique(top_specific_markers1 %>% pull(gene_id))

plot_genes_by_group(monocle_traj,
                    top_specific_marker_ids1,
                    group_cells_by="treat.time1",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=3)) 



#order_cells()needs you to specify the root nodes of the trajectory graph. 
monocle_traj_order <- order_cells(monocle_traj)

plot_cells(monocle_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


marker_clusters_cls2 <- top_markers(monocle_traj_order, group_cells_by="pseudotime",  reference_cells=1000, cores=8)
write.table(marker_clusters_cls2,file="markers_integrated_traj_clus2_pseudotime.txt", sep='\t',  quote = F,row.names = FALSE)



#The data frame marker_clusters contains a number of metrics for how specifically expressed each gene is in each partition


#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers1 <- marker_clusters1 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids1 <- unique(top_specific_markers1 %>% pull(gene_id))

plot_genes_by_group(monocle_traj,
                    top_specific_marker_ids1,
                    group_cells_by="treat.time1",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=3)) 



#Differential gene expression analysis based on clusters

all_genes <- unique(rowData(monocle_traj)$gene_short_name)
all_genes 

all_gene_subset <- monocle_traj[rowData(monocle_traj)$gene_short_name %in% all_genes,]


all_gene_fits_cls <- fit_models(all_gene_subset, model_formula_str = "~cluster")
head(all_gene_fits_cls)

#all_gene_fits <- fit_models(all_gene_subset, model_formula_str = "~treat.time1")

all_gene_fit_coefs <- coefficient_table(all_gene_fits_cls)
#all_gene_fit_coefs <- coefficient_table(all_gene_fits_cls)

all_gene_fit_coefs
#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:


all_gene_time_terms <- all_gene_fit_coefs %>% filter(term != "(Intercept)")
all_gene_time_terms1 <- all_gene_time_terms %>% filter(status == "OK")

all_gene_time_terms <- all_gene_fit_coefs %>% filter(term == "cluster*")
all_gene_time_terms1

#We can see that five of the 7 genes significantly vary as a function of time.

sig_genes <- all_gene_time_terms1 %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
sig_genes_df  <- as.data.frame(sig_genes)
sig_genes_df
write.table(sig_genes_df,file="all_sig_genes_df.txt", sep='\t',  quote = F,row.names = FALSE)


########################## Visualization ##################


sig <- c("Peg3", "Fau", "Plod2", "Hint1", "Sox17", "Foxq1", "Dab2", "Pdcd4", "Hist1h1d", "Gata6")

colData( monocle_traj)$cluster <- as.character(clusters( monocle_traj))
top_sig_set1 <- monocle_traj[rowData(monocle_traj)$gene_short_name %in% sig,]


colData(monocle_cluster)$cluster <- as.character(clusters(monocle_cluster))
top_sig_set <- monocle_cluster[rowData(monocle_cluster)$gene_short_name %in% sig,]

#Plot genes

plot_genes_violin(top_sig_set, group_cells_by="cluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

######## Differential Gene expression analysis based on partition

all_genes <- unique(rowData(monocle_traj)$gene_short_name)
all_genes 



all_gene_fits_cls1 <- fit_models(all_gene_subset1, model_formula_str = "~partition")
head(all_gene_fits_cls)

#all_gene_fits <- fit_models(all_gene_subset, model_formula_str = "~treat.time1")

all_gene_fit_coefs1 <- coefficient_table(all_gene_fits_cls1)

all_gene_fit_coefs1
#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:

all_gene_time_terms2 <- all_gene_fit_coefs1 %>% filter(term != "(Intercept)")
all_gene_time_terms3 <- all_gene_time_terms2 %>% filter(status == "OK")
all_gene_time_terms3 
#Extract significant genes based on q-values
sig_genes1 <- all_gene_time_terms3 %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
sig_genes_df1  <- as.data.frame(sig_genes1)
sig_genes_df1
write.table(sig_genes_df1,file="all_sig_genes_df_between_partition.txt", sep='\t',  quote = F,row.names = FALSE)

#Visualize top genes
sig1 <- c("Flrt3", "Retreg1", "Nostrin", "Lifr", "2200002D01Rik", "Cubn", "Pga5", "Klf6", "Pdlim4", "Gjb5", "Sox17")

colData( monocle_traj)$partitions<- as.character(partitions( monocle_traj))

#colData(monocle_cluster)$cluster <- as.character(partitions(monocle_cluster))


top_sig_set1 <- monocle_traj[rowData(monocle_traj)$gene_short_name %in% sig1,]
#top_sig_set <- monocle_cluster[rowData(monocle_cluster)$gene_short_name %in% sig,]

plot_genes_violin(top_sig_set1, group_cells_by="cluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))













####################################################
#order_cells()needs you to specify the root nodes of the trajectory graph. 
monocle_traj_order <- order_cells(monocle_traj)

plot_cells(monocle_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)














###################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)

pt.matrix <- exprs(monocle_traj)
pt.matrix1 <- pt.matrix[match(sig,rownames(rowData(monocle_traj))) ,order(pseudotime(monocle_traj))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix1

pt.matrix2 <- t(apply(pt.matrix1,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 6,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

heatmap <- plot_pseudotime_heatmap(top_sig_set1, cluster_rows = FALSE, show_rownames = TRUE, return_heatmap = T)
heatmap

 

plot_multiple_branches_heatmap(top_sig_set)
top_sig_set_m <- as.matrix(top_sig_set)
pheatmap::pheatmap(top_sig_set, cluster_cols = F, cluster_rows = F,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(10))

###########################################################################
############ automated assignment of cell types ###################

## Install the monocle3 branch of garnett
BiocManager::install(c("org.Mm.eg.db", "org.Hs.eg.db"))
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")





########################################################################

