library(monocle3)
library(dplyr) # imported for some downstream data manipulation

expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/cao_l2_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


cds_p <- preprocess_cds(cds, num_dim = 50)


plot_pc_variance_explained(cds_p)


cds_b = align_cds(cds_p, num_dim = 50, alignment_group = "plate")


cds_umap <- reduce_dimension(cds_b)

plot_cells(cds_umap)

plot_cells(cds_umap, color_cells_by="cao_cell_type")


cds_cluster <-  cluster_cells(cds_umap, resolution=1e-5)
plot_cells(cds_cluster)

plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster")
