library(monocle3)
library(dplyr)


expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(urls("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
gene_annotation <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_rowData.rds"))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

#Pre-process the data
cds_p <- preprocess_cds(cds, num_dim = 50)

#batch correction: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
cds_b <- align_cds(cds_p, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")


#dimensionality reduction: by default UMAP
cds_umap <- reduce_dimension(cds_b)

plot_cells(cds_umap, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")


#Let's look at some genes with interesting patterns of expression in ciliated neurons:
ciliated_genes <- c("che-1", "hlh-17", "nhr-6", "dmd-6","ceh-36","ham-1")

plot_cells(cds_umap, genes=ciliated_genes,label_cell_groups=FALSE, show_trajectory_graph=FALSE)


#Cluster your cells: each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
cds_cluster <- cluster_cells(cds_umap)
plot_cells(cds_cluster, color_cells_by = "partition")


#Learn the trajectory graph: we will fit a principal graph within each parition using the learn_graph() function:
cds_traj <- learn_graph(cds_cluster)

plot_cells(cds_traj, color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)

#Order the cells in pseudotime
#Once we've learned a graph, 
#we are ready to order the cells according to their progress through the developmental program. Monocle measures this progress in pseudotime
#Pseudotime is a measure of how much progress an individual cell has made through a process such as cell differentiation.

plot_cells(cds_traj, color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,label_leaves=TRUE,
           label_branch_points=TRUE, graph_label_size=1.5)

#order_cells()needs you to specify the root nodes of the trajectory graph. 
cds_traj_order <- order_cells(cds_traj)

plot_cells(cds_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)


## It's often desirable to specify the root of the trajectory programmatically, rather than manually picking it. 
#The function below does so by first grouping the cells according to which trajectory graph node they are nearest to. 
#Then, it calculates what fraction of the cells at each node come from the earliest time point. 
#Then it picks the node that is most heavily occupied by early cells and returns that as the root.


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds_traj, time_bin="130-170")
  {
  cell_ids <- which(colData(cds_traj)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_traj), ])
  root_pr_nodes <- igraph::V(principal_graph(cds_traj)[["UMAP"]])$name[as.numeric(names (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}



# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds_traj, time_bin="130-170")
  {
  cell_ids <- which(colData(cds_traj)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-cds_traj@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds_traj), ])
  root_pr_nodes <-igraph::V(principal_graph(cds_traj)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

cds_traj_order1 <- order_cells(cds_traj, root_pr_nodes=get_earliest_principal_node(cds_traj))

#Passing the programatically selected root node to order_cells() via the root_pr_nodeargument yields:
  
plot_cells(cds_traj_order1,color_cells_by = "pseudotime",
             label_cell_groups=FALSE,label_leaves=FALSE,
             label_branch_points=FALSE, graph_label_size=1.5)


#Subset cells by branch: It is often useful to subset cells based on their branch in the trajectory.
cds_sub <- choose_graph_segments(cds_traj_order1)


#Working with 3D trajectories
cds_3d <- reduce_dimension(cds_traj_order1, max_components = 3)
cds_3d1 <- cluster_cells(cds_3d)
cds_3d2 <- learn_graph(cds_3d1)
cds_3d3 <- order_cells(cds_3d2, root_pr_nodes=get_earliest_principal_node(cds_traj_order1))
plot_cells_3d(cds_3d3, color_cells_by="partition")
cds_3d_plot_obj <- plot_cells_3d(cds_3d3, color_cells_by="partition")
