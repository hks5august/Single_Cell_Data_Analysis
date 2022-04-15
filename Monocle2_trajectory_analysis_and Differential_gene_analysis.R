library(monocle3)
library(dplyr)


expression_matrix <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_expression.rds"))
cell_metadata <- readRDS(url("http://staff.washington.edu/hpliner/data/packer_embryo_colData.rds"))
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


#############################################################################
################# Differential Gene expression analysis

#link - https://cole-trapnell-lab.github.io/monocle3/docs/differential/

#Let's begin with a small set of genes that we know are important in ciliated neurons to demonstrate Monocle's capabilities:
ciliated_genes <- c("che-1","hlh-17","nhr-6", "dmd-6", "ceh-36","ham-1")




cds_subset <- cds_traj_order[rowData(cds_traj_order)$gene_short_name %in% ciliated_genes,]
#cds_subset <- cds[rowData(cds)$gene_short_name %in% ciliated_genes,]


#identify genes that vary over time by fitting this model to each one, and then testing whether its 
#?? is significantly different from zero. To do so, we first call the fit_models() function:

gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time")

#Now let's see which of these genes have time-dependent expression. First, we extract a table of coefficients from each model using the coefficient_table() function:
fit_coefs <- coefficient_table(gene_fits)
fit_coefs

#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:
emb_time_terms <- fit_coefs %>% filter(term == "embryo.time")


#Now, let's pull out the genes that have a significant time component.
#coefficient_table() tests whether each coefficient differs significantly from zero under the Wald test
#adjusted values can be found in the q_value column. We can filter the results and control the false discovery rate as follows:


emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

#visualize the differences revealed by the tests above. One type of plot is a "violin" plot.
plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


#Controlling for batch effects and other factors
gene_fits <- fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs <- coefficient_table(gene_fits)
fit_coefs %>% filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate)


#Evaluating models of gene expression. We can evaluate the fits of each model using the evaluate_fits() function:

evaluate_fits(gene_fits)


#Should we include the batch term in our model of gene expression or not? 
#Monocle provides a function compare_models() that can help you decide. 
#Compare models takes two models and returns the result of a likelihood ratio test between them. 
#Any time you add terms to a model, it will improve the fit. 
#But we should always to use the simplest model we can to explain our data. 
#The likelihood ratio test helps us decide whether the improvement in fit is large enough to justify the complexity our extra terms introduce. You run compare_models() like this:


#The first of the two models is called the full model. 
#This model is essentially a way of predicting the expression value of each gene in a given cell knowing both what time it was collected and which batch of cells it came from. 
time_batch_models <- fit_models(cds_subset,
                                model_formula_str = "~embryo.time + batch",
                                expression_family="negbinomial")
time_batch_models

#second model, called the reduced model, does the same thing, 
#but it only knows about the time each cell was collected.
time_models <- fit_models(cds_subset,
                          model_formula_str = "~embryo.time",
                          expression_family="negbinomial")
compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)

#As we can see, all of the genes' likelihood ratio tests are significant, indicating that there are substantial batch effects in the data. 
#We are therefore justified in adding the batch term to our model.

#Choosing a distribution for modeling gene expression


#######################################################
#Finding genes that change as a function of pseudotime
#Identifying the genes that change as cells progress along a trajectory is a core objective of this type of analysis. 
#Knowing the order in which genes go on and off can inform new models of development



plot_cells(cds_traj_order, color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE)
#How do we find the genes that are differentially expressed on the different paths through the trajectory? How do we find the ones that are restricted to the beginning of the trajectory? Or excluded from it?

#Let's perform graph_test(), this time passing it neighbor_graph="principal_graph", which tells it to test whether cells at similar positions on the trajectory have correlated expression:
ciliated_cds_pr_test_res<- graph_test(cds_traj_order, neighbor_graph="principal_graph", cores=8)

pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
pr_deg_ids 
#Here are a couple of interesting genes that score as highly significant according to graph_test():
plot_cells(cds_traj_order, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

#we can collect the trajectory-variable genes into modules:

#gene_module_df <- find_gene_modules(cds_traj_order[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
gene_module_df <- find_gene_modules(cds_traj_order[pr_deg_ids,])

#Here we plot the aggregate module scores within each group of cell types as annotated by Packer & Zhu et al:

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_traj_order)), 
                                cell_group=colData(cds_traj_order)$cell.type)

agg_mat <- aggregate_gene_expression(cds_traj_order, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")


#We can also pass gene_module_df to plot_cells() as we did when we compared clusters in the L2 data above.

plot_cells(cds_traj_order,
           genes=gene_module_df %>% filter(module %in% c(1, 2,5, 6, 7, 8)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


#Monocle offers another plotting function that can sometimes give a clearer view of a gene's dynamics along a single path. 
#You can select a path with choose_cells() or by subsetting the cell data set by cluster, cell type, or other annotation that's restricted to the path. Let's pick one such path, the AFD cells:

AFD_genes <- c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds <- cds_traj_order[rowData(cds_traj_order)$gene_short_name %in% AFD_genes,
                                  colData(cds_traj_order)$cell.type %in% c("AFD")]


#The function plot_genes_in_pseudotime() takes a small set of genes 
#and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)


#Analyzing branches in single-cell trajectories
#Analyzing the genes that are regulated around trajectory branch points often provides insights into the genetic circuits that control cell fate decisions

cds_subset <- choose_cells(cds_traj_order)

#This will identify genes with interesting patterns of expression that fall only within the region of the trajectory you selected, giving you a more refined and relevant set of genes.

subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)

pr_deg_ids1 <- row.names(subset(subset_pr_test_res, q_value < 0.05))

#Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:

gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids1,], resolution=0.001)

#We will organize the modules by their similarity (using hclust) over the trajectory so it's a little easier to see which ones come on before others:

agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module, 
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
