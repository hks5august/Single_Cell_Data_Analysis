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

#Load Data

data1 <- Read10X(data.dir = "E1_S1/filtered_feature_bc_matrix/")

#Create Seurat Object
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

#Load data
data4 <- Read10X(data.dir = "E4_S4/filtered_feature_bc_matrix/")

#Create Seurat Object
seurat_ob_S4 <- CreateSeuratObject(counts = data4, project = "E4_S4", min.cells = 3, min.features = 200)

head(seurat_ob_S2@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S4[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S4, pattern = "^MT-")

seurat_ob_S4

# Show QC metrics for the first 5 cells
head(seurat_ob_S4@meta.data, 5)

#jpeg('QC_EB_2D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

#Add Desired columns to the Seurat object
#Sample Type
seurat_ob_S1$type <- "EB_Ctrl"
seurat_ob_S4$type <- "Encap_ctrl"

# filter genes with low values (Here, we are removing empty cells/droplet)
EB_ctrl <- subset(seurat_ob_S1, subset = nFeature_RNA > 250)
Encap_EB_ctrl <- subset(seurat_ob_S4, subset = nFeature_RNA > 250)

#normalize data
EB_ctrl_norm<- NormalizeData(EB_ctrl, verbose = FALSE)
Encap_EB_ctrl_norm<- NormalizeData(Encap_EB_ctrl, verbose = FALSE)



#select highly variable features from each data
EB_ctrl_ftr<- FindVariableFeatures(EB_ctrl_norm, selection.method = "vst", nfeatures = 5000)
Encap_EB_ctrl_ftr<- FindVariableFeatures(Encap_EB_ctrl_norm, selection.method = "vst", nfeatures = 5000)


#let's first analyze each individual sample
# Convert seurat object to monocle dataset
EB_ctrl_cds <- as.cell_data_set(EB_ctrl_ftr)
#Non-linear dimensionality reduction, by default method is uMAP (one can also use t-SNE)


#Preprocess data, linear dimensionality reduction - by deafulat PCA
EB_cds_p<- preprocess_cds(EB_ctrl_cds, num_dim = 100)

#Visualize PCA components
plot_pc_variance_explained(EB_cds_p)

#Non-linear dimensionality reduction, by default method is uMAP (one can also use t-SNE)
EB_cds_umap <- reduce_dimension(EB_cds_p)

#Scatterplot - UMAP
plot_cells(EB_cds_umap)

EB_cds_clus
#Clustering
EB_cds_clus = cluster_cells(EB_cds_umap, resolution=1e-5)

#Plot cells based on Clustering
plot_cells(EB_cds_clus, color_cells_by="cluster" ) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
#Plot cells based on type
plot_cells(EB_cds_clus, color_cells_by="type" ) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


#save plot
jpeg('EB_ctrl_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(EB_cds_clus, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

#Plot based on Sample type
plot_cells(EB_cds_clus, color_cells_by="type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(EB_cds_clus, color_cells_by="partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")



EB_mat <- as(as.matrix(GetAssayData(seurat_ob_S1, assay = "RNA")),  'sparseMatrix')

#create metadata dataframe 
pd <- data.frame(seurat_ob_S1@meta.data)
pd
#keep important columns in metadata that are relevant 
pData <- pd %>% select(nCount_RNA, nFeature_RNA, type)

#create features dataframe
fData <- data.frame(gene_short_name = row.names(EB_mat), row.names = row.names(EB_mat))

#Construct monocle dataset
monocle_EB_ctrl <- new_cell_data_set(expression_data = EB_mat, cell_metadata = pData, gene_metadata = fData)


#preprocess - Linear dimensionality reduction
monocle_EB_p = preprocess_cds(monocle_EB_ctrl, num_dim = 50)

#PCA plot
plot_pc_variance_explained(monocle_EB_p)

#batch correction: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
#cds_b <- align_cds(cds_p, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
#monocle_object_p1 <- align_cds(monocle_object_p, alignment_group = "type")
#Non-linear dimensionality reduction - UMAP
#monocle_umap1 <- reduce_dimension(monocle_object_p1)
#plot_cells(monocle_umap1, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
#  theme(legend.text=element_text(size=6)) + #set the size of the text
# theme(legend.position="right")

#Non-linear dimensionality reduction - UMAP
monocle_EB_umap <- reduce_dimension(monocle_EB_p)

#scatterplot uMAP
plot_cells(monocle_EB_umap, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Cluster your cells: each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
monocle_EB_cluster <- cluster_cells(monocle_EB_umap)

monocle_EB_cluster1 <- cluster_cells(monocle_EB_umap,  
                                     cluster_method =  "louvain", k=120)

plot_cells(monocle_EB_cluster1, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

monocle_Encap_cluster1 <- cluster_cells(monocle_Encap_umap,  
                                     cluster_method =  "louvain", k=120)
plot_cells(monocle_Encap_cluster1, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#partition based
plot_cells(monocle_EB_cluster, color_cells_by = "partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Clusters based
plot_cells(monocle_EB_cluster , color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(monocle_EB_cluster , color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


#Known markers
markers_ESC <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4") ##Esc
markers_stem <- c("CD146", "CD45", "PDGFR-beta") ##stem markers
markers_ectodem <- c("Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6") ##Ectoderm
markers_endoderm <- c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1") ##Endoderm
markers_mesoderm <- c("Mest" , "Bmp2","Eomes", "Hand1" , "Isl1", "Kdr", "Mesdc1" , "Mesdc2", "Myf5", "Myod1", "Nkx2-5", "T", "Tbx2", "Prdm1", "Tbx1") ##mesoderm
markers_neural_stem <- c("Apoe", "Blbp", "Gfap", "Slc1a3", "Sox2", "Sox9", "Thrsp") ## neural stem cells
markers_microglial <- c("1700112E06Rik", "4632428N05Rik", "Abca9", "Abi3", "Adap2", "AF251705", "Aif1", "Apbb1ip", "Arsb", "Bmp2k", "C1qa", "C1qb", "C1qc", "Ccr5", "Cd14", "Cd37", "Cd53", "Cd68", "Csf1r", "Ctsb", "Ctss", "Cx3cr1", "Cyth4", "Dock2", "Dock8", "Emr1", "Entpd1", "Fcer1g", "Fcgr3", "Fcrls", "Fli1", "Fyb", "Gpr34", "Hexb", "Hmha1", "Hpgd", "Hpgds", "Ikzf1", "Il6ra", "Inpp5d", "Irf8", "Itgam", "Itgb5", "Lair1", "Laptm5", "Lcp2", "Lgmn", "Lpcat2", "Ltc4s", "Ly86", "Lyn", "Mafb", "Mertk", "Mpeg1", "Myo1f", "Ncf1", "Nckap1l", "Olfml3", "P2ry12", "P2ry13", "Pik3ap1", "Pld4", "Pros1", "Ptgs1", "Ptprc", "Rasgrp3", "Rnase4", "Rreb1", "Runx1", "Selplg", "Serinc3", "Siglech", "Sirpa", "Skap2", "Slco2b1", "Tbxas1", "Tgfbr1", "Tgfbr2", "Tmem119", "Trem2", "Tyrobp", "Unc93b1", "Vav1", "Zfp710", "0610040J01Rik")
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
marker_genes <- markers_ESC #### change marker genes according to cells we are looking for


exprs(monocle_EB_cluster)

marker_genes <- markers_typeII_ganglion

plot_cells(monocle_EB_cluster,  
           min_expr= 1.5, 
           norm_method = "log",
           label_groups_by_cluster=FALSE, 
           genes=marker_genes,
           label_cell_groups=FALSE, show_trajectory_graph=FALSE)

#actual expression
plot_cells(monocle_Encap_cluster,  
           min_expr= 1.5, 
           norm_method = "log",
           label_groups_by_cluster=FALSE, 
           genes=marker_genes,
           label_cell_groups=FALSE,  show_trajectory_graph=FALSE)


plot_cells(monocle_EB_cluster, genes=marker_genes,label_cell_groups=FALSE, show_trajectory_graph=FALSE)





###############################


#let's first analyze each individual sample
# Convert seurat object to monocle dataset
Encap_ctrl_cds <- as.cell_data_set(Encap_EB_ctrl_ftr)
#Non-linear dimensionality reduction, by default method is uMAP (one can also use t-SNE)


#Preprocess data, linear dimensionality reduction - by deafulat PCA
Encap_cds_p<- preprocess_cds(Encap_ctrl_cds, num_dim = 100)

#Visualize PCA components
plot_pc_variance_explained(EB_cds_p)

#Non-linear dimensionality reduction, by default method is uMAP (one can also use t-SNE)
Encap_cds_umap <- reduce_dimension(Encap_cds_p)

#Scatterplot - UMAP
plot_cells(Encap_cds_umap)


#Clustering
Encap_cds_clus = cluster_cells(Encap_cds_umap, resolution=1e-5)

#Plot cells based on Clustering
plot_cells(Encap_cds_clus, color_cells_by="cluster" ) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
#Plot cells based on type
plot_cells(Encap_cds_clus, color_cells_by="type" ) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Plot cells based on type
plot_cells(Encap_cds_clus, color_cells_by="partition" ) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")





#save plot
jpeg('EB_ctrl_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(EB_cds_clus, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

#Plot based on Sample type
plot_cells(EB_cds_clus, color_cells_by="type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(EB_cds_clus, color_cells_by="partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")




Encap_mat <- as(as.matrix(GetAssayData(seurat_ob_S2, assay = "RNA")),  'sparseMatrix')

#create metadata dataframe 
Encap_pd <- data.frame(seurat_ob_S2@meta.data)
Encap_pd
#keep important columns in metadata that are relevant 
Encap_pData <- Encap_pd %>% select(nCount_RNA, nFeature_RNA, type)

#create features dataframe
Encap_fData <- data.frame(gene_short_name = row.names(Encap_mat), row.names = row.names(Encap_mat))

#Construct monocle dataset
monocle_Encap_ctrl <- new_cell_data_set(expression_data = Encap_mat, cell_metadata = Encap_pData, gene_metadata = Encap_fData)


#preprocess - Linear dimensionality reduction
monocle_Encap_p = preprocess_cds(monocle_Encap_ctrl, num_dim = 50)

#PCA plot
plot_pc_variance_explained(monocle_Encap_p)

#batch correction: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
#cds_b <- align_cds(cds_p, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
#monocle_object_p1 <- align_cds(monocle_object_p, alignment_group = "type")
#Non-linear dimensionality reduction - UMAP
#monocle_umap1 <- reduce_dimension(monocle_object_p1)
#plot_cells(monocle_umap1, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
#  theme(legend.text=element_text(size=6)) + #set the size of the text
# theme(legend.position="right")

#Non-linear dimensionality reduction - UMAP
monocle_Encap_umap <- reduce_dimension(monocle_Encap_p)

#scatterplot uMAP
plot_cells(monocle_Encap_umap, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Cluster your cells: each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
# monocle_Encap_cluster <- cluster_cells(monocle_Encap_umap)

#partition based
plot_cells(monocle_Encap_cluster, color_cells_by = "partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#partition based
plot_cells(monocle_Encap_cluster, color_cells_by = "partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Clusters based
plot_cells(monocle_Encap_cluster , color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")



# change marker genes
marker_genes <- markers_typeII_ganglion #### change marker genes according to cells we are looking for

plot_cells(monocle_Encap_cluster, genes=marker_genes,label_cell_groups=FALSE, show_trajectory_graph=FALSE)




#### EB trajectory ##############
plot_cells(monocle_EB_cluster, genes=marker_genes,label_cell_groups=FALSE, show_trajectory_graph=FALSE)


monocle_EB_traj <- learn_graph(monocle_EB_cluster)

plot_cells(monocle_EB_traj, color_cells_by = "type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

monocle_EB_traj_order <- order_cells(monocle_EB_traj)

plot_cells(monocle_EB_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)



#Looks at pattern of specific cell type marker as pseudotime function
#First let's see ESC markers
ESC_EB_cds <- monocle_EB_traj_order[rowData(monocle_EB_traj_order)$gene_short_name %in% markers_ESC,]
ESC_EB_cds

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(ESC_EB_cds, color_cells_by="pseudotime", min_expr=0.5)



Neural_stem_EB_cds <- monocle_EB_traj_order[rowData(monocle_EB_traj_order)$gene_short_name %in% markers_neural_stem,]
Neural_stem_EB_cds

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(Neural_stem_EB_cds, color_cells_by="pseudotime", min_expr=0.5)



############# Encap Trajectory ##################


#### EB trajectory ##############

monocle_Encap_traj <- learn_graph(Encap_cds_clus)

plot_cells(monocle_Encap_traj, color_cells_by = "type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE)+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

monocle_Encap_traj_order <- order_cells(monocle_Encap_traj)

plot_cells(monocle_Encap_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)



#Looks at pattern of specific cell type marker as pseudotime function
#First let's see ESC markers
ESC_Encap_cds <- monocle_Encap_traj_order[rowData(monocle_Encap_traj_order)$gene_short_name %in% markers_ESC,]
ESC_Encap_cds

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(ESC_Encap_cds, color_cells_by="pseudotime", min_expr=0.5)



Neural_stem_Encap_cds <- monocle_Encap_traj_order[rowData(monocle_Encap_traj_order)$gene_short_name %in% markers_neural_stem,]
Neural_stem_Encap_cds

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(Neural_stem_Encap_cds, color_cells_by="pseudotime", min_expr=0.5)





monocle_Encap_traj <- learn_graph(monocle_Encap_cluster)






#identify marker specific for cluster
marker_test_EB<- top_markers(EB_cds_umap, group_cells_by="cluster", reference_cells=1000, cores=8)

#filter top 20 markers
top_specific_markers<- marker_test_EB %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(20, pseudo_R2)

top_specific_markers

#Get ids of markers
top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

#create dataframe
marker_df <- as.data.frame(top_specific_marker_ids)

#Write markers into a file
write.table(marker_df,file="marker_list_EB_ctrl.txt", sep='\t',  quote = F,row.names = FALSE)

#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:
#since gene_short name is required to plot
rowData(EB_cds_clus)$gene_short_name <- row.names(rowData(EB_cds_clus))
#plot
plot_genes_by_group(EB_cds_clus, top_specific_marker_ids,
                    group_cells_by="partition", ordering_type="maximal_on_diag",max.size=2)

# Analysis of 2nd sample: EB-2D
Encap_EB_ctrl_cds <- as.cell_data_set(Encap_EB_ctrl_ftr)
Encap_EB_ctrl_p<- preprocess_cds(Encap_EB_ctrl_cds, num_dim = 100)
plot_pc_variance_explained(Encap_EB_ctrl_p)
Encap_EB_ctrl_cds_umap <- reduce_dimension(Encap_EB_ctrl_p)
plot_cells(Encap_EB_ctrl_cds_umap)
#clustering
Encap_EB_ctrl_cds_clus = cluster_cells(Encap_EB_ctrl_cds_umap, resolution=1e-5)
#plot
plot_cells(Encap_EB_ctrl_cds_clus, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Save plot
jpeg('EB_2D_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(EB_2D_cds_clus, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

#Find marker genes expressed for each cluster
marker_test_2D <- top_markers(EB_2D_cds_clus, group_cells_by="partition", reference_cells=1000, cores=8)

# top markers
top_specific_markers_2d <- marker_test_2D %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(10, pseudo_R2)

top_specific_markers_2d
#marker ids
top_specific_marker_ids_2d <- unique(top_specific_markers_2d %>% pull(gene_id))
# convert into a data frame
marker_df_2d <- as.data.frame(top_specific_marker_ids_2d)

#Write markers into a file
write.table(marker_df_2d,file="marker_list_2D.txt", sep='\t',  quote = F,row.names = FALSE)

#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:
rowData(EB_2D_cds_clus)$gene_short_name <- row.names(rowData(EB_2D_cds_clus))
#plot
plot_genes_by_group(EB_2D_cds_clus, top_specific_marker_ids_2d,
                    group_cells_by="cluster", ordering_type="maximal_on_diag",max.size=2)


#analysis of 3rd sample: EB-4D
EB_4D_cds <- as.cell_data_set(EB_4D_ftr)
EB_4D_cds_p<- preprocess_cds(EB_4D_cds, num_dim = 100)
plot_pc_variance_explained(EB_4D_cds_p)
EB_4D_cds_umap <- reduce_dimension(EB_4D_cds_p)
plot_cells(EB_4D_cds_umap)

#clustering
EB_4D_cds_clus = cluster_cells(EB_4D_cds_umap, resolution=1e-5)
#UMAP Scatter plot
plot_cells(EB_4D_cds_clus, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
#save plot
jpeg('EB_4D_clusters.jpg', units="in", width=10, height=8, res=300)
plot_cells(EB_4D_cds_clus, color_cells_by="cluster")
dev.off()

#Find marker genes expressed by each cluster
marker_test_4D <- top_markers(EB_4D_cds_clus, group_cells_by="cluster", reference_cells=1000, cores=8)
#top markers
top_specific_markers_4d <- marker_test_4D %>% filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>% top_n(10, pseudo_R2)
#marker ids
top_specific_marker_ids_4d <- unique(top_specific_markers_4d %>% pull(gene_id))
#convert into a dataframe
marker_df_4d <- as.data.frame(top_specific_marker_ids_4d)

#Write markers into a file
write.table(marker_df_4d,file="marker_list_4D.txt", sep='\t',  quote = F,row.names = FALSE)

#plot the expression and fraction of cells that express each marker in each group with the plot_genes_by_group function:
rowData(EB_4D_cds_clus)$gene_short_name <- row.names(rowData(EB_4D_cds_clus))
#Dot pSlot
plot_genes_by_group(EB_4D_cds_clus, top_specific_marker_ids_4d,
                    group_cells_by="cluster", ordering_type="maximal_on_diag",max.size=2)


################## integrate samples ####################
integ_anchors_ctrls <- FindIntegrationAnchors(object.list = list(EB_ctrl_ftr, Encap_EB_ctrl_ftr), dims = 1:20)
integ_anchors_ctrls

#Integrate data based on anchors
Integ_data_ctrl <- IntegrateData(anchorset = integ_anchors_ctrls, dims = 1:100)

#define integrated assay type as default assay
DefaultAssay(Integ_data_ctrl) <- "integrated"

#create dataset from seurat object:
#create matrix
Integ_data_cds_ctrl <- as(as.matrix(GetAssayData(Integ_data_ctrl, assay = "integrated")),  'sparseMatrix')

#create metadata dataframe 
pd_ctrl <- data.frame(Integ_data_ctrl@meta.data)

#keep important columns in metadata that are relevant 
pData_ctrl <- pd_ctrl %>% select(orig.ident, nCount_RNA, nFeature_RNA, type, treat, treat.time1)

#create features dataframe
fData_ctrl <- data.frame(gene_short_name = row.names(Integ_data_cds_ctrl), row.names = row.names(Integ_data_cds_ctrl))

#Construct monocle dataset
monocle.object_ctrl <- new_cell_data_set(expression_data = Integ_data_cds_ctrl, cell_metadata = pData_ctrl, gene_metadata = fData_ctrl)



#preprocess - Linear dimensionality reduction
monocle_object_p_ctrl = preprocess_cds(monocle.object_ctrl, num_dim = 50)

#PCA plot
plot_pc_variance_explained(monocle_object_p_ctrl)

#batch correction: Your data will not have the loading batch information demonstrated here, you will correct batch using your own batch information.
#cds_b <- align_cds(cds_p, alignment_group = "batch", residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
#monocle_object_p1_ctrl <- align_cds(monocle_object_p_ctrl, alignment_group = "type")
#Non-linear dimensionality reduction - UMAP
#monocle_umap1_ctrl <- reduce_dimension(monocle_object_p1_ctrl)
#plot_cells(monocle_umap1_ctrl, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
#  theme(legend.text=element_text(size=6)) + #set the size of the text
# theme(legend.position="right")

#Non-linear dimensionality reduction - UMAP
monocle_umap_ctrl <- reduce_dimension(monocle_object_p_ctrl)

#scatterplot uMAP
plot_cells(monocle_umap_ctrl, label_groups_by_cluster=FALSE,  color_cells_by = "type")+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Cluster your cells: each cell is assigned not only to a cluster but also to a partition. When you are learning trajectories, each partition will eventually become a separate trajectory. We run cluster_cells()as before.
monocle_cluster_ctrl <- cluster_cells(monocle_umap_ctrl)

#partition based
plot_cells(monocle_cluster_ctrl, color_cells_by = "partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Clusters based
plot_cells(monocle_cluster, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Cell type based
monocle_cluster <- cluster_cells(monocle_umap)
plot_cells(monocle_cluster, color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

jpeg('Integrated_partition.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_cluster, color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

jpeg('Integrated_type.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_cluster, color_cells_by = "type")+
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

jpeg('Integrated_treatment_time.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_cluster, color_cells_by = "treat.time1") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()

############# Count cells ######################


# cluster distribution in each group
#table(Idents(monocle_cluster), monocle_cluster$type)
monocle_cluster$type

group_dist_p  <- table(partitions(monocle_cluster), monocle_cluster$type)
group_dist  <- table(clusters(monocle_cluster), monocle_cluster$type)


group_dist

#Write results into a file
write.table(group_dist,file="cell_distribution_in_clusters.txt", sep='\t',  quote = F,row.names = FALSE)


# What proportion of cells are in each cluster?
prop.table(table(partitions(monocle_cluster)))

table(partitions(monocle_cluster), monocle_cluster$type)
table(clusters(monocle_cluster), monocle_cluster$type)

#proportion of clusters in each group
group_dist_p <- prop.table(table(partitions(monocle_cluster), monocle_cluster$type), margin = 2)


#Write results into a file
write.table(group_dist_p,file="cell_dis_proporton.txt", sep='\t',  quote = F,row.names = FALSE)

#Known markers
markers_ESC <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4") ##Esc
markers_stem <- c("CD146", "CD45", "PDGFR-beta") ##stem markers
markers_ectodem <- c("Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6") ##Ectoderm
markers_endoderm <- c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1") ##Endoderm
markers_mesoderm <- c("Mest" , "Bmp2","Eomes", "Hand1" , "Isl1", "Kdr", "Mesdc1" , "Mesdc2", "Myf5", "Myod1", "Nkx2-5", "T", "Tbx2", "Prdm1", "Tbx1") ##mesoderm
markers_neural_stem <- c("Apoe", "Blbp", "Gfap", "Slc1a3", "Sox2", "Sox9", "Thrsp") ## neural stem cells
markers_microglial <- c("1700112E06Rik", "4632428N05Rik", "Abca9", "Abi3", "Adap2", "AF251705", "Aif1", "Apbb1ip", "Arsb", "Bmp2k", "C1qa", "C1qb", "C1qc", "Ccr5", "Cd14", "Cd37", "Cd53", "Cd68", "Csf1r", "Ctsb", "Ctss", "Cx3cr1", "Cyth4", "Dock2", "Dock8", "Emr1", "Entpd1", "Fcer1g", "Fcgr3", "Fcrls", "Fli1", "Fyb", "Gpr34", "Hexb", "Hmha1", "Hpgd", "Hpgds", "Ikzf1", "Il6ra", "Inpp5d", "Irf8", "Itgam", "Itgb5", "Lair1", "Laptm5", "Lcp2", "Lgmn", "Lpcat2", "Ltc4s", "Ly86", "Lyn", "Mafb", "Mertk", "Mpeg1", "Myo1f", "Ncf1", "Nckap1l", "Olfml3", "P2ry12", "P2ry13", "Pik3ap1", "Pld4", "Pros1", "Ptgs1", "Ptprc", "Rasgrp3", "Rnase4", "Rreb1", "Runx1", "Selplg", "Serinc3", "Siglech", "Sirpa", "Skap2", "Slco2b1", "Tbxas1", "Tgfbr1", "Tgfbr2", "Tmem119", "Trem2", "Tyrobp", "Unc93b1", "Vav1", "Zfp710", "0610040J01Rik")
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
marker_genes <- markers_type1_ganglion  #### change marker genes according to cells we are looking for

plot_cells(monocle_cluster, genes=marker_genes,label_cell_groups=FALSE, show_trajectory_graph=FALSE)


#Find marker genes for each partition
marker_clusters <- top_markers(monocle_cluster, group_cells_by="partition",  reference_cells=1000, cores=8)
write.table(marker_clusters,file="markers_integrated_data_partition.txt", sep='\t',  quote = F,row.names = FALSE)

#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers <- marker_clusters %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=4)) 


plot_cells(monocle_cluster, color_cells_by = "partition") +
  #scale_color_manual(values = cluster) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right") #put the color legend on the right

#Find marker genes for each cluster
marker_clusters1 <- top_markers(monocle_cluster, group_cells_by="cluster",  reference_cells=1000, cores=8)
write.table(marker_clusters1,file="markers_integrated_data_cluster.txt", sep='\t',  quote = F,row.names = FALSE)

#We can rank markers according to pseudo_R2, a specificity metrics 
top_specific_markers_c <- marker_clusters1 %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(10, pseudo_R2)

top_specific_marker_ids_c <- unique(top_specific_markers_c %>% pull(gene_id))

plot_genes_by_group(monocle_cluster,
                    top_specific_marker_ids_c,
                    group_cells_by="cluster",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=4)) 


plot_cells(monocle_cluster, color_cells_by = "cluster") +
  #scale_color_manual(values = cluster) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right") #put the color legend on the right

#Learn the trajectory graph: we will fit a principal graph within each parition using the learn_graph() function:
monocle_traj <- learn_graph(monocle_cluster)

#Trajectory plot
plot_cells(monocle_traj, color_cells_by = "type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#save plot
jpeg('Integrated_trajectory_type.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_traj, color_cells_by = "type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()


#Trajectory plot treatment time wise
plot_cells(monocle_traj, color_cells_by = "treat.time1",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

#Save plot
jpeg('Integrated_trajectory_time.jpg', units="in", width=10, height=8, res=300)
plot_cells(monocle_traj, color_cells_by = "treat.time1",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,label_branch_points=FALSE) +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")
dev.off()


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



plot_genes_by_group(monocle_traj,
                    top_specific_marker_ids1,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=1.5) +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=5), axis.text.y=element_text(angle=0, hjust=1, size=3)) 


#Differential gene expression analysis between partitions
all_genes <- unique(rowData(monocle_traj)$gene_short_name)
all_genes 

all_gene_subset <- monocle_traj[rowData(monocle_traj)$gene_short_name %in% all_genes,]


all_gene_fits_cls <- fit_models(all_gene_subset, model_formula_str = "~partition")
head(all_gene_fits_cls)

#all_gene_fits <- fit_models(all_gene_subset, model_formula_str = "~treat.time1")

all_gene_fit_coefs <- coefficient_table(all_gene_fits_cls)
all_gene_fit_coefs <- coefficient_table(all_gene_fits_cls)

all_gene_fit_coefs
#Note that the table includes one row for each term of each gene's model. We generally don't care about the intercept term 
#??0, so we can easily just extract the time terms:

all_gene_time_terms <- all_gene_fit_coefs %>% filter(term != "(Intercept)")
all_gene_time_terms1 <- all_gene_time_terms %>% filter(status == "OK")

#Extract significant genes based on q-values
sig_genes <- all_gene_time_terms1 %>% filter (q_value < 0.05) %>% select(gene_short_name, term, q_value, estimate)
sig_genes_df  <- as.data.frame(sig_genes)
sig_genes_df
write.table(sig_genes_df,file="all_sig_genes_df.txt", sep='\t',  quote = F,row.names = FALSE)

#Visualize top genes
sig <- c("H2afj", "Malat1", "Prelid2", "Bcar3", "Ifitm1", "Sec61g", "Fn1", "Minos1", "AC160336.1", "Rpl37")

colData( monocle_traj)$cluster <- as.character(partitions( monocle_traj))

#colData(monocle_cluster)$cluster <- as.character(partitions(monocle_cluster))


top_sig_set <- monocle_traj[rowData(monocle_traj)$gene_short_name %in% sig,]
#top_sig_set <- monocle_cluster[rowData(monocle_cluster)$gene_short_name %in% sig,]

plot_genes_violin(top_sig_set, group_cells_by="cluster", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))


################### Create Trajectory for specific cluster   ########################
#order_cells()needs you to specify the root nodes of the trajectory graph. 
monocle_traj_order <- order_cells(monocle_traj)

#Plot trajectory plot for specific cluster
plot_cells(monocle_traj_order, color_cells_by = "pseudotime",
           label_cell_groups=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,graph_label_size=1.5)

#Looks at pattern of specific cell type marker as pseudotime function
#First let's see ESC markers
ESC_cds <- monocle_traj_order[rowData(monocle_traj_order)$gene_short_name %in% markers_ESC,]
ESC_cds

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(ESC_cds, color_cells_by="pseudotime", min_expr=0.5)

#Neural stem cells markers

neural_stem_lineage_cds <- monocle_traj_order[rowData(monocle_traj_order)$gene_short_name %in% markers_neural_stem,]

#The function plot_genes_in_pseudotime() takes a small set of genes and shows you their dynamics as a function of pseudotime:

plot_genes_in_pseudotime(neural_stem_lineage_cds, color_cells_by="pseudotime", min_expr=0.5)

######################################################################

#Choose subset cluster2
cluster2_subset <- choose_cells(monocle_traj_order)


#Now we have a smaller cell_data_set object that contains just the cells from the partition we'd like to drill into. 
#We can use graph_test() to identify genes that are differentially expressed in different subsets of cells from this partition:
pr_graph_test_clus2 <- graph_test(cluster2_subset, neighbor_graph="knn", cores=8)
pr_deg_ids_clus2 <- row.names(subset(pr_graph_test_clus2, morans_I > 0.05 & q_value < 0.05))

pr_graph_test_clus2
pr_deg_ids_clus2
#We will learn more about graph_test() in the differential expression analysis section later. We can take all the genes that vary across this set of cells and group those that have similar patterns of expression into modules:

gene_module_df_clus2 <- find_gene_modules(cluster2_subset[pr_deg_ids_clus2,])


#Plotting these modules' aggregate expression values reveals which cells express which modues.

plot_cells(cluster2_subset, genes=gene_module_df_clus2, 
           show_trajectory_graph=FALSE, 
           label_cell_groups=FALSE)





################
## integrate samples
integ_anchors_ctrls <- FindIntegrationAnchors(object.list = list(Encap_EB_ftr, Encap_EB_2D_ftr, Encap_EB_4D_ftr))
integ_anchors_encap

Integ_data_encap <- IntegrateData(anchorset = integ_anchors_encap, dims = 1:100)
DefaultAssay(Integ_data_encap) <- "integrated"


#Integ_data_cds <- as(as.matrix(GetAssayData(Integ_data, assay = "integrated", slot = "scale.data")), 'sparseMatrix')

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

#UMAP plot
plot_cells(monocle_umap_encap1, label_groups_by_cluster=FALSE,  color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

# clustering
monocle_cluster_encap <- cluster_cells(monocle_umap_encap1)

#UMAP plot based on sample
plot_cells(monocle_cluster_encap , label_groups_by_cluster=FALSE,  color_cells_by = "type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


#UMAP plot based on cluster
plot_cells(monocle_cluster_encap, label_groups_by_cluster=FALSE,  color_cells_by = "cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

