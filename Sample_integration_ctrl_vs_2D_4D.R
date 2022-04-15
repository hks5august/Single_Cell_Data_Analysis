#install.packages('cowplot')

#install.packages("harmony")

remotes::install_github('satijalab/seurat-wrappers')
library(SeuratWrappers)

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

seurat_ob_S1

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S1[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S1, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(seurat_ob_S1@meta.data, 5)

jpeg('QC_sample1.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

data2 <- Read10X(data.dir = "E2_S2/filtered_feature_bc_matrix/")

seurat_ob_S2 <- CreateSeuratObject(counts = data2, project = "E2_S2", min.cells = 3, min.features = 200)


head(seurat_ob_S2@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S2[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S2, pattern = "^MT-")

seurat_ob_S2

# Show QC metrics for the first 5 cells
head(seurat_ob_S2@meta.data, 5)

jpeg('QC_EB_2D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


data3 <- Read10X(data.dir = "E3_S3/filtered_feature_bc_matrix/")

seurat_ob_S3 <- CreateSeuratObject(counts = data3, project = "E3_S3", min.cells = 3, min.features = 200)


head(seurat_ob_S3@meta.data, 5)

# To check the content of dead cells, we can check mitocondrial percentage; The [[ operator can add columns to object metadata. This is a great place to stash QC stats
seurat_ob_S3[["percent.mt"]] <- PercentageFeatureSet(seurat_ob_S3, pattern = "^MT-")

seurat_ob_S3

# Show QC metrics for the first 5 cells
head(seurat_ob_S3@meta.data, 5)

jpeg('QC_EB_4D.jpg', units="in", width=10, height=8, res=300)
# Visualize QC metrics as a violin plot
VlnPlot(seurat_ob_S3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()





seurat_ob_S1$type <- "EB_Ctrl"
seurat_ob_S2$type <- "EB_2D"
seurat_ob_S3$type <- "EB_4D"


# filter genes with low values (Here, we are removing empty cells/droplet)
EB <- subset(seurat_ob_S1, subset = nFeature_RNA > 200)
EB_2D <- subset(seurat_ob_S2, subset = nFeature_RNA > 200)
EB_4D <- subset(seurat_ob_S3, subset = nFeature_RNA > 200)


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


## integrate samples
integ_anchors <- FindIntegrationAnchors(object.list = list( EB_2D_ftr, EB_ftr, EB_4D_ftr), dims = 1:20)
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

jpeg('UMAP_plot_EB_2D_4D_comb.jpg', units="in", width=10, height=8, res=300)
p1
dev.off()

jpeg('UMAP_plot_EB_2D_4D_plot_cluster.jpg', units="in", width=10, height=8, res=300)
p2
dev.off()

jpeg('UMAP_plot_cluster_comb.jpg', units="in", width=10, height=8, res=300)
plot_grid(p1, p2)

dev.off()

#visualize both samples side by side

jpeg('UMAP_plot_EB_EB_2D_comb2_side_by_side.jpg', units="in", width=10, height=8, res=300)
DimPlot(Integ_data_clus, reduction = "umap", split.by = "type")
dev.off()


DefaultAssay(Integ_data_clus) <- "RNA"


Integ_data_clus1 <- RenameIdents(Integ_data_clus, `0` = "cell_type1", `1` = "cell_type2", `2` = "cell_type3", 
                                 `3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9", `9` = "cell_type10")

jpeg('UMAP_plot_EB_vs_2D_with_diff_cell_types.jpg', units="in", width=10, height=8, res=300)

DimPlot(Integ_data_clus1, label = TRUE)

dev.off()

jpeg('UMAP_plot_EB_vs_2D_cell_types_overlapping.jpg', units="in", width=10, height=8, res=300)

DimPlot(Integ_data_clus1, label = TRUE, group.by = "type")
dev.off()

# Visualization
pp1 <- DimPlot(Integ_data_clus1, reduction = "umap", group.by = "type")
pp1

pp2 <- DimPlot(Integ_data_clus1, reduction = "umap", label = TRUE)
pp2

jpeg('UMAP_plot_by_type_cell_type_side_by_side.jpg', units="in", width=10, height=8, res=300)
plot_grid(pp1, pp2)

dev.off()



jpeg('UMAP_plot_EB_EB_2D_comb2_side_by_side_with_cell_type.jpg', units="in", width=10, height=8, res=300)
DimPlot(Integ_data_clus1, reduction = "umap", split.by = "type")
dev.off()


# cluster distribution in each group
table(Idents(Integ_data_clus1), Integ_data_clus1$type)

group_dist = as.data.frame(table(Idents(Integ_data_clus1), Integ_data_clus1$type))


#Write results into a file
write.table(group_dist,file="EB_ctrl_vs_EB_2D_4D_cell_dis.txt", sep='\t',  quote = F,row.names = FALSE)


# What proportion of cells are in each cluster?
prop.table(table(Idents(Integ_data_clus1)))

#proportion of clusters in each group
round(prop.table(table(Idents(Integ_data_clus1), Integ_data_clus1$type), margin = 2),2)


#Write results into a file
write.table(group_dist,file="EB_ctrl_vs_EB_2D_4D_cell_dis_proporton.txt", sep='\t',  quote = F,row.names = FALSE)














#find conserved markers: calculated the genes that are conserved markers irrespective of stimulation condition in cluster 0 (large cluster).
clus0_markers <- FindConservedMarkers(Integ_data_clus, ident.1 = 0, grouping.var = "type", verbose = FALSE)
head(clus0_markers, 10)

dim(clus0_markers)

clus0_markers_df <- as.data.frame(clus0_markers)
#Write results into a file
write.table(clus0_markers_df,file="cluster0_markers.txt", sep='\t',  quote = F,row.names = T)


######### 

jpeg('Marker_feature_plot_ED.jpg', units="in", width=10, height=8, res=300)

#We can explore marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(Integ_data_clus1, features = c("Hspd1", "Ncl", "Supt16", "Set", "Nasp","Hnrnpf", "mt-Co1", "Gm49359", "mt-Cytb",  "Gm42418", "Gm47283", "Gm26917"), min.cutoff = "q9")
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

pdf('plot_genes_marker_genes.pdf', width=20, height=16)

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
                                 `3` = "cell_type4", `4` = "cell_type5", `5` = "cell_type6", `6` = "cell_type7", `7` = "cell_type8", `8` = "cell_type9", `9` = "cell_type10")

DimPlot(Integ_data_clus1, label = TRUE)
DimPlot(Integ_data_clus1, label = TRUE, group.by = "type")

# Visualization
pp1 <- DimPlot(Integ_data_clus1, reduction = "umap", group.by = "type")
pp2 <- DimPlot(Integ_data_clus1, reduction = "umap", label = TRUE)
plot_grid(pp1, pp2)

# Dotplot 

Idents(Integ_data_clus1) <- factor(Idents(Integ_data_clus1), levels = c("cell_type1", "cell_type2", 
                                                                        "cell_type3", "cell_type4", "cell_type5", "cell_type6", "cell_type7", "cell_type8", "cell_type9", "cell_type10"))
#markers.to.plot <- c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4") ##Esc
#markers.to.plot <- c("CD146", "CD45", "PDGFR-beta") ##stem markers
#markers.to.plot <- c("Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6") ##Ectoderm
#markers.to.plot <- c("Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1") ##Endoderm
#markers.to.plot <- c("Mest" , "Bmp2","Eomes", "Hand1" , "Isl1", "Kdr", "Mesdc1" , "Mesdc2", "Myf5", "Myod1", "Nkx2-5", "T", "Tbx2", "Prdm1", "Tbx1") ##mesoderm
#markers.to.plot <- c("Apoe", "Blbp", "Gfap", "Slc1a3", "Sox2", "Sox9", "Thrsp") ## neural stem cells
#markers.to.plot <- c("1700112E06Rik", "4632428N05Rik", "Abca9", "Abi3", "Adap2", "AF251705", "Aif1", "Apbb1ip", "Arsb", "Bmp2k", "C1qa", "C1qb", "C1qc", "Ccr5", "Cd14", "Cd37", "Cd53", "Cd68", "Csf1r", "Ctsb", "Ctss", "Cx3cr1", "Cyth4", "Dock2", "Dock8", "Emr1", "Entpd1", "Fcer1g", "Fcgr3", "Fcrls", "Fli1", "Fyb", "Gpr34", "Hexb", "Hmha1", "Hpgd", "Hpgds", "Ikzf1", "Il6ra", "Inpp5d", "Irf8", "Itgam", "Itgb5", "Lair1", "Laptm5", "Lcp2", "Lgmn", "Lpcat2", "Ltc4s", "Ly86", "Lyn", "Mafb", "Mertk", "Mpeg1", "Myo1f", "Ncf1", "Nckap1l", "Olfml3", "P2ry12", "P2ry13", "Pik3ap1", "Pld4", "Pros1", "Ptgs1", "Ptprc", "Rasgrp3", "Rnase4", "Rreb1", "Runx1", "Selplg", "Serinc3", "Siglech", "Sirpa", "Skap2", "Slco2b1", "Tbxas1", "Tgfbr1", "Tgfbr2", "Tmem119", "Trem2", "Tyrobp", "Unc93b1", "Vav1", "Zfp710", "0610040J01Rik")
#markers.to.plot <- c("Sox11", "Tbr2") #Early intermediate precursor cell
#markers.to.plot <- c("Acvr2a", "Fgf10", "Fgf3", "Fgf8", "Hes5", "Neurog1", "Notch1", "Six1") #Early_neuroblast
#markers.to.plot <- c("Bcar3", "Celsr3", "Chrna4", "Cpox", "Dlgap3", "Dok7", "Egfr", "Kcnh3", "Lrrc3b", "Miat", "Npy2r", "Nr2e1", "Ostf1", "Ptms", "Rasd2", "Rgs16", "Rprm", "Sema5b", "Slc12a4", "St3gal1", "Tpbg", "Unc5a", "Wscd2", "2010300C02Rik")#Interneuron-selective cell
#markers.to.plot <- c("1700086L19Rik" , "A930011G23Rik" , "Adamts16" , "Arhgdig" , "Arhgef15" , "Arl4a" , "Ass1" , "B2m" , "Bace2" , "Carhsp1" , "Ccdc109b" , "Chodl" , "Clic1" , "Cnih3" , "Col11a1" , "Cort" , "Ctxn2" , "D430019H16Rik" , "Dbpht2" , "Dpy19l1" , "Egfl7" , "Fam46a" , "Fndc1" , "Gabra2" , "Gabrg1" , "Gm11744" , "Gpr126" , "Gpr88" , "Hey1" , "Htr1a" , "Kcnmb4" , "N28178" , "Ndrg1" , "Ndst4" , "Nnat" , "Nos1" , "Nptx1" , "Ntn1" , "Opn3" , "Oxtr" , "Patl2" , "Pnma3" , "Pou3f2" , "Prex1" , "Ptn" , "Ptpru" , "Rasgef1b" , "Rbp1" , "S100a10" , "Sc5d" , "Sfrp1" , "Slc5a5" , "Slc5a7" , "Slc7a3" , "Sst" , "St6galnac2" , "Sv2b" , "Tppp3" , "Tspan17" , "Wnt2" , "Wwc1" , "1700058G18Rik") #Long-projecting GABAergic cell
#markers.to.plot <- c("Abcb1a", "Abcg2", "Acvrl1", "Apcdd1", "Car4", "Ccdc141", "CD34", "Cd93", "Cdh5", "Cgnl1", "Cldn5", "Clic5", "Ctla2a", "Cxcl12", "Cyyr1", "Egfl7", "Eltd1", "Emcn", "Eng", "Epas1", "Erg", "Esam", "Ets1", "Fli1", "Flt1", "Fn1", "Foxq1", "Fzd6", "Gpr116", "Hmcn1", "Hspb1", "Hspg2", "Ifitm3", "Ifltd1", "Igfbp7", "Itga1", "Itga4", "Itm2a", "Kank3", "Kdr", "Kitl", "Klf2", "Lama4", "Lsr", "Ly6a", "Ly6c1", "Ly6c2", "Ly75", "Mecom", "Mfsd2a", "Nos3", "Nostrin", "Palmd", "Paqr5", "Pcp4l1", "Pecam1", "Pglyrp1", "Pltp", "Podxl", "Ptprb", "Ramp2", "Rasgrp3", "Rassf9", "Rbpms", "Rgs5", "Rhoj", "Sdpr", "Sema3c", "Sgms1", "Slc16a1", "Slc22a8", "Slc2a1", "Slc40a1", "Slc7a1", "Slc7a5", "Slc9a3r2", "Slco1a4", "Slco1c1", "Slco2b1", "Sox17", "Sparc", "Srgn", "St3gal6", "Tek", "Tm4sf1", "Vwf", "Wfdc1", "Wwtr1", "Zfp366", "9430020K01Rik") ###endothelial cells
#markers.to.plot <- c("6330403K07Rik", "A030009H04Rik", "Agap2", "Ahi1", "Arf3", "Atp1a3", "Atp1b1", "Bex2", "Camk2a", "Camk2b", "Camk2n1", "Celf4", "Chd3", "Cnih2", "Cplx2", "Ctxn1", "Dynll2", "Fbxl16", "Gap43", "Gng13", "Hap1", "Hpcal4", "Kif5a", "Kif5c", "Ly6h", "Map1a", "Map1b", "Map2", "Myt1l", "Ncdn", "Nefl", "Nrgn", "Nrsn1", "Pcp4", "Pcsk1n", "Rab3c", "Rasgef1a", "Rph3a", "Rtn1", "S100b", "Scg2", "Sepw1", "Sez6", "Snap25", "Snrpn", "Stmn2", "Stmn3", "Syt1", "Ttc3", "Unc13c") ##neurons
#markers.to.plot <- c("Ccnd2", "Crmp1", "Dbn1", "Dchs1", "Dlx1", "Dlx2", "Dlx5", "Dlx6", "Evf1", "Foxc1", "Fxyd6", "Klf5a", "Olfm1", "Pfn2", "Pou3f3", "Sox5", "Trn1", "Btg1", "Celf4", "Dcx", "Meis2", "Stmn2") #neuroblast
#markers.to.plot <- c("2810417H13Rik", "Cdk1", "Cenpf", "Hist1h2ap", "Hmgb2", "Top2a")##Late activated neural stem cell
#markers.to.plot <- c("Eya2", "Fzd10", "Isl1", "Jag2", "Neurod1", "Neurog2") #Late neuroblast
#markers.to.plot <- c("Csf3", "Fbxo2", "Fxyd1", "Gfap", "Id3", "Thbs4") #Activated neural stem cell
#markers.to.plot <- c("Neurod1", "Neurog1") #precusor cells
#markers.to.plot <- c("Omp") #Mature olfactory sensory neuron
#markers.to.plot <- c("Ascl1") #Progenitor cell
#markers.to.plot <- c("Gap43", "Gng8") #immature olfactory sensory neuron

markers.to.plot <- c( "3300002A11Rik" , "4930529M08Rik", "9330101J02Rik", "Acta2", "Adamts20", "Ak7", "Ak9", "Armc3", "Ccdc108", "Ccdc114", "Ccdc146", "Ccdc153", "Ccdc162", "Ccdc170", "Ccdc60", "Clu", "Daw1", "Dnah10", "Dnah11", "Dnah12", "Dnah3", "Dnah5", "Dnah6", "E230008N13Rik", "Efhb", "Enkur", "Fhad1", "Foxj1", "Gja1", "Gm973", "Hydin", "Ifltd1", "Iqck", "Kif6", "Kif9", "Lrguk", "Lrrc48", "Lrriq1", "Mia", "Nek11", "Rarres2", "Rgs22", "S100b", "Sox9", "Spag16", "Spag17", "Spata17", "Spef2", "Tm4sf1", "Tmem212", "Tmem232", "Ttc21a", "Ttc29", "Ttll8", "Vwa3a", "Wdr52", "Wdr63", "Wdr96", "Zbbx", "1700007G11Rik") ##Ependymal cells


jpeg('Dotplot_Ependymal_markers.jpg', units="in", width=10, height=8, res=300)

DotPlot(Integ_data_clus1, features = markers.to.plot, cols = c("blue", "red", "green"), dot.scale = 8, split.by = "type") + RotatedAxis()

dev.off()



jpeg('Marker_feature_plot.jpg', units="in", width=10, height=8, res=300)

#We can explore marker genes for each cluster and use them to annotate our clusters as specific cell types.
FeaturePlot(Integ_data_clus, features = c("Hspd1", "Ncl", "Supt16", "Set", "Nasp","Hnrnpf", "mt-Co1", "Gm49359", "mt-Cytb",  "Gm42418", "Gm47283", "Gm26917"), min.cutoff = "q9")
dev.off()



jpeg('plot_genes_marker_genes.jpg', units="in", width=20, height=16, res=300)

FeaturePlot(Integ_data_clus, features = c("Alcam", "Isl1", "N-cadherin", "Nanog", "Oct4", "Hoxb1", "Lhx5", "Nes", "Neurod1", "Otx1", "Pax6", "Mest", "Foxa2", "Sox17",  "Gata4", "Gata5", "Gata6", "Onecut1", "Sox10"), min.cutoff = "q9")


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