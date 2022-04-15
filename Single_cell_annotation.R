
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("celldex")



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SingleR")
BiocManager::install("scran")
BiocManager::install("scRNAseq")
BiocManager::install("scuttle")
BiocManager::install("GSEABase")
BiocManager::install("AUCell")

library(scran)
library(GSEABase)
library(AUCell)
library(SingleR)
library(celldex)
library(scuttle)
ref <- MouseRNAseqData()
ref

colData(EB_4D_cds_clus)$cluster <- as.character(clusters(EB_4D_cds_clus))
colData(Encap_4D_cds_clus)$cluster <- as.character(clusters(Encap_4D_cds_clus))

colData(EB_2D_cds_clus)$cluster <- as.character(clusters(EB_2D_cds_clus))
colData(Encap_2D_cds_clus)$cluster <- as.character(clusters(Encap_2D_cds_clus))


library(SingleR)
pred <- SingleR(test=EB_4D_cds_clus, ref=ref, labels=ref$label.fine)

pred

table(pred$labels)

plotScoreHeatmap(pred)

tab <- table(Assigned=pred$pruned.labels, Cluster=clusters(EB_4D_cds_clus))

tab
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))




pred1_m <- SingleR(test=monocle_EB_4D_cluster, ref=ref, assay.type.test=1, labels=ref$label.fine)

pred1_m

table(pred1_m$labels)
table(pred1_m$first.labels)
table(pred1_m$pruned.labels)
plotScoreHeatmap(pred1_m)

tab1_1_m <- table(Assigned=pred1_m$first.labels, Cluster=clusters(monocle_EB_4D_cluster))

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab1_1_m+10), color=colorRampPalette(c("white", "blue"))(50))

tab1_1_m <- table(Assigned=pred1_m$labels, Cluster=clusters(monocle_EB_4D_cluster))


tab1_1_m



colData(monocle_EB_4D_cluster)$cell_type <- as.character(pred1_m$pruned.labels)

plot_cells(monocle_EB_4D_cluster, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(monocle_EB_4D_cluster, color_cells_by="partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


plot_cells(monocle_EB_4D_cluster, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")



#### EB 2D #####
pred1_2D_EB <- SingleR(test=monocle_EB_2D_cluster, ref=ref, assay.type.test=1, labels=ref$label.fine)

pred1_2D_EB 

table(pred1_2D_EB$labels)
table(pred1_2D_EB$first.labels)
table(pred1_2D_EB$pruned.labels)
plotScoreHeatmap(pred1_2D_EB)

tab1_EB_2D <- table(Assigned=pred1_2D_EB$labels, Cluster=clusters(monocle_EB_2D_cluster))
tab1_EB_2D 

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab1_EB_2D+10), color=colorRampPalette(c("white", "blue"))(50))

tab1_EB_2D <- table(Assigned=pred1_2D_EB$labels, Cluster=clusters(monocle_EB_2D_cluster))


tab1_EB_2D



colData(monocle_EB_2D_cluster)$cell_type <- as.character(pred1_2D_EB$labels)

plot_cells(monocle_EB_2D_cluster, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(monocle_EB_2D_cluster, color_cells_by="partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


plot_cells(monocle_EB_2D_cluster, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")



#Encap


#pred1 <- SingleR(test=Encap_4D_cds_clus, ref=ref, labels=ref$label.main)


pred1 <- SingleR(test=Encap_4D_cds_clus, ref=ref, assay.type.test=1, labels=ref$label.fine)



table(colnames(pred1$scores))
unique(pred1$scores)
plotScoreHeatmap(pred1)

tab1_1 <- table(Assigned=pred1$labels, Cluster=clusters(Encap_4D_cds_clus))
tab1_1

tab1 <- table(Assigned=pred1$pruned.labels, Cluster=clusters(Encap_4D_cds_clus))

tab1

# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab1+10), color=colorRampPalette(c("white", "blue"))(50))




pred_2D_en <- SingleR(test=monocle_Encap_2D_cluster, ref=ref, assay.type.test=1, labels=ref$label.fine)

pred_2D_en 

table(pred_2D_en$labels)
table(pred_2D_en$first.labels)
table(pred_2D_en$pruned.labels)
plotScoreHeatmap(pred_2D_en )

tab1_2D_en <- table(Assigned=pred_2D_en$labels, Cluster=clusters(monocle_Encap_2D_cluster))
tab1_2D_en 
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
library(pheatmap)
pheatmap(log2(tab1_2D_en+10), color=colorRampPalette(c("white", "blue"))(50))

tab1_2D_en<- table(Assigned=pred_2D_en$labels, Cluster=clusters(monocle_EB_2D_cluster))


colData(monocle_Encap_2D_cluster)$cell_type <- as.character(pred_2D_en$labels)

plot_cells(monocle_Encap_2D_cluster, color_cells_by="cluster") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")

plot_cells(monocle_Encap_2D_cluster, color_cells_by="partition") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")


plot_cells(monocle_Encap_2D_cluster, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")




plot_cells(monocle_Encap_4D_cluster, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")




plot_cells(monocle_Encap_4D_cluster, color_cells_by="cell_type") +
  theme(legend.text=element_text(size=6)) + #set the size of the text
  theme(legend.position="right")





library(scran)

library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel,   id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))


library(org.Mm.eg.db)
rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")

unfiltered <- sce.zeisel

stats <- perCellQCMetrics(sce.zeisel, subsets=list(  Mt=rowData(sce.zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", "subsets_Mt_percent"))
sce.zeisel <- sce.zeisel[,!qc$discard]

colData(unfiltered) <- cbind(colData(unfiltered), stats)
unfiltered$discard <- qc$discard


library(scran)
set.seed(1000)
clusters <- quickCluster(sce.zeisel)
sce.zeisel <- computeSumFactors(sce.zeisel, cluster=clusters) 
sce.zeisel <- logNormCounts(sce.zeisel)


dec.zeisel <- modelGeneVarWithSpikes(sce.zeisel, "ERCC")
top.hvgs <- getTopHVGs(dec.zeisel, prop=0.1)


library(BiocSingular)

sce.zeisel <- denoisePCA(sce.zeisel, technical=dec.zeisel, subset.row=top.hvgs)
sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA")


wilcox.z <- pairwiseWilcox(sce.zeisel, sce.zeisel$level1class, lfc=1, direction="up")

markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, pairwise=FALSE, n=50)
lengths(markers.z)



library(GSEABase)
all.sets <- lapply(names(markers.z), function(x) {GeneSet(markers.z[[x]], setName=x)})
all.sets <- GeneSetCollection(all.sets)
all.sets 
library(AUCell)
rankings <- AUCell_buildRankings(counts(monocle_cluster_EB), plotStats=FALSE, verbose=FALSE)
cell.aucs <- AUCell_calcAUC(all.sets, rankings)

par(mfrow=c(3,3))
AUCell_exploreThresholds(cell.aucs, plotHist=TRUE, assign=TRUE) 

results <- t(assay(cell.aucs))
head(results)

new.labels <- colnames(results)[max.col(results)]
new.labels

tab <- table(new.labels, monocle_cluster_EB$cluster)
tab

library(pheatmap)

fullheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)
fullheat

subheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)

gridExtra::grid.arrange(fullheat[[4]], subheat[[4]])


# Downloading the signatures and caching them locally.
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
scsig.path <- bfcrpath(bfc, file.path("http://software.broadinstitute.org", "gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"))
scsigs <- getGmt(scsig.path)

#to get all dataset names
names(scsigs)

list_data <- as.data.frame(names(scsigs))
list_data
write.table(list_data,file="single_cell_dataset_list.txt", sep='\t',  quote = F,row.names = T)

# Restricting to the subset of gene sets:
scsigs.neurons <- scsigs[grep("neuron", names(scsigs))]
sub.aucs <- AUCell_calcAUC(scsigs.neurons, rankings)
sub.results <- t(assay(sub.aucs))
sub.labels <- colnames(sub.results)[max.col(sub.results)]
tab <- table(sub.labels, sce.muraro$label)
subheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)

gridExtra::grid.arrange(fullheat[[4]], subheat[[4]])



markers.mam <- scoreMarkers(sce.mam, lfc=1)

chosen <- "2"
cur.markers <- markers.mam[[chosen]]
is.de <- order(cur.markers$median.logFC.cohen, decreasing=TRUE)[1:100]
cur.markers[is.de,1:4]










EB_4D_cds_clus


wilcox.z <- pairwiseWilcox(sce.zeisel, sce.zeisel$level1class, lfc=1, direction="up")

wilcox.z <- pairwiseWilcox(EB_4D_cds_clus, EB_4D_cds_clus$cluster, lfc=1, direction="up")

markers.z <- getTopMarkers(wilcox.z$statistics, wilcox.z$pairs, pairwise=FALSE, n=50)
lengths(markers.z)


###############
### cell type annotation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GSVA")

install.packages("clustermole")
library(clustermole)

Type <- c()
clusters <- c()
for (i in seq(1,2, by=1)) {
  markers.filt <- as.data.frame(top_specific_markers_2d %>% filter(cell_group == as.character(i)) %>% top_n(n = 25))
  my_overlaps <- clustermole_overlaps(genes = markers.filt$gene, species = "mm")
  type <- as.character(my_overlaps[1,5])
  Type <- c(Type, type)
  clusters <- c(clusters, i)
}


clust_map <- as.data.frame(cbind(clusters,Type))
head(clust_map)



#############
#install.packages("remotes")
remotes::install_github("Irrationone/cellassign")

library(cellassign)
