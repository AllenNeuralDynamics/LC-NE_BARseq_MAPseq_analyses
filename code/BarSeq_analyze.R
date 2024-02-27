####### Analyze transcriptomic data from pre-processed rds files #######
# script order 03
# barseq-r4 environment to run

setwd('~/capsule/results')

#perform primary analyses to  the BarSeq data (old params  k_umap=100, k_cluster=50)
analyze_barseq = function(barseq, output_name, n_pca=30, k_umap=100, k_cluster=50) {
  output_dir = file.path("analysis", make.names(output_name))
  dir.create(output_dir, showWarnings = FALSE, recursive=TRUE)
  
  barseq = scater::runPCA(barseq, ncomponents=n_pca, exprs_values="logcounts")
  barseq = scater::runUMAP(barseq, dimred="PCA", n_neighbors=k_umap)
  write_csv(as.data.frame(reducedDim(barseq, "UMAP")), file.path(output_dir, "umap.csv"))
  
  clusters = cluster_data(barseq, k_cluster, "jaccard")
  write_csv(clusters, file.path(output_dir, "cluster.csv")) 
  
  # create default cluster annotation
  labels = as.factor(clusters$label)
  cluster_annotation = tibble(cluster_id = levels(labels)) %>%
    mutate(cluster_name = paste0(output_name, "_", cluster_id))
  write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
  return(list(barseq,clusters))
}

cluster_data = function(barseq, k=10, type="rank") {
  # type -> "rank" (preservation of neighbor rank),
  # "number" (number of shared neighbors), "jaccard"
  snn = scran::buildSNNGraph(barseq, k=k, use.dimred = "PCA", type=type,
                             BNPARAM = BiocNeighbors::AnnoyParam(),
                             BPPARAM = BiocParallel::MulticoreParam(12))
  #label = leiden::leiden(snn)
  #label = igraph::cluster_walktrap(g)$membership
  label = igraph::cluster_louvain(snn)$membership
  result = tibble(sample = colnames(barseq), label)
  # only for walktrap
  #table(igraph::cut_at(community.walktrap, n=5))
  return(result)
}

# #run all clustering methods concurrently
# cluster_data = function(barseq, k=10, type="rank") {
#   snn = scran::buildSNNGraph(barseq, k=k, use.dimred = "PCA", type=type,
#                              BNPARAM = BiocNeighbors::AnnoyParam(),
#                              BPPARAM = BiocParallel::MulticoreParam(12))
#   
#   label_leiden = leiden::leiden(snn)
#   label_walktrap = igraph::cluster_walktrap(snn)$membership
#   label_louvain = igraph::cluster_louvain(snn)$membership
#   
#   result_leiden = tibble(sample = colnames(barseq), label = label_leiden)
#   result_walktrap = tibble(sample = colnames(barseq), label = label_walktrap)
#   result_louvain = tibble(sample = colnames(barseq), label = label_louvain)
#   
#   list(leiden = result_leiden, walktrap = result_walktrap, louvain = result_louvain)
# }

#load the relevant data files
barseq <- readRDS('~/capsule/data/filt_neurons-fullbrain_cpm_log.rds') #keeps only cells with min_genes=5, min_counts=20

#subset only barcoded cells for clustering analyses on the local machine
barseq <- barseq[, barseq@colData$barcode == 1]

v<-analyze_barseq(barseq, "barseq_output")
new_barseq <- v[[1]]
clusters <- v[[2]]

#extract the UMAP coordinates, create plotting data frame and plot UMAP
umap_data <- reducedDim(new_barseq, "UMAP")
plot_data <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

#add gene expression levels to the data frame and create gene expression plots
genes <- c("Dbh", "Slc17a7", "Gad1")
for (gene in genes) {
  plot_data[[gene]] <- logcounts(new_barseq)[gene, ]
}
png("BarSeq_louvain_clust.png")
plots <- list(p1)
for (gene in genes) {
  p <- ggplot(plot_data, aes_string(x = "UMAP1", y = "UMAP2", color = gene)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(title = paste("UMAP plot of", gene, "expression"), x = "UMAP1", y = "UMAP2", color = "Expression")
  plots[[length(plots) + 1]] <- p
}
grid.arrange(grobs = plots, ncol = 2)
dev.off()

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LC_barseq)) # should be true
colData(LC_barseq)$leiden_cluster <- as.factor(clusters[["label"]])

#subset cells belonging to LC cluster and cluster them again
LC_barseq <- barseq[, barseq@colData$leiden_cluster == 2]

v<-analyze_barseq(LC_barseq, "LC_output")
new_barseq <- v[[1]]
clusters <- v[[2]]

#extract the UMAP coordinates, create plotting data frame and plot UMAP
umap_data <- reducedDim(new_barseq, "UMAP")
plot_data <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
p1 <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal() +
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

#add gene expression levels to the data frame and create gene expression plots
genes <- c("Dbh", "Th", "Slc18a2")
for (gene in genes) {
  plot_data[[gene]] <- logcounts(new_barseq)[gene, ]
}
plots <- list(p1)
for (gene in genes) {
  p <- ggplot(plot_data, aes_string(x = "UMAP1", y = "UMAP2", color = gene)) +
    geom_point() +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(title = paste("UMAP plot of", gene, "expression"), x = "UMAP1", y = "UMAP2", color = "Expression")
  plots[[length(plots) + 1]] <- p
}
grid.arrange(grobs = plots, ncol = 2)
quartz.save("BarSeq_LC_louvain_clust.png", type = "png", dpi = 300, bg = "white")
dev.off()

#find DE genes for LC clusters identified
colData(new_barseq)$cluster <- factor(clusters[["label"]])

# Compute size factors to account for differences in library sizes
new_barseq <- scran::computeSumFactors(new_barseq, cluster=colData(new_barseq)$cluster)

# Estimate the log-fold changes and adjusted p-values for each gene in each cluster
diff_expr_res <- scran::findMarkers(new_barseq, colData(new_barseq)$cluster)

# 'diff_expr_res' is a list of data frames, one for each cluster, with the log-fold changes and adjusted p-values
# You can access the results for a specific cluster like this:
cluster2_results <- diff_expr_res[["2"]]

#generate a Seurat object out of barseq data
rownames(LC_barseq) <- make.unique(rownames(LC_barseq))
LC_barseq_seurat <- as.Seurat(LC_barseq, counts = "counts", data = "logcounts")
LC_barseq_seurat@assays$originalexp@counts
LC_barseq_seurat@assays$originalexp@data
VlnPlot(LC_barseq_seurat, features = c("nFeature_originalexp", "nCount_originalexp"))
quartz.save("BarSeq_nFeat_nCount_violins.png", type = "png", dpi = 300, bg = "white")
dev.off()

#load scRNAseq data to which BarSeq data is to be compared to
#LC <- readRDS('/Users/polina.kosillo/Seurat_projects/LC_integrate_5/logCPM_integrated_LC_NA_5datasets.rds')
LC<- readRDS('/Users/polina.kosillo/Seurat_projects/LC-NA_large_dataset_integration/integrated_LC_NA_raw_CPM.rds')
#LC<- readRDS('/Users/polina.kosillo/Seurat_projects/LC-NA_large_dataset_integration/harmonized_merged_raw_CPM.rds')

Idents(LC) = 'integrated_snn_res.1'

#label transfer, based on Seurat V3
#use all barseq sampled genes as markers
all_genes <- rownames(LC_barseq_seurat)
anchors <- FindTransferAnchors(reference = LC, query = LC_barseq_seurat, features = all_genes)
#label transfer
predictions <- TransferData(anchorset = anchors, refdata = LC$integrated_snn_res.0.1, k.weight=22)
LC_barseq_seurat <- AddMetaData(LC_barseq_seurat, metadata = predictions)
#check which predicted clusters cells fall into
LC_barseq_seurat$prediction.match <- LC_barseq_seurat$predicted.id
table(LC_barseq_seurat$prediction.match)

#SEURAT V4 map quiery implementation
LC <- RunUMAP(LC, reduction = "pca", dims = 1:20, return.model=TRUE)
DimPlot(LC, reduction = "umap", label = TRUE, 
        label.size = 5) + ggtitle('Integrated clusters UMAP, dims=20, res=0.1')

LC_barseq_seurat <- ScaleData(LC_barseq_seurat)
LC_barseq_seurat <- RunPCA(LC_barseq_seurat, features = rownames(LC_barseq_seurat))
LC_barseq_seurat <- RunUMAP(object = LC_barseq_seurat, reduction = "pca", dims = 1:20, return.model=TRUE)
DimPlot(LC_barseq_seurat, reduction = "umap", label = TRUE,
        label.size = 5) + ggtitle('BarSeq clusters UMAP, dims=20, res=0.1')

all_genes <- rownames(LC_barseq_seurat)
anchors <- FindTransferAnchors(reference = LC, query = LC_barseq_seurat, features = all_genes, 
                               normalization.method = 'LogNormalize', reference.reduction = 'pca')
# LC_barseq_seurat <- MapQuery(anchorset = anchors, reference = LC, query = LC_barseq_seurat, refdata = LC$seurat_clusters, 
#                              reference.reduction = 'pca', reduction.model = 'umap') #there is a bug in this step, perform computations one by one
LC_barseq_seurat <- TransferData(
  anchorset = anchors, 
  reference = LC,
  query = LC_barseq_seurat,
  refdata = LC$integrated_snn_res.0.5,
  k.weight=16
)
LC_barseq_seurat <- IntegrateEmbeddings(
  anchorset = anchors,
  reference = LC,
  query = LC_barseq_seurat, 
  new.reduction.name = "ref.spca"
)
LC_barseq_seurat <- ProjectUMAP(
  query = LC_barseq_seurat, 
  query.reduction = "ref.spca", 
  reference = LC, 
  reference.reduction = "pca", 
  reduction.model = "umap"
)

table(LC_barseq_seurat$predicted.id)



#########################MetaNeighbour#################################
barseq <- LC_barseq_seurat
LC <- LC

type_barseq <- barseq$BarSeq_subclass
type_LC <- LC$integrated_snn_res.0.6

#metaneighbor
new_colData = data.frame(
  study_id = rep(c('barseq','LC'), c(ncol(barseq), ncol(LC))),
  cell_type = c(as.character(type_barseq), as.character(type_LC))
)

common_genes <- intersect(rownames(barseq), rownames(LC))
subset_barseq <- subset(barseq, features = common_genes)
subset_LC <- subset(LC, features = common_genes)

combined_data <- SingleCellExperiment(
  Matrix(cbind(subset_barseq@assays$originalexp@counts,subset_LC@assays$RNA@counts), sparse = TRUE),
  colData = new_colData
)

# use all genes for metaneighbor
var_genes = common_genes
celltype_NV = MetaNeighborUS(var_genes = var_genes,
                             dat = combined_data,
                             study_id = combined_data$study_id,
                             cell_type = combined_data$cell_type,
                             fast_version = TRUE)
cols = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlBu"))(100))
breaks = seq(0, 1, length=101)
sidebarcolors = rownames(celltype_NV)
sidebarcolors[grepl("^barseq",sidebarcolors)]="coral"
sidebarcolors[grepl("^LC",sidebarcolors)]="aquamarine3"
celltype_NV[is.nan(celltype_NV)] = 0

gplots::heatmap.2(celltype_NV,
                  margins=c(8,8),
                  keysize=1,
                  key.xlab="AUROC",
                  key.title="NULL",
                  trace = "none",
                  density.info = "none",
                  col = cols,
                  breaks = breaks,
                  offsetRow=0.1,
                  offsetCol=0.1,
                  cexRow = 1.5,
                  cexCol = 1.5,
                  RowSideColors = sidebarcolors,
                  ColSideColors=sidebarcolors,
                  hclustfun = function(x) hclust(x, method="ward.D")
)
par(lend = 1)
legend("topright",
       legend = c("barseq","LC"),
       col=c("coral","aquamarine3"),
       lty=1,
       lwd=10,
       cex=0.9
)
quartz.save("BarSeq_MetaNeighbour_largedatasetsLC_NA_res06.png", type = "png", dpi = 300, bg = "white")
dev.off()
