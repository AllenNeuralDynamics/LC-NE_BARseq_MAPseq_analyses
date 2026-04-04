################# BARseq Analysis Pipeline ##################
# Normalization, clustering, and spatial coherence analyses to isolate LC-NE neurons 
# have been performed on brain#3 and brain#4 datasets independently
# now bringing the data together to perform transcriptomic group identification, harmonize normalized matrices if necessary

# Load functions
source("~/capsule/code/1_BARseq_analyses_functions_brain3-4_combined.R")

# Set working directory
setwd(BARSEQ_OUTPUT_DIR)

##########################################################################################################
# Isolating LC-NE population fromFULL dataset where brain3 and brain4 cells are combined at post-QC step
##########################################################################################################
############################################################################################################################################################################################################
# load normalized gene expression data, concatenate the matrices - this is all cells which pass QC filter regardless of identity
brain3 <- readRDS("/results/BARseq_780345/combined_neurons_clust_CCFv2_uid_cpm_log.rds") 
dim (brain3)

brain4 <- readRDS("/results/BARseq_780346/combined_neurons_clust_CCFv2_uid_cpm_log.rds") 
dim(brain4)

# Examine why gene vectors are different length, only keep shared genes
genes_brain3 <- rownames(brain3)
genes_brain4 <- rownames(brain4)

unique_to_brain3 <- setdiff(genes_brain3, genes_brain4)
unique_to_brain4 <- setdiff(genes_brain4, genes_brain3)

# Define the shared gene set AND the order you want to enforce.
# Use brain3's order (stable, reproducible) for the common genes:
common_genes <- genes_brain3[genes_brain3 %in% genes_brain4]

cat("Genes unique to brain3:", length(unique_to_brain3), "\n")
cat("Genes unique to brain4:", length(unique_to_brain4), "\n")
cat("Common genes:", length(common_genes), "\n")

head(unique_to_brain3)
head(unique_to_brain4)

# Subset BOTH objects using the SAME ordered vector
brain3_subset <- brain3[common_genes, ]
brain4_subset <- brain4[common_genes, ]

# Hard checks: same genes, same order, no NAs
stopifnot(
  length(common_genes) > 0,
  !anyNA(common_genes),
  identical(rownames(brain3_subset), rownames(brain4_subset))
)

if (any(duplicated(genes_brain3)) || any(duplicated(genes_brain4))) {
  stop("Duplicate gene names detected in rownames(); fix before intersect/subsetting.")
}
# Concatenate along columns (cells)
combined_sce <- cbind(brain3_subset, brain4_subset)
# Verify dimensions
dim(combined_sce) 
# Add a batch column
colData(combined_sce)$batch <- c(rep("brain3", ncol(brain3_subset)), rep("brain4", ncol(brain4)))

####################################### perform initial clustering and check for batch effects ###############################################################################################
# Check if clustering analysis already exists given this is large data set and takes a while to run
# Paths: check pre-computed data asset first, fall back to writing new results
clustering_data_dir <- file.path(BARSEQ_CLUSTERING_DIR, "barseq_all_QCed_cells")
clustering_results_dir <- file.path(BARSEQ_OUTPUT_DIR, "analysis/barseq_all_QCed_cells")

umap_file    <- file.path(clustering_data_dir, "umap.csv")
cluster_file <- file.path(clustering_data_dir, "cluster.csv")
annot_file   <- file.path(clustering_data_dir, "cluster_annotation.csv")

if (dir.exists(clustering_data_dir) && file.exists(umap_file) && file.exists(cluster_file) && file.exists(annot_file)) {
  cat("All analysis results already exist in data asset. Skipping clustering.\n")
} else if (dir.exists(clustering_data_dir) && file.exists(umap_file)) {
  cat("UMAP exists but clustering incomplete. Loading UMAP and running PCA-based clustering...\n")
  umap_mat <- as.matrix(read_csv(umap_file))
  rownames(umap_mat) <- colnames(combined_sce)
  reducedDim(combined_sce, "UMAP") <- umap_mat
  v <- analyze_barseq_pca_only(combined_sce, "barseq_all_QCed_cells")
  # update paths to newly written results
  umap_file    <- file.path(clustering_results_dir, "umap.csv")
  cluster_file <- file.path(clustering_results_dir, "cluster.csv")
} else {
  cat("Running full clustering analysis for barseq_all_QCed_cells...\n")
  v <- analyze_barseq(combined_sce, "barseq_all_QCed_cells")
  # update paths to newly written results
  umap_file    <- file.path(clustering_results_dir, "umap.csv")
  cluster_file <- file.path(clustering_results_dir, "cluster.csv")
}

# Load umap and cluster for plotting — paths already set correctly above
umap_data <- read.csv(umap_file)
clusters  <- read.csv(cluster_file)
x <- umap_data[['UMAP1']]
y <- umap_data[['UMAP2']]

#visualize UMAP of samples 
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)
umap_data$row.names <- rownames(umap_data)
clusters$row.names <- rownames(clusters)
centroid_data <- umap_data %>%
  dplyr::left_join(clusters, by = "row.names") %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

plot <- ggplot(umap_data, aes(x=x, y=y, color=as.factor(clusters[["label"]]))) +
  geom_point(size=0.02) +
  geom_text(data = centroid_data, aes(x=x, y=y, label=as.factor(label)), colour="black", vjust=1.6, hjust=0.5, size=3.5) +
  xlab("UMAP1") +
  ylab("UMAP2") +
  scale_color_manual(values=color_palette, guide = guide_legend(override.aes = list(size=4))) +
  ggtitle("All neurons") +
  theme_minimal() +
  theme(panel.grid = element_blank()) #+ theme(legend.position = "none")
print(plot)
ggsave("all_QCed_cells_clust.pdf", plot = plot, device = "pdf", width = 10, height = 8)


# Plot parameters of interest onto a UMAP 
total_genes <- colSums(counts(combined_sce))
# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[['UMAP1']], UMAP2 = umap_data[['UMAP2']], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[['UMAP1']], UMAP2 = umap_data[['UMAP2']], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[['UMAP1']], UMAP2 = umap_data[['UMAP2']], batch = colData(combined_sce)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.01) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = as.factor(label)), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP Clusters", x = "UMAP1", y = "UMAP2", color = "Cluster")

p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.01) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.01, alpha = 0.1) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("all_QCed_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[['UMAP1']], UMAP2 = umap_data[['UMAP2']], cluster = factor(clusters[["label"]]), batch = colData(combined_sce)$batch)
# Subset for brain3 and brain4
plot_data_brain3 <- plot_data_cluster[plot_data_cluster$batch == "brain3", ]
plot_data_brain4 <- plot_data_cluster[plot_data_cluster$batch == "brain4", ]

# Overall plots
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = as.factor(label)), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(title = "All Clusters", x = "UMAP1", y = "UMAP2")

p2 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  # Smaller size
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "All Batch", x = "UMAP1", y = "UMAP2")

# Brain3 subset plots
p3 <- ggplot(plot_data_brain3, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain3 Clusters", x = "UMAP1", y = "UMAP2")

p4 <- ggplot(plot_data_brain3, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  
  scale_color_manual(values = c("brain3" = "blue")) + 
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain3 Batch", x = "UMAP1", y = "UMAP2")

# Brain4 subset plots
p5 <- ggplot(plot_data_brain4, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain4 Clusters", x = "UMAP1", y = "UMAP2")

p6 <- ggplot(plot_data_brain4, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  
  scale_color_manual(values = c("brain4" = "green")) +  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain4 Batch", x = "UMAP1", y = "UMAP2")

p_combined2 <- grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6), ncol = 2, nrow = 3)
print(p_combined2)
ggsave("all_QCed_cells_clusters_batch_subsets.pdf", plot = p_combined2, device = "pdf", width = 8, height = 16)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(combined_sce)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
print(p_combined3)
ggsave("all_QCed_cells_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

############################################################# check gene expression, add clusters to the sce object, save it ###############################################################
# Check select gene expression across clusters
gene_expression_summary <- plot_data_cluster %>%
  group_by(cluster) %>%
  summarise(
    Dbh_mean = mean(Dbh, na.rm = TRUE),
    Th_mean = mean(Th, na.rm = TRUE),
    Ddc_mean = mean(Ddc, na.rm = TRUE),
    Slc18a2_mean = mean(Slc18a2, na.rm = TRUE)
  )
print(gene_expression_summary)
write.csv(gene_expression_summary, file = "all_QCed_cells_gene_expression_summary_by_cluster.csv", row.names = FALSE)

# Calculate mean expression for each gene across clusters
cluster_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(combined_sce)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "all_QCed_cells_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()

# Calculate variance of gene expression across clusters
gene_variance <- apply(mean_expression_matrix, 2, var)
# Sort genes by variance
highly_variable_genes <- sort(gene_variance, decreasing = TRUE)
# Select top N highly variable genes
top_genes <- names(head(highly_variable_genes, 30))  # Adjust the number as needed
# Subset the mean expression matrix for these genes
top_genes_matrix <- mean_expression_matrix[, top_genes]
# Plot heatmap
pheatmap(top_genes_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Top Highly Variable Genes by Cluster")
dev.copy(pdf, "all_QCed_cells_gene_expr_HVG_heatmap.pdf", width = 14, height = 8)
dev.off()

# Check the distribution of cluster labels
table(clusters$label)
# Define a threshold for gene expression (e.g., logcounts > 0)
threshold <- 0
# Subset cells from cluster 11 
cluster_11_cells <- which(clusters$label == 11)
# Extract logcounts for the genes of interest
genes_of_interest <- c("Dbh", "Th", "Slc18a2", "Ddc")
gene_expression <- logcounts(combined_sce)[genes_of_interest, cluster_11_cells]
# Create a logical matrix: TRUE if expression > threshold, FALSE otherwise
expression_matrix <- gene_expression > threshold
# Count cells expressing each gene alone and in combination
# Convert to a data frame for easier manipulation
expression_df <- as.data.frame(t(expression_matrix))
colnames(expression_df) <- genes_of_interest
# Add a column for the combination of expressed genes
expression_df$combination <- apply(expression_df, 1, function(row) paste(names(row)[row], collapse = ", "))
# Count the number of cells for each combination
combination_counts <- expression_df %>%
  group_by(combination) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))
# View the results
print(combination_counts)
write.csv(combination_counts, file = "all_QCed_cells_gene_coexpression_summary_LC_cluster.csv", row.names = FALSE)

# Ensure that the sample names in clusters match the column names of barseq
if (!all(clusters[['sample']] == colnames(combined_sce))) {
  # Align clusters to match the column names of barseq
  clusters <- clusters[match(colnames(combined_sce), clusters[['sample']]), ]
  # Check if there are any unmatched samples
  if (any(is.na(clusters[['sample']]))) {
    stop("Error: Some samples in barseq do not have matching entries in clusters.")
  }
}
# Add cluster assignment data to the SingleCellExperiment object
colData(combined_sce)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(combined_sce, "combined_allQCed_neurons_clust_CCFv2_uid_cpm_log_clust.rds")

# Subset out LC cluster such that only LC cluster neurons are retained for normalization
dim(combined_sce)
X <- 11 # specify which cluster to subset
LC_barseq <- combined_sce[, colData(combined_sce)$louvain_cluster == X]
dim(LC_barseq)

# Save the LC cluster containing object
saveRDS(LC_barseq, file = "LCcluster_neurons_CCFv2_uid.rds")

# Normalize and save LC cluster object
LC <- load_barseq(filename="LCcluster_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(LC)
summary(Matrix::colSums(assay(LC, "cpm")))
anyDuplicated(colnames(LC))
colnames(LC)[duplicated(colnames(LC))][1:10]
# sanity: metadata rows correspond to columns
stopifnot(identical(rownames(colData(LC)), colnames(LC)))
stopifnot(!anyNA(colData(LC)$batch))

# keep original IDs
colData(LC)$cell_id_old <- colnames(LC)

# build new IDs using batch + uid (or batch + existing colname)
new_ids <- paste0(colData(LC)$batch, "|", colData(LC)$uid)   # uses your 'uid' column
new_ids <- make.unique(new_ids)

colnames(LC) <- new_ids
logcounts(LC) = log1p(assay(LC, "cpm"))/log(2)
saveRDS(LC,file.path(BARSEQ_OUTPUT_DIR, "LCcluster_neurons_CCFv2_uid_cpm_log.rds"))


############################################################# cluster cells from putative LC group after they have been re-normalized ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
#take cells belonging to LC cluster and cluster them again after they have been subset and normalized again in script #2
LC_barseq <- readRDS("LCcluster_neurons_CCFv2_uid_cpm_log.rds")
dim(LC_barseq)
colData(LC_barseq)
assayNames(LC_barseq)

v<-analyze_barseq(LC_barseq, "barseq_LCcluster_cells")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LC_barseq)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.01, alpha = 0.1) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCcluster_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LC_barseq)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
print(p_combined3)
ggsave("LCcluster_cells_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LC_barseq)) # should be true
colData(LC_barseq)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LC_barseq, "LCcluster_neurons_CCFv2_uid_cpm_log_clust.rds")

# Plot clustered cells to visualize their locations to ensure only proper LC-NE are retained for the subsequent steps
# Extract colData as a data.frame
LC_barseq_df <- as.data.frame(colData(LC_barseq))
# Ensure louvain_cluster is a factor
LC_barseq_df$louvain_cluster <- as.factor(LC_barseq_df$louvain_cluster)
# Ensure slice is ordered properly as a numeric factor
LC_barseq_df$slice <- factor(LC_barseq_df$slice, levels = sort(as.numeric(unique(LC_barseq_df$slice))))

# Create and save the combined plot
combined_plot <- create_plot(LC_barseq_df)
print(combined_plot)
ggsave("LC_cluster_cells_slices.pdf", plot = combined_plot, device = "pdf", width = 20, height = 12)

# Create and save the plot for brain3
brain3_df <- LC_barseq_df[LC_barseq_df$batch == "brain3", ]
if (nrow(brain3_df) > 0) {
  brain3_plot <- create_plot(brain3_df, " - Brain3")
  print(brain3_plot)
  ggsave("LC_cluster_cells_slices_brain3.pdf", plot = brain3_plot, device = "pdf", width = 20, height = 12)
}

# Create and save the plot for brain4
brain4_df <- LC_barseq_df[LC_barseq_df$batch == "brain4", ]
if (nrow(brain4_df) > 0) {
  brain4_plot <- create_plot(brain4_df, " - Brain4")
  print(brain4_plot)
  ggsave("LC_cluster_cells_slices_brain4.pdf", plot = brain4_plot, device = "pdf", width = 20, height = 12)
}

# Subset out LC-NE proper cells
table(clusters$label)
LCNE_barseq <- LC_barseq[, colData(LC_barseq)$louvain_cluster %in% c(1,2,6)]
dim(LCNE_barseq)
saveRDS(LCNE_barseq, "LCNE_cluster_neurons_CCFv2_uid.rds")

# Normalize and save LC-NE cluster object
LCNE_barseq <- load_barseq(filename="LCNE_cluster_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(LCNE_barseq)
logcounts(LCNE_barseq) = log1p(assay(LCNE_barseq, "cpm"))/log(2)
summary(Matrix::colSums(assay(LCNE_barseq, "logcounts")))
saveRDS(LCNE_barseq,file.path(BARSEQ_OUTPUT_DIR, "LCNE_cluster_neurons_CCFv2_uid_cpm_log.rds"))

############################################################# cluster LC-NE cells after they have been re-normalized ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
# Take LC-NE cells and cluster them again
LCNE_barseq <- readRDS("LCNE_cluster_neurons_CCFv2_uid_cpm_log.rds")
dim(LCNE_barseq)
colData(LCNE_barseq)
assayNames(LCNE_barseq)

v<-analyze_barseq(LCNE_barseq, "barseq_LC_NE_cells")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LCNE_barseq)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.05, alpha = 0.6) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNEcluster_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LCNE_barseq)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 4, nrow = 1)
print(p_combined3)
ggsave("LCNEcluster_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LCNE_barseq)) # should be true
colData(LCNE_barseq)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LCNE_barseq, "LCNE_neurons_CCFv2_uid_cpm_log_clust.rds")


#plot gene expression hetamap sorted by cluster
# Extract log-normalized counts
expr_matrix <- logcounts(LCNE_barseq)
# Extract cluster labels
cluster_labels <- colData(LCNE_barseq)$louvain_cluster
# Sort the expression matrix by cluster labels
sorted_indices <- order(cluster_labels)
expr_matrix_sorted <- expr_matrix[, sorted_indices]
# Update column names to reflect sorted cluster labels
colnames(expr_matrix_sorted) <- cluster_labels[sorted_indices]
# Add gene names as row names (if not already present)
rownames(expr_matrix_sorted) <- rownames(LCNE_barseq)
# Create a data frame for column annotations
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices]))  # Convert to factor
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))  # Unique row names
# Update column names of expr_matrix_sorted to match annotation_col row names
colnames(expr_matrix_sorted) <- rownames(annotation_col)
# Map the custom colors to the cluster levels
cluster_colors <- setNames(color_palette[1:length(levels(annotation_col$Cluster))], 
                           levels(annotation_col$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Plot heatmap using pheatmap with annotations
pheatmap(expr_matrix_sorted,
         cluster_rows = FALSE,  # Disable row clustering to preserve order
         cluster_cols = FALSE,  # Disable column clustering to preserve sorting
         scale = "none",        # No scaling to match the original behavior
         color = viridis::viridis(50),  # Use viridis color palette
         show_rownames = TRUE,  # Show gene names as row labels
         show_colnames = FALSE, # Hide dense column labels
         annotation_col = annotation_col,  # Add cluster annotations
         annotation_colors = annotation_colors,  # Use custom cluster-specific colors
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,      # Adjust font size for rows
         legend = TRUE)         # Ensure legend is displayed
dev.copy(pdf, "Log-count gene expression in LCNE clusters.pdf", width = 8, height = 10)
dev.off()

# Calculate mean expression for each gene across clusters
clusters_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(LCNE_barseq)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "LCNE_cluster_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()

# Plot clustered cells to visualize their locations to ensure only proper LC-NE are retained for the subsequent steps
# Extract colData as a data.frame
LCNE_barseq_df <- as.data.frame(colData(LCNE_barseq))
# Ensure louvain_cluster is a factor
LCNE_barseq_df$louvain_cluster <- as.factor(LCNE_barseq_df$louvain_cluster)
# Ensure slice is ordered properly as a numeric factor
LCNE_barseq_df$slice <- factor(LCNE_barseq_df$slice, levels = sort(as.numeric(unique(LCNE_barseq_df$slice))))

# Create and save the combined plot
combined_plot <- create_plot(LCNE_barseq_df)
print(combined_plot)
ggsave("LCNE_cluster_cells_slices.pdf", plot = combined_plot, device = "pdf", width = 20, height = 12)

# Create and save the plot for brain3
brain3_df <- LCNE_barseq_df[LCNE_barseq_df$batch == "brain3", ]
if (nrow(brain3_df) > 0) {
  brain3_plot <- create_plot(brain3_df, " - Brain3")
  print(brain3_plot)
  ggsave("LCNE_cluster_cells_slices_brain3.pdf", plot = brain3_plot, device = "pdf", width = 20, height = 12)
}

# Create and save the plot for brain4
brain4_df <- LCNE_barseq_df[LCNE_barseq_df$batch == "brain4", ]
if (nrow(brain4_df) > 0) {
  brain4_plot <- create_plot(brain4_df, " - Brain4")
  print(brain4_plot)
  ggsave("LCNE_cluster_cells_slices_brain4.pdf", plot = brain4_plot, device = "pdf", width = 20, height = 12)
}

############################################################## Check spatial density and coherence of LC-NE clusters filtered cells ########################################################################
# Clean up between processing steps
clear_objects_except_functions()
# Load data and stored UMAP info
LCNE_barseq <- readRDS("LCNE_neurons_CCFv2_uid_cpm_log_clust.rds")

umap_path <- "/results/BARseq_780345-780346_combined/analysis/barseq_LC_NE_cells/umap.csv"
umap_df <- read.csv(umap_path, header = TRUE)

#Compute kNN-based cluster purity: proportion of same-cluster neighbors.
# Distance weighted spatial coherence score calculation
# Run the function
LCNE_barseq <- calculate_spatial_coherence_3D(LCNE_barseq, cluster_col = "louvain_cluster", k = 10, slice_thickness = 20)
colData(LCNE_barseq)

# Extract metadata
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(louvain_cluster = as.factor(louvain_cluster))

# Violin plot for 3D coherence
n_clusters <- length(unique(LCNE_barseq$louvain_cluster))
color_palette <- get_cluster_colors(n_clusters)
# Violin plot
p <- ggplot(meta_df, aes(x = louvain_cluster, y = spatial_coherence_3D_weighted, fill = louvain_cluster)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.4) +
  scale_fill_manual(values = color_palette) +
  labs(title = "3D Spatial Coherence by Louvain Cluster",
       x = "Louvain Cluster", y = "3D Coherence Score (Weighted)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
ggsave("spatial_coherence_violin.pdf", plot=p, width = 8, height = 8, units = "in")

# Examine distribution and possible cut off points
summary(meta_df$spatial_coherence_3D_weighted)
p <- ggplot(meta_df, aes(x = spatial_coherence_3D_weighted)) +
  geom_histogram(bins = 100, fill = "lightgray", color = "black") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "Histogram proposed Coherence Cutoff",
       x = "Spatial Coherence (Weighted)", y = "Cell count") +
  theme_minimal()
sum(meta_df$spatial_coherence_3D_weighted < 0.05)
print(p)
ggsave("spatial_coherence_hist.pdf", plot=p, width = 8, height = 8, units = "in")


# Extract metadata for UMAP space plotting
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(
    UMAP1 = umap_df$UMAP1,
    UMAP2 = umap_df$UMAP2,
    louvain_cluster = as.factor(louvain_cluster),
    low_coherence_3D_flag = spatial_coherence_3D_weighted < 0.05
  )
# Plot 3D coherence
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = spatial_coherence_3D_weighted)) +
  geom_point(size = 0.5) +
  scale_color_gradient(low = "lightgray", high = "royal blue") +
  theme_minimal() +
  labs(title = "3D Spatial Coherence per Cell", color = "Coherence Score (3D)")
print(p)
ggsave("spatial_coherence_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

# Highlight low 3D-coherence cells
p <-ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = low_coherence_3D_flag)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  theme_minimal() +
  labs(title = "Low Spatial Coherence Cells", color = "Low Coherence")
print(p)
ggsave("spatial_coherence_UMAP_cutoff.pdf", plot=p, width = 8, height = 8, units = "in")


# Calculate average kNN distance in 3D space, and inverts this distance → spatial density score (larger = denser).
LCNE_barseq <- calculate_spatial_density(LCNE_barseq, k = 10, slice_thickness = 20)
colData(LCNE_barseq)

# Extract metadata
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(louvain_cluster = as.factor(louvain_cluster))
# Violin plot for 3D density
p <- ggplot(meta_df, aes(x = louvain_cluster, y = spatial_density_3D, fill = louvain_cluster)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.4) +
  scale_fill_manual(values = color_palette) +
  labs(title = "3D Spatial Density by Louvain Cluster",
       x = "Louvain Cluster", y = "3D Spatial Density Score") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
ggsave("spatial_density_violin.pdf", plot=p, width = 8, height = 8, units = "in")

# Examine distribution and possible cut off points
summary(meta_df$spatial_density_3D)
cutoff <- quantile(meta_df$spatial_density_3D, probs = 0.01)  # bottom 1%
p <- ggplot(meta_df, aes(x = spatial_density_3D)) +
  geom_histogram(bins = 100, fill = "lightgray", color = "black") +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed") +
  labs(title = "Histogram with Bottom 1% Density Cutoff",
       x = "3D Spatial Density Score", y = "Cell count") +
  theme_minimal()
print(p)
sum(meta_df$spatial_density_3D < cutoff)
ggsave("spatial_density_hist_5perc_cut.pdf", plot=p, width = 8, height = 8, units = "in")

# Extract metadata for UMAP space plotting
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(
    UMAP1 = umap_df$UMAP1,
    UMAP2 = umap_df$UMAP2,
    louvain_cluster = as.factor(louvain_cluster),
    low_density_flag = spatial_density_3D < quantile(spatial_density_3D, 0.01) 
  )
# Density color map
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = spatial_density_3D)) +
  geom_point(size = 0.5) +
  scale_color_gradient(low = "gray", high = "royal blue") +
  theme_minimal() +
  labs(title = "Spatial Density Score (3D)", color = "Density Score")
print(p)
ggsave("spatial_density_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

# Highlight low-density cells
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = low_density_flag)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "orange")) +
  theme_minimal() +
  labs(title = "Low Spatial Density Cells", color = "Low Density")
print(p)
ggsave("spatial_density_UMAP_cut.pdf", plot=p, width = 8, height = 8, units = "in")

# Combine density and coherence measurements for exclusion criterion of suspect cells
# spatial_coherence_3D: ranges from 0 (neighbors from other clusters) to 1 (all same cluster). Higher = more spatially coherent.
# spatial_density_3D: higher = denser local neighborhood (from 1 / mean kNN distance).
# Set absolute coherence cutoff
coherence_cutoff <- 0.05  # Or any other fixed value you prefer
# Set quantile-based density cutoff (bottom 5%)
density_cutoff <- quantile(meta_df$spatial_density_3D, 0.01, na.rm = TRUE)
# Apply filtering
meta_df <- meta_df %>%
  mutate(
    high_coherence = spatial_coherence_3D_weighted >= coherence_cutoff,
    high_density   = spatial_density_3D >= density_cutoff,
    keep_cell = high_coherence & high_density
  )
table(meta_df$keep_cell)
# Visualize in UMAP space the cells to be excluded 
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = keep_cell)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "orange", "TRUE" = "darkgreen")) +
  theme_minimal() +
  labs(title = "Filtered Cells (High Coherence + Density)", color = "Retained")
print(p)
ggsave("combined_excluded_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

table(meta_df$louvain_cluster)
table(meta_df$keep_cell, meta_df$louvain_cluster, meta_df$batch)

# Filter out cells with low density and coherence scores
retained_cells <- LCNE_barseq[, meta_df$keep_cell]
dim(LCNE_barseq)
dim(retained_cells)
saveRDS(retained_cells, "LCNE_clusters_filtered_coherence_filtered.rds")

# Normalize and save LC-NE clusters filtered coherence filtered sce object
LCNE_barseq_clusters_filtered_coherence_filtered <- load_barseq(filename="LCNE_clusters_filtered_coherence_filtered.rds", from_output = TRUE)
colData(LCNE_barseq_clusters_filtered_coherence_filtered)
logcounts(LCNE_barseq_clusters_filtered_coherence_filtered) = log1p(assay(LCNE_barseq_clusters_filtered_coherence_filtered, "cpm"))/log(2)
saveRDS(LCNE_barseq_clusters_filtered_coherence_filtered,file.path(BARSEQ_OUTPUT_DIR, "LCNE_clusters_filtered_coherence_filtered_cpm_log.rds"))

############################################################# re-cluster LC-NE cluster cleaned up cells ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
# Take cleaned up LC-NE cells and cluster them again
LCNE <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log.rds")
dim(LCNE)
colData(LCNE)
assayNames(LCNE)

# Drop columns from previous processing
cols_to_drop <- c("louvain_cluster", "spatial_coherence_3D_weighted", "spatial_density_3D", "cell_id_old")
colData(LCNE) <- colData(LCNE)[, !(colnames(colData(LCNE)) %in% cols_to_drop)]

# Run clustering on the data and save results to a named folder
v<-analyze_barseq(LCNE, "barseq_LC_NE_clusters_filtered_coherence_filtered")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LCNE)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.05, alpha = 0.6) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNE_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LCNE)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 4, nrow = 1)
print(p_combined3)
ggsave("LCNE_cells_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LCNE)) # should be true
colData(LCNE)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LCNE, "LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")

#plot gene expression hetamap sorted by cluster
# Extract log-normalized counts
expr_matrix <- logcounts(LCNE)
# Extract cluster labels
cluster_labels <- colData(LCNE)$louvain_cluster
# Sort the expression matrix by cluster labels
sorted_indices <- order(cluster_labels)
expr_matrix_sorted <- expr_matrix[, sorted_indices]
# Update column names to reflect sorted cluster labels
colnames(expr_matrix_sorted) <- cluster_labels[sorted_indices]
# Add gene names as row names (if not already present)
rownames(expr_matrix_sorted) <- rownames(LCNE)
# Create a data frame for column annotations
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices]))  # Convert to factor
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))  # Unique row names
# Update column names of expr_matrix_sorted to match annotation_col row names
colnames(expr_matrix_sorted) <- rownames(annotation_col)
# Map the custom colors to the cluster levels
cluster_colors <- setNames(color_palette[1:length(levels(annotation_col$Cluster))], 
                           levels(annotation_col$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Plot heatmap using pheatmap with annotations
pheatmap(expr_matrix_sorted,
         cluster_rows = FALSE,  # Disable row clustering to preserve order
         cluster_cols = FALSE,  # Disable column clustering to preserve sorting
         scale = "none",        # No scaling to match the original behavior
         color = viridis::viridis(50),  # Use viridis color palette
         show_rownames = TRUE,  # Show gene names as row labels
         show_colnames = FALSE, # Hide dense column labels
         annotation_col = annotation_col,  # Add cluster annotations
         annotation_colors = annotation_colors,  # Use custom cluster-specific colors
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,      # Adjust font size for rows
         legend = TRUE)         # Ensure legend is displayed
dev.copy(pdf, "Log-count gene expression in LCNE clusters filtered coherence filtered.pdf", width = 10, height = 12)
dev.off()

# Sort gene expression within a cluster for heatmap plotting
# Calculate total expression per cell 
total_expr <- colSums(expr_matrix)
sorted_indices_within <- c()
for (clust in levels(cluster_labels)) {
  cells_in_cluster <- which(cluster_labels == clust)
  
  # Sort cells inside cluster by total expression (descending)
  order_in_cluster <- cells_in_cluster[order(total_expr[cells_in_cluster], decreasing = TRUE)]
  
  sorted_indices_within <- c(sorted_indices_within, order_in_cluster)
}
# Subset and reorder expression matrix and annotation
expr_matrix_sorted <- expr_matrix[, sorted_indices_within]
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices_within]))
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))
colnames(expr_matrix_sorted) <- rownames(annotation_col)
pheatmap(expr_matrix_sorted,
         cluster_rows = TRUE,   # Let genes be clustered to reveal patterns
         cluster_cols = FALSE,  # Keep your cell ordering intact
         scale = "none",
         color = viridis::viridis(50),
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,
         legend = TRUE)
dev.copy(pdf, "Log-count gene expression in LCNE  clusters filtered coherence filtered sorted.pdf", width = 10, height = 12)
dev.off()

# Calculate mean expression for each gene across clusters
clusters_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(LCNE)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "LCNE_cells_clusters_filtered_coherence_filtered_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()



################################################################################################################################################################################################################
######################################################## compare the LC-NE final overlap between combined and separate processing ##############################################################################
################################################################################################################################################################################################################
fromFULL <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(fromFULL)
summary(Matrix::colSums(assay(fromFULL, "cpm")))
summary(Matrix::colSums(assay(fromFULL, "logcounts")))

# Load separately processed datasets and concatenate them
brain3 <- readRDS("/results/BARseq_780345/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(brain3)
summary(Matrix::colSums(assay(brain3, "cpm")))
summary(Matrix::colSums(assay(brain3, "logcounts")))

brain4 <- readRDS("/results/BARseq_780346/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(brain4)
summary(Matrix::colSums(assay(brain4, "cpm")))
summary(Matrix::colSums(assay(brain4, "logcounts")))

# Examine why gene vectors are different length, only keep shared genes
genes_brain3 <- rownames(brain3)
genes_brain4 <- rownames(brain4)
unique_to_brain3 <- setdiff(genes_brain3, genes_brain4)
unique_to_brain4 <- setdiff(genes_brain4, genes_brain3)

# Define the shared gene set AND the order to enforce.
# Use brain3's order (stable, reproducible) for the common genes:
common_genes <- genes_brain3[genes_brain3 %in% genes_brain4]

cat("Genes unique to brain3:", length(unique_to_brain3), "\n")
cat("Genes unique to brain4:", length(unique_to_brain4), "\n")
cat("Common genes:", length(common_genes), "\n")

head(unique_to_brain3)
head(unique_to_brain4)

# Subset BOTH objects using the SAME ordered vector
brain3_subset <- brain3[common_genes, ]
brain4_subset <- brain4[common_genes, ]

# Hard checks: same genes, same order, no NAs
stopifnot(
  length(common_genes) > 0,
  !anyNA(common_genes),
  identical(rownames(brain3_subset), rownames(brain4_subset))
)

# Concatenate along columns (cells)
combined_sce <- cbind(brain3_subset, brain4_subset)
# Verify dimensions
dim(combined_sce)
# Add a batch column
colData(combined_sce)$batch <- c(rep("brain3", ncol(brain3_subset)), rep("brain4", ncol(brain4)))

dim(fromFULL)
colnames(fromFULL)
colData(fromFULL)

fromLCNE <- combined_sce
rm(combined_sce)
dim(fromLCNE)
colnames(fromLCNE)
colData (fromLCNE)
anyDuplicated(colnames(fromLCNE))
colnames(fromLCNE)[duplicated(colnames(fromLCNE))][1:10]
# sanity: metadata rows correspond to columns
stopifnot(identical(rownames(colData(fromLCNE)), colnames(fromLCNE)))
stopifnot(!anyNA(colData(fromLCNE)$batch))
# build new IDs using batch + uid (or batch + existing colname)
new_ids <- paste0(colData(fromLCNE)$batch, "|", colData(fromLCNE)$uid)   # uses your 'uid' column
new_ids <- make.unique(new_ids)
colnames(fromLCNE) <- new_ids
saveRDS(fromLCNE,file.path(BARSEQ_OUTPUT_DIR, "fromLCNE_combined_neurons_CCFv2_uid.rds"))

# Normalize and save fromLCNE combined object
combined_LCNE_barseq <- load_barseq(filename="fromLCNE_combined_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(combined_LCNE_barseq)
logcounts(combined_LCNE_barseq) = log1p(assay(combined_LCNE_barseq, "cpm"))/log(2)
summary(Matrix::colSums(assay(combined_LCNE_barseq, "cpm")))
summary(Matrix::colSums(assay(combined_LCNE_barseq, "logcounts")))
saveRDS(combined_LCNE_barseq,file.path(BARSEQ_OUTPUT_DIR, "fromLCNE_combined_neurons_CCFv2_uid_cpm_log.rds"))

####################################### perform initial clustering of the combined dataset and check for batch effects ###############################################################################################
v<-analyze_barseq(combined_LCNE_barseq, "fromLCNE_combined")
new_barseq <- v[[1]]
clusters <- v[[2]]
#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)
# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))
# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(combined_LCNE_barseq)$batch)
# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)
# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.01) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster),  # Fixed: Use 'cluster' instead of 'label'
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP Clusters", x = "UMAP1", y = "UMAP2", color = "Cluster")

p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.01) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Total Count")

p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.01, alpha = 0.5) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("fromLCNE_combined_QCed_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]), batch = colData(combined_LCNE_barseq)$batch)
# Subset for brain3 and brain4
plot_data_brain3 <- plot_data_cluster[plot_data_cluster$batch == "brain3", ]
plot_data_brain4 <- plot_data_cluster[plot_data_cluster$batch == "brain4", ]

# Overall plots
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = as.factor(cluster)), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(title = "All Clusters", x = "UMAP1", y = "UMAP2")

p2 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  # Smaller size
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "All Batch", x = "UMAP1", y = "UMAP2")

# Brain3 subset plots
p3 <- ggplot(plot_data_brain3, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain3 Clusters", x = "UMAP1", y = "UMAP2")

p4 <- ggplot(plot_data_brain3, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  
  scale_color_manual(values = c("brain3" = "blue")) + 
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain3 Batch", x = "UMAP1", y = "UMAP2")

# Brain4 subset plots
p5 <- ggplot(plot_data_brain4, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.005) +  
  scale_color_manual(values = color_palette) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain4 Clusters", x = "UMAP1", y = "UMAP2")

p6 <- ggplot(plot_data_brain4, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.005, alpha = 0.5) +  
  scale_color_manual(values = c("brain4" = "green")) +  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Brain4 Batch", x = "UMAP1", y = "UMAP2")

p_combined2 <- grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6), ncol = 2, nrow = 3)
print(p_combined2)
ggsave("fromLCNE_combined_QCed_cells_clusters_batch_subsets.pdf", plot = p_combined2, device = "pdf", width = 8, height = 16)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(combined_LCNE_barseq)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
print(p_combined3)
ggsave("fromLCNE_combined_QCed_cells_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

############################################################# check gene expression, add clusters to the sce object, save it ###############################################################
# Check select gene expression across clusters
gene_expression_summary <- plot_data_cluster %>%
  group_by(cluster) %>%
  summarise(
    Dbh_mean = mean(Dbh, na.rm = TRUE),
    Th_mean = mean(Th, na.rm = TRUE),
    Ddc_mean = mean(Ddc, na.rm = TRUE),
    Slc18a2_mean = mean(Slc18a2, na.rm = TRUE)
  )
print(gene_expression_summary)
write.csv(gene_expression_summary, file = "fromLCNE_combined_gene_expression_summary_by_cluster.csv", row.names = FALSE)

# Calculate mean expression for each gene across clusters
cluster_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(combined_LCNE_barseq)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "fromLCNE_combined_gene_QCed_cells_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()

# Calculate variance of gene expression across clusters
gene_variance <- apply(mean_expression_matrix, 2, var)
# Sort genes by variance
highly_variable_genes <- sort(gene_variance, decreasing = TRUE)
# Select top N highly variable genes
top_genes <- names(head(highly_variable_genes, 30))  # Adjust the number as needed
# Subset the mean expression matrix for these genes
top_genes_matrix <- mean_expression_matrix[, top_genes]
# Plot heatmap
pheatmap(top_genes_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Top Highly Variable Genes by Cluster")
dev.copy(pdf, "fromLCNE_combined_gene_QCed_cells_gene_expr_HVG_heatmap.pdf", width = 14, height = 8)
dev.off()

# Check the distribution of cluster labels
table(clusters$label)
# Define a threshold for gene expression (e.g., logcounts > 0)
threshold <- 0
# Subset cells from clusters with good Dbh expression
cluster_LCNE_cells <- which(clusters$label %in% c(1, 2, 3, 4, 6))
# Extract logcounts for the genes of interest
genes_of_interest <- c("Dbh", "Th", "Slc18a2", "Ddc")
gene_expression <- logcounts(combined_LCNE_barseq)[genes_of_interest, cluster_LCNE_cells]
# Create a logical matrix: TRUE if expression > threshold, FALSE otherwise
expression_matrix <- gene_expression > threshold
# Count cells expressing each gene alone and in combination
# Convert to a data frame for easier manipulation
expression_df <- as.data.frame(t(expression_matrix))
colnames(expression_df) <- genes_of_interest
# Add a column for the combination of expressed genes
expression_df$combination <- apply(expression_df, 1, function(row) paste(names(row)[row], collapse = ", "))
# Count the number of cells for each combination
combination_counts <- expression_df %>%
  group_by(combination) %>%
  summarise(cell_count = n()) %>%
  arrange(desc(cell_count))
# View the results
print(combination_counts)
write.csv(combination_counts, file = "fromLCNE_combined_gene_coexpression_summary_LC_cluster.csv", row.names = FALSE)

# Ensure that the sample names in clusters match the column names of barseq
if (!all(clusters[['sample']] == colnames(combined_LCNE_barseq))) {
  # Align clusters to match the column names of barseq
  clusters <- clusters[match(colnames(combined_LCNE_barseq), clusters[['sample']]), ]
  # Check if there are any unmatched samples
  if (any(is.na(clusters[['sample']]))) {
    stop("Error: Some samples in barseq do not have matching entries in clusters.")
  }
}
# Add cluster assignment data to the SingleCellExperiment object
colData(combined_LCNE_barseq)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(combined_LCNE_barseq, "fromLCNE_combined_neurons_clust_CCFv2_uid_cpm_log_clust.rds")

# Subset out LC cluster such that only LC cluster neurons are retained for normalization
dim(combined_LCNE_barseq)
X <- c(1,2,3,4,6) # specify which clusters to subset
LC_barseq <- combined_LCNE_barseq[, colData(combined_LCNE_barseq)$louvain_cluster %in% X]
dim(LC_barseq)

# Save the LC cluster containing object
saveRDS(LC_barseq, file = "fromLCNE_combined_LCcluster_neurons_CCFv2_uid.rds")

# Normalize and save LC cluster object
LC <- load_barseq(filename="fromLCNE_combined_LCcluster_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(LC)
summary(Matrix::colSums(assay(LC, "cpm")))
logcounts(LC) = log1p(assay(LC, "cpm"))/log(2)
summary(Matrix::colSums(assay(LC, "logcounts")))
saveRDS(LC,file.path(BARSEQ_OUTPUT_DIR, "fromLCNE_combined_LCcluster_neurons_CCFv2_uid_cpm_log.rds"))

############################################################# cluster cells from putative LC group after they have been re-normalized ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
#take cells belonging to LC cluster and cluster them again after they have been subset and normalized again in script #2
LC_barseq <- readRDS("fromLCNE_combined_LCcluster_neurons_CCFv2_uid_cpm_log.rds")
dim(LC_barseq)
colData(LC_barseq)
assayNames(LC_barseq)

v<-analyze_barseq(LC_barseq, "fromLCNE_combined_barseq_LCcluster_cells")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LC_barseq)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Total Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.01, alpha = 0.1) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("fromLCNE_combined_LCcluster_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LC_barseq)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
print(p_combined3)
ggsave("fromLCNE_combined_LCcluster_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 14)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LC_barseq)) # should be true
colData(LC_barseq)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LC_barseq, "fromLCNE_combined_LCcluster_neurons_CCFv2_uid_cpm_log_clust.rds")

# Plot clustered cells to visualize their locations to ensure only proper LC-NE are retained for the subsequent steps
# Extract colData as a data.frame
LC_barseq_df <- as.data.frame(colData(LC_barseq))
# Ensure louvain_cluster is a factor
LC_barseq_df$louvain_cluster <- as.factor(LC_barseq_df$louvain_cluster)
# Ensure slice is ordered properly as a numeric factor
LC_barseq_df$slice <- factor(LC_barseq_df$slice, levels = sort(as.numeric(unique(LC_barseq_df$slice))))
write.csv(LC_barseq_df, file = "fromLCNE_combined_LCcluster_neurons_CCFv2_uid_cpm_log_clust_colData.csv", row.names = FALSE)

# Create and save the combined plot
combined_plot <- create_plot(LC_barseq_df)
print(combined_plot)
ggsave("fromLCNE_combined_LC_cluster_cells_slices.pdf", plot = combined_plot, device = "pdf", width = 20, height = 12)

# Create and save the plot for brain3
brain3_df <- LC_barseq_df[LC_barseq_df$batch == "brain3", ]
if (nrow(brain3_df) > 0) {
  brain3_plot <- create_plot(brain3_df, " - Brain3")
  print(brain3_plot)
  ggsave("fromLCNE_combined_LC_cluster_cells_slices_brain3.pdf", plot = brain3_plot, device = "pdf", width = 20, height = 12)
}

# Create and save the plot for brain4
brain4_df <- LC_barseq_df[LC_barseq_df$batch == "brain4", ]
if (nrow(brain4_df) > 0) {
  brain4_plot <- create_plot(brain4_df, " - Brain4")
  print(brain4_plot)
  ggsave("fromLCNE_combined_LC_cluster_cells_slices_brain4.pdf", plot = brain4_plot, device = "pdf", width = 20, height = 12)
}












# Subset out LC-NE proper cells
table(clusters$label)
LCNE_barseq <- LC_barseq[, colData(LC_barseq)$louvain_cluster %in% c(1,2)]
dim(LCNE_barseq)
saveRDS(LCNE_barseq, "LCNE_cluster_neurons_CCFv2_uid.rds")

# Normalize and save LC-NE cluster object
LCNE_barseq <- load_barseq(filename="LCNE_cluster_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(LCNE_barseq)
logcounts(LCNE_barseq) = log1p(assay(LCNE_barseq, "cpm"))/log(2)
summary(Matrix::colSums(assay(LCNE_barseq, "logcounts")))
saveRDS(LCNE_barseq,file.path(BARSEQ_OUTPUT_DIR, "LCNE_cluster_neurons_CCFv2_uid_cpm_log.rds"))

############################################################# cluster LC-NE cells after they have been re-normalized ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
# Take LC-NE cells and cluster them again
LCNE_barseq <- readRDS("LCNE_cluster_neurons_CCFv2_uid_cpm_log.rds")
dim(LCNE_barseq)
colData(LCNE_barseq)
assayNames(LCNE_barseq)

v<-analyze_barseq(LCNE_barseq, "barseq_LC_NE_cells")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LCNE_barseq)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.05, alpha = 0.6) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNEcluster_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LCNE_barseq)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 4, nrow = 1)
print(p_combined3)
ggsave("LCNEcluster_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LCNE_barseq)) # should be true
colData(LCNE_barseq)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LCNE_barseq, "LCNE_neurons_CCFv2_uid_cpm_log_clust.rds")

#plot gene expression hetamap sorted by cluster
# Extract log-normalized counts
expr_matrix <- logcounts(LCNE_barseq)
# Extract cluster labels
cluster_labels <- colData(LCNE_barseq)$louvain_cluster
# Sort the expression matrix by cluster labels
sorted_indices <- order(cluster_labels)
expr_matrix_sorted <- expr_matrix[, sorted_indices]
# Update column names to reflect sorted cluster labels
colnames(expr_matrix_sorted) <- cluster_labels[sorted_indices]
# Add gene names as row names (if not already present)
rownames(expr_matrix_sorted) <- rownames(LCNE_barseq)
# Create a data frame for column annotations
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices]))  # Convert to factor
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))  # Unique row names
# Update column names of expr_matrix_sorted to match annotation_col row names
colnames(expr_matrix_sorted) <- rownames(annotation_col)
# Map the custom colors to the cluster levels
cluster_colors <- setNames(color_palette[1:length(levels(annotation_col$Cluster))], 
                           levels(annotation_col$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Plot heatmap using pheatmap with annotations
pheatmap(expr_matrix_sorted,
         cluster_rows = FALSE,  # Disable row clustering to preserve order
         cluster_cols = FALSE,  # Disable column clustering to preserve sorting
         scale = "none",        # No scaling to match the original behavior
         color = viridis::viridis(50),  # Use viridis color palette
         show_rownames = TRUE,  # Show gene names as row labels
         show_colnames = FALSE, # Hide dense column labels
         annotation_col = annotation_col,  # Add cluster annotations
         annotation_colors = annotation_colors,  # Use custom cluster-specific colors
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,      # Adjust font size for rows
         legend = TRUE)         # Ensure legend is displayed
dev.copy(pdf, "Log-count gene expression in LCNE clusters.pdf", width = 8, height = 10)
dev.off()

# Calculate mean expression for each gene across clusters
clusters_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(LCNE_barseq)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "LCNE_cluster_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()

# Plot clustered cells to visualize their locations to ensure only proper LC-NE are retained for the subsequent steps
# Extract colData as a data.frame
LCNE_barseq_df <- as.data.frame(colData(LCNE_barseq))
# Ensure louvain_cluster is a factor
LCNE_barseq_df$louvain_cluster <- as.factor(LCNE_barseq_df$louvain_cluster)
# Ensure slice is ordered properly as a numeric factor
LCNE_barseq_df$slice <- factor(LCNE_barseq_df$slice, levels = sort(as.numeric(unique(LCNE_barseq_df$slice))))

# Create and save the combined plot
combined_plot <- create_plot(LCNE_barseq_df)
print(combined_plot)
ggsave("LCNE_cluster_cells_slices.pdf", plot = combined_plot, device = "pdf", width = 20, height = 12)

# Create and save the plot for brain3
brain3_df <- LCNE_barseq_df[LCNE_barseq_df$batch == "brain3", ]
if (nrow(brain3_df) > 0) {
  brain3_plot <- create_plot(brain3_df, " - Brain3")
  print(brain3_plot)
  ggsave("LCNE_cluster_cells_slices_brain3.pdf", plot = brain3_plot, device = "pdf", width = 20, height = 12)
}

# Create and save the plot for brain4
brain4_df <- LCNE_barseq_df[LCNE_barseq_df$batch == "brain4", ]
if (nrow(brain4_df) > 0) {
  brain4_plot <- create_plot(brain4_df, " - Brain4")
  print(brain4_plot)
  ggsave("LCNE_cluster_cells_slices_brain4.pdf", plot = brain4_plot, device = "pdf", width = 20, height = 12)
}

############################################################## Check spatial density and coherence of LC-NE clusters filtered cells ########################################################################
# Clean up between processing steps
clear_objects_except_functions()
# Load data and stored UMAP info
LCNE_barseq <- readRDS("LCNE_neurons_CCFv2_uid_cpm_log_clust.rds")

umap_path <- "/results/BARseq_780345-780346_combined/analysis/barseq_LC_NE_cells/umap.csv"
umap_df <- read.csv(umap_path, header = TRUE)

#Compute kNN-based cluster purity: proportion of same-cluster neighbors.
# Distance weighted spatial coherence score calculation
# Run the function
LCNE_barseq <- calculate_spatial_coherence_3D(LCNE_barseq, cluster_col = "louvain_cluster", k = 10, slice_thickness = 20)
colData(LCNE_barseq)

# Extract metadata
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(louvain_cluster = as.factor(louvain_cluster))

# Violin plot for 3D coherence
n_clusters <- length(unique(LCNE_barseq$louvain_cluster))
color_palette <- get_cluster_colors(n_clusters)
# Violin plot
p <- ggplot(meta_df, aes(x = louvain_cluster, y = spatial_coherence_3D_weighted, fill = louvain_cluster)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.4) +
  scale_fill_manual(values = color_palette) +
  labs(title = "3D Spatial Coherence by Louvain Cluster",
       x = "Louvain Cluster", y = "3D Coherence Score (Weighted)") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
ggsave("spatial_coherence_violin.pdf", plot=p, width = 8, height = 8, units = "in")

# Examine distribution and possible cut off points
summary(meta_df$spatial_coherence_3D_weighted)
p <- ggplot(meta_df, aes(x = spatial_coherence_3D_weighted)) +
  geom_histogram(bins = 100, fill = "lightgray", color = "black") +
  geom_vline(xintercept = 0.05, color = "red", linetype = "dashed") +
  labs(title = "Histogram proposed Coherence Cutoff",
       x = "Spatial Coherence (Weighted)", y = "Cell count") +
  theme_minimal()
sum(meta_df$spatial_coherence_3D_weighted < 0.05)
print(p)
ggsave("spatial_coherence_hist.pdf", plot=p, width = 8, height = 8, units = "in")


# Extract metadata for UMAP space plotting
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(
    UMAP1 = umap_df$UMAP1,
    UMAP2 = umap_df$UMAP2,
    louvain_cluster = as.factor(louvain_cluster),
    low_coherence_3D_flag = spatial_coherence_3D_weighted < 0.05
  )
# Plot 3D coherence
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = spatial_coherence_3D_weighted)) +
  geom_point(size = 0.5) +
  scale_color_gradient(low = "lightgray", high = "royal blue") +
  theme_minimal() +
  labs(title = "3D Spatial Coherence per Cell", color = "Coherence Score (3D)")
print(p)
ggsave("spatial_coherence_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

# Highlight low 3D-coherence cells
p <-ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = low_coherence_3D_flag)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
  theme_minimal() +
  labs(title = "Low Spatial Coherence Cells", color = "Low Coherence")
print(p)
ggsave("spatial_coherence_UMAP_cutoff.pdf", plot=p, width = 8, height = 8, units = "in")


# Calculate average kNN distance in 3D space, and inverts this distance → spatial density score (larger = denser).
LCNE_barseq <- calculate_spatial_density(LCNE_barseq, k = 10, slice_thickness = 20)
colData(LCNE_barseq)

# Extract metadata
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(louvain_cluster = as.factor(louvain_cluster))
# Violin plot for 3D density
p <- ggplot(meta_df, aes(x = louvain_cluster, y = spatial_density_3D, fill = louvain_cluster)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.4) +
  scale_fill_manual(values = color_palette) +
  labs(title = "3D Spatial Density by Louvain Cluster",
       x = "Louvain Cluster", y = "3D Spatial Density Score") +
  theme_minimal() +
  theme(legend.position = "none")
print(p)
ggsave("spatial_density_violin.pdf", plot=p, width = 8, height = 8, units = "in")

# Examine distribution and possible cut off points
summary(meta_df$spatial_density_3D)
cutoff <- quantile(meta_df$spatial_density_3D, probs = 0.05)  # bottom 5%
p <- ggplot(meta_df, aes(x = spatial_density_3D)) +
  geom_histogram(bins = 100, fill = "lightgray", color = "black") +
  geom_vline(xintercept = cutoff, color = "red", linetype = "dashed") +
  labs(title = "Histogram with Bottom 5% Density Cutoff",
       x = "3D Spatial Density Score", y = "Cell count") +
  theme_minimal()
print(p)
sum(meta_df$spatial_density_3D < cutoff)
ggsave("spatial_density_hist_5perc_cut.pdf", plot=p, width = 8, height = 8, units = "in")

# Extract metadata for UMAP space plotting
meta_df <- as.data.frame(colData(LCNE_barseq)) %>%
  mutate(
    UMAP1 = umap_df$UMAP1,
    UMAP2 = umap_df$UMAP2,
    louvain_cluster = as.factor(louvain_cluster),
    low_density_flag = spatial_density_3D < quantile(spatial_density_3D, 0.05) 
  )
# Density color map
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = spatial_density_3D)) +
  geom_point(size = 0.5) +
  scale_color_gradient(low = "gray", high = "royal blue") +
  theme_minimal() +
  labs(title = "Spatial Density Score (3D)", color = "Density Score")
print(p)
ggsave("spatial_density_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

# Highlight low-density cells
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = low_density_flag)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "orange")) +
  theme_minimal() +
  labs(title = "Low Spatial Density Cells", color = "Low Density")
print(p)
ggsave("spatial_density_UMAP_cut.pdf", plot=p, width = 8, height = 8, units = "in")

# Combine density and coherence measurements for exclusion criterion of suspect cells
# spatial_coherence_3D: ranges from 0 (neighbors from other clusters) to 1 (all same cluster). Higher = more spatially coherent.
# spatial_density_3D: higher = denser local neighborhood (from 1 / mean kNN distance).
# Set absolute coherence cutoff
coherence_cutoff <- 0.05  # Or any other fixed value you prefer
# Set quantile-based density cutoff (bottom 5%)
density_cutoff <- quantile(meta_df$spatial_density_3D, 0.05, na.rm = TRUE)
# Apply filtering
meta_df <- meta_df %>%
  mutate(
    high_coherence = spatial_coherence_3D_weighted >= coherence_cutoff,
    high_density   = spatial_density_3D >= density_cutoff,
    keep_cell = high_coherence & high_density
  )
table(meta_df$keep_cell)
# Visualize in UMAP space the cells to be excluded 
p <- ggplot(meta_df, aes(x = UMAP1, y = UMAP2, color = keep_cell)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "orange", "TRUE" = "darkgreen")) +
  theme_minimal() +
  labs(title = "Filtered Cells (High Coherence + Density)", color = "Retained")
print(p)
ggsave("combined_excluded_UMAP.pdf", plot=p, width = 8, height = 8, units = "in")

table(meta_df$louvain_cluster)
table(meta_df$keep_cell, meta_df$louvain_cluster, meta_df$batch)

# Filter out cells with low density and coherence scores
retained_cells <- LCNE_barseq[, meta_df$keep_cell]
dim(LCNE_barseq)
dim(retained_cells)
saveRDS(retained_cells, "LCNE_clusters_filtered_coherence_filtered.rds")

# Normalize and save LC-NE clusters filtered coherence filtered sce object
LCNE_barseq_clusters_filtered_coherence_filtered <- load_barseq(filename="LCNE_clusters_filtered_coherence_filtered.rds", from_output = TRUE)
colData(LCNE_barseq_clusters_filtered_coherence_filtered)
logcounts(LCNE_barseq_clusters_filtered_coherence_filtered) = log1p(assay(LCNE_barseq_clusters_filtered_coherence_filtered, "cpm"))/log(2)
saveRDS(LCNE_barseq_clusters_filtered_coherence_filtered,file.path(BARSEQ_OUTPUT_DIR, "LCNE_clusters_filtered_coherence_filtered_cpm_log.rds"))

############################################################# re-cluster LC-NE cluster cleaned up cells ###########################################################
# Clean up between processing steps
clear_objects_except_functions()
# Take cleaned up LC-NE cells and cluster them again
LCNE <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log.rds")
dim(LCNE)
colData(LCNE)
assayNames(LCNE)

# Drop columns from previous processing
cols_to_drop <- c("louvain_cluster", "spatial_coherence_3D_weighted", "spatial_density_3D", "cell_id_old")
colData(LCNE) <- colData(LCNE)[, !(colnames(colData(LCNE)) %in% cols_to_drop)]

# Run clustering on the data and save results to a named folder
v<-analyze_barseq(LCNE, "barseq_LC_NE_clusters_filtered_coherence_filtered")
new_barseq <- v[[1]]
clusters <- v[[2]]

#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(LCNE)$batch)  # Use colData directly

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster), 
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +  # Smaller point size
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.05, alpha = 0.6) +  # Smaller point size, alpha for density
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
  theme_minimal() +
  theme(panel.grid = element_blank()) +  # Remove grids
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNE_cells_clusters_genes_batch.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Dbh", "Th", "Ddc", "Slc18a2")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(LCNE)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 4, nrow = 1)
print(p_combined3)
ggsave("LCNE_cells_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(LCNE)) # should be true
colData(LCNE)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(LCNE, "LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")

#plot gene expression hetamap sorted by cluster
# Extract log-normalized counts
expr_matrix <- logcounts(LCNE)
# Extract cluster labels
cluster_labels <- colData(LCNE)$louvain_cluster
# Sort the expression matrix by cluster labels
sorted_indices <- order(cluster_labels)
expr_matrix_sorted <- expr_matrix[, sorted_indices]
# Update column names to reflect sorted cluster labels
colnames(expr_matrix_sorted) <- cluster_labels[sorted_indices]
# Add gene names as row names (if not already present)
rownames(expr_matrix_sorted) <- rownames(LCNE)
# Create a data frame for column annotations
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices]))  # Convert to factor
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))  # Unique row names
# Update column names of expr_matrix_sorted to match annotation_col row names
colnames(expr_matrix_sorted) <- rownames(annotation_col)
# Map the custom colors to the cluster levels
cluster_colors <- setNames(color_palette[1:length(levels(annotation_col$Cluster))], 
                           levels(annotation_col$Cluster))
annotation_colors <- list(Cluster = cluster_colors)
# Plot heatmap using pheatmap with annotations
pheatmap(expr_matrix_sorted,
         cluster_rows = FALSE,  # Disable row clustering to preserve order
         cluster_cols = FALSE,  # Disable column clustering to preserve sorting
         scale = "none",        # No scaling to match the original behavior
         color = viridis::viridis(50),  # Use viridis color palette
         show_rownames = TRUE,  # Show gene names as row labels
         show_colnames = FALSE, # Hide dense column labels
         annotation_col = annotation_col,  # Add cluster annotations
         annotation_colors = annotation_colors,  # Use custom cluster-specific colors
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,      # Adjust font size for rows
         legend = TRUE)         # Ensure legend is displayed
dev.copy(pdf, "Log-count gene expression in LCNE clusters filtered coherence filtered.pdf", width = 10, height = 12)
dev.off()

# Sort gene expression within a cluster for heatmap plotting
# Calculate total expression per cell 
total_expr <- colSums(expr_matrix)
sorted_indices_within <- c()
for (clust in levels(cluster_labels)) {
  cells_in_cluster <- which(cluster_labels == clust)
  
  # Sort cells inside cluster by total expression (descending)
  order_in_cluster <- cells_in_cluster[order(total_expr[cells_in_cluster], decreasing = TRUE)]
  
  sorted_indices_within <- c(sorted_indices_within, order_in_cluster)
}
# Subset and reorder expression matrix and annotation
expr_matrix_sorted <- expr_matrix[, sorted_indices_within]
annotation_col <- data.frame(Cluster = factor(cluster_labels[sorted_indices_within]))
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))
colnames(expr_matrix_sorted) <- rownames(annotation_col)
pheatmap(expr_matrix_sorted,
         cluster_rows = TRUE,   # Let genes be clustered to reveal patterns
         cluster_cols = FALSE,  # Keep your cell ordering intact
         scale = "none",
         color = viridis::viridis(50),
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors,
         main = "Log-count Gene Expression by Cluster",
         fontsize_row = 6,
         legend = TRUE)
dev.copy(pdf, "Log-count gene expression in LCNE  clusters filtered coherence filtered sorted.pdf", width = 10, height = 12)
dev.off()

# Calculate mean expression for each gene across clusters
clusters_labels <- clusters[["label"]]
logcounts_matrix <- logcounts(LCNE)  # Extract log-normalized counts
# Create a data frame with cluster labels and logcounts
expression_data <- as.data.frame(t(logcounts_matrix))  # Transpose to make cells rows
expression_data$cluster <- cluster_labels
# Calculate mean expression for each gene by cluster
mean_expression <- expression_data %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
# Convert mean_expression to a matrix for heatmap
mean_expression_matrix <- as.matrix(mean_expression[,-1])  # Remove cluster column
rownames(mean_expression_matrix) <- mean_expression$cluster
# Plot heatmap
pheatmap(mean_expression_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         main = "Mean Gene Expression by Cluster")
dev.copy(pdf, "LCNE_cells_clusters_filtered_coherence_filtered_gene_expr_heatmap.pdf", width = 14, height = 8)
dev.off()









###################################################################################################################################################################################
# Troubleshooting comparisons between cells for LCNE isolation combined at the outset vs each sample processed separately first and MAPseq/BARseq overlap
###################################################################################################################################################################################
################################## Check what is the divergence between unique and overlapping cells for fromFULL and fromLCNE subsetting #########################################
# Extract cell IDs
fromFULL <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(fromFULL)
cells_fromFULL <- colnames(fromFULL)
fromLCNE <- readRDS("fromLCNE_combined_LCcluster_neurons_CCFv2_uid_cpm_log_clust.rds")
dim(fromLCNE)
cells_fromLCNE <- colnames(fromLCNE)
# Find shared cells
shared_cells <- intersect(cells_fromFULL, cells_fromLCNE)
# Find unique to fromFULL
unique_to_fromFULL <- setdiff(cells_fromFULL, cells_fromLCNE)
# Find unique to fromLCNE
unique_to_fromLCNE <- setdiff(cells_fromLCNE, cells_fromFULL)

cat("Shared cells:", length(shared_cells), "\n")
cat("Unique to fromFULL:", length(unique_to_fromFULL), "\n")
cat("Unique to fromLCNE:", length(unique_to_fromLCNE), "\n")

# Check batch for unique cells in fromFULL
if (length(unique_to_fromFULL) > 0) {
  coldata_full_unique <- colData(fromFULL)[unique_to_fromFULL, ]
  batch_counts_full <- table(coldata_full_unique$batch)
  cat("Batch distribution for unique cells in fromFULL:\n")
  print(batch_counts_full)
}

# Check batch for unique cells in fromLCNE
if (length(unique_to_fromLCNE) > 0) {
  coldata_lcne_unique <- colData(fromLCNE)[unique_to_fromLCNE, ]
  batch_counts_lcne <- table(coldata_lcne_unique$batch)
  cat("Batch distribution for unique cells in fromLCNE:\n")
  print(batch_counts_lcne)
}

# Prepare data for plotting unique cells
# Extract colData for unique cells and add source
if (length(unique_to_fromFULL) > 0) {
  df_full_unique <- as.data.frame(colData(fromFULL)[unique_to_fromFULL, ])
  df_full_unique$source <- "FULL"
} else {
  df_full_unique <- data.frame()
}

if (length(unique_to_fromLCNE) > 0) {
  df_lcne_unique <- as.data.frame(colData(fromLCNE)[unique_to_fromLCNE, ])
  df_lcne_unique$source <- "LCNE"
} else {
  df_lcne_unique <- data.frame()
}
df_lcne_unique <- df_lcne_unique[, colnames(df_lcne_unique) != "cell_id_old"]
# Combine into one dataframe
unique_df <- rbind(df_full_unique, df_lcne_unique)
# Create a combined color variable for source and batch
unique_df$source_batch <- paste(unique_df$source, unique_df$batch, sep = "_")
# Ensure slice is ordered properly as a numeric factor
if (nrow(unique_df) > 0) {
  unique_df$slice <- factor(unique_df$slice, levels = sort(as.numeric(unique(unique_df$slice))))
  
  # Define color palette: two shades of pink for FULL, two shades of green for LCNE
  color_palette <- c(
    "FULL_brain3" = "#FFB6C1",  # Light pink
    "FULL_brain4" = "#DC143C",  # Dark pink
    "LCNE_brain3" = "#90EE90",  # Light green
    "LCNE_brain4" = "#006400"   # Dark green
  )
  
  # Create the scatter plot with color by source_batch
  plot <- ggplot(unique_df, aes(x = CCF_ML, y = -CCF_DV, color = source_batch)) +  # Use -CCF_DV to invert the y-axis
    geom_point(alpha = 0.9, size = 0.3) +  # Same size dots
    scale_color_manual(values = color_palette) +  # Custom color palette
    facet_wrap(~ slice, ncol = 11, scales = "fixed") +  # Arrange plots by slice, fix axes ranges
    labs(
      x = "CCF ML",
      y = "CCF DV",
      color = "Source and Batch"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 6, face = "bold"),  # Reduce slice label font size
      axis.title = element_text(size = 7),  # Reduce x and y axis label font size
      legend.position = "bottom",
      legend.title = element_text(size = 10),  # Increase legend title size
      legend.text = element_text(size = 8),    # Increase legend text size
      panel.grid = element_blank(),  # Remove grid lines (wireframes)
      panel.spacing = unit(0.1, "lines"),  # Reduce spacing between facets
      plot.margin = margin(2, 2, 2, 2)  # Reduce plot margins
    )
  print(plot)
  ggsave("unique_cells_localization_by_batch_shades.pdf", plot = plot, device = "pdf", width = 20, height = 12)
} else {
  cat("No unique cells to plot.\n")
}

############################################################# check how MAPseq cells relate to combined versus separate processing ###############################################################
# Read the CSV file into MAPseq_cells
MAPseq_cells <- read.csv("cell_top_projections_with_coords.csv", header = TRUE)
head(MAPseq_cells)

# First, get the base uids from MAPseq_cells by removing the trailing ".number"
mapseq_uids <- sub("\\..*", "", MAPseq_cells$cell_id)

# Get cell IDs (uids) from the SCE objects
cells_fromFULL <- colData(fromFULL)$uid
cells_fromLCNE <- colData(fromLCNE)$uid

# Find shared and unique cells as before
shared_cells <- intersect(cells_fromFULL, cells_fromLCNE)
unique_to_fromFULL <- setdiff(cells_fromFULL, cells_fromLCNE)
unique_to_fromLCNE <- setdiff(cells_fromLCNE, cells_fromFULL)

# Check which MAPseq uids are in each category
mapseq_in_shared <- mapseq_uids %in% shared_cells
mapseq_in_unique_FULL <- mapseq_uids %in% unique_to_fromFULL
mapseq_in_unique_LCNE <- mapseq_uids %in% unique_to_fromLCNE
num_shared <- sum(mapseq_in_shared)
num_unique_FULL <- sum(mapseq_in_unique_FULL)
num_unique_LCNE <- sum(mapseq_in_unique_LCNE)
num_not_in_any <- length(mapseq_uids) - (num_shared + num_unique_FULL + num_unique_LCNE)

cat("MAPseq cells in shared cells:", num_shared, "\n")
cat("MAPseq cells unique to FULL:", num_unique_FULL, "\n")
cat("MAPseq cells unique to LCNE:", num_unique_LCNE, "\n")
cat("MAPseq cells not in any dataset:", num_not_in_any, "\n")

# Identify MAPseq uids in each category
mapseq_in_unique_FULL <- mapseq_uids %in% unique_to_fromFULL
mapseq_in_unique_LCNE <- mapseq_uids %in% unique_to_fromLCNE

# Get the unique MAPseq uids
mapseq_unique_FULL <- mapseq_uids[mapseq_in_unique_FULL]
mapseq_unique_LCNE <- mapseq_uids[mapseq_in_unique_LCNE]

# Identify indices of MAPseq cells unique to LCNE
indices_unique_LCNE <- which(mapseq_uids %in% mapseq_unique_LCNE)
# Create dataframe for MAPseq cells unique to LCNE
mapseq_unique_LCNE_df <- MAPseq_cells[indices_unique_LCNE, ]

# Now, get batches for unique to FULL (match by uid column, not rownames)
if (length(mapseq_unique_FULL) > 0) {
  uid_full <- colData(fromFULL)$uid
  batch_full <- colData(fromFULL)$batch
  
  idx_full <- match(mapseq_unique_FULL, uid_full)     # positions in fromFULL
  batches_unique_FULL <- batch_full[idx_full]
  
  cat("Batches for MAPseq cells unique to FULL:\n")
  print(table(batches_unique_FULL, useNA = "ifany"))
} else {
  cat("No MAPseq cells unique to FULL.\n")
}

# Get batches for unique to LCNE
if (length(mapseq_unique_LCNE) > 0) {
  uid_lcne <- colData(fromLCNE)$uid
  batch_lcne <- colData(fromLCNE)$batch
  
  idx_lcne <- match(mapseq_unique_LCNE, uid_lcne)     # positions in fromLCNE
  batches_unique_LCNE <- batch_lcne[idx_lcne]
  
  cat("Batches for MAPseq cells unique to LCNE:\n")
  print(table(batches_unique_LCNE, useNA = "ifany"))
} else {
  cat("No MAPseq cells unique to LCNE.\n")
}

# Add batch information to the dataframe
mapseq_unique_LCNE_df$batch <- batches_unique_LCNE

# Subset for brain3
mapseq_brain3_df <- mapseq_unique_LCNE_df[mapseq_unique_LCNE_df$batch == "brain3", ]

# Subset for brain4
mapseq_brain4_df <- mapseq_unique_LCNE_df[mapseq_unique_LCNE_df$batch == "brain4", ]

# View or print
cat("MAPseq data for unique LCNE cells in brain3:\n")
print(mapseq_brain3_df)

cat("MAPseq data for unique LCNE cells in brain4:\n")
print(mapseq_brain4_df)

####################################################################################################################################################################################### 
# check DE genes between brain3 and 4 for combined and separate processing 
#########################################################################################################################################################################################

# Function to perform DE
perform_de <- function(sce, assay_name = "logcounts") {
  # Extract expression matrix
  expr <- assay(sce, assay_name)
  # Design matrix: batch as factor
  design <- model.matrix(~ batch, data = colData(sce))
  # Fit linear model
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  # Get top DE genes (brain4 vs brain3)
  top <- topTable(fit, coef = "batchbrain4", number = Inf, sort.by = "P")
  return(top)
}

# For fromFULL 
de_fromFULL <- perform_de(fromFULL, "logcounts")
head(de_fromFULL)
de_fromFULL$gene <- rownames(de_fromFULL)
significant <- de_fromFULL[de_fromFULL$P.Value < 0.05, ]
ggplot(de_fromFULL, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = significant, aes(x = AveExpr, y = logFC), color = "red") +
  geom_text_repel(data = significant, aes(label = gene), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "fromFULL", x = "Average Expression", y = "logFC") +
  theme_minimal()

# For fromLCNE
de_fromLCNE <- perform_de(fromLCNE, "logcounts")
head(de_fromLCNE)
de_fromLCNE$gene <- rownames(de_fromLCNE)
significant_lcne <- de_fromLCNE[de_fromLCNE$P.Value < 0.05, ]
ggplot(de_fromLCNE, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = significant_lcne, aes(x = AveExpr, y = logFC), color = "red") +
  geom_text_repel(data = significant_lcne, aes(label = gene), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "fromLCNE", x = "Average Expression", y = "logFC") +
  theme_minimal()

# Function to get stats
get_gene_stats <- function(sce, genes) {
  stats_list <- list()
  for (gene in genes) {
    raw_expr <- assay(sce, "counts")[gene, ]
    raw_by_batch <- colData(sce) %>% 
      as.data.frame() %>% 
      mutate(expr = raw_expr) %>% 
      group_by(batch) %>% 
      summarise(
        mean_raw = mean(expr, na.rm = TRUE),
        median_raw = median(expr, na.rm = TRUE),
        sd_raw = sd(expr, na.rm = TRUE)
      )
    if ("logcounts" %in% assayNames(sce)) {
      log_expr <- assay(sce, "logcounts")[gene, ]
      log_by_batch <- colData(sce) %>% 
        as.data.frame() %>% 
        mutate(expr = log_expr) %>% 
        group_by(batch) %>% 
        summarise(
          mean_log = mean(expr, na.rm = TRUE),
          median_log = median(expr, na.rm = TRUE),
          sd_log = sd(expr, na.rm = TRUE)
        )
      combined <- left_join(raw_by_batch, log_by_batch, by = "batch")
    } else {
      combined <- raw_by_batch
    }
    stats_list[[gene]] <- combined
  }
  return(stats_list)
}

genes <- c("Dbh", "Th", "Tacr3", "Slc18a2")

stats_fromFULL <- get_gene_stats(fromFULL, genes)
print(stats_fromFULL)

stats_fromLCNE <- get_gene_stats(fromLCNE, genes)
print(stats_fromLCNE)

combine_stats <- function(stats_list, dataset_name) {
  bind_rows(stats_list, .id = "gene") %>%
    mutate(dataset = dataset_name)
}

df_fromFULL <- combine_stats(stats_fromFULL, "fromFULL")
df_fromLCNE <- combine_stats(stats_fromLCNE, "fromLCNE")
combined_df <- bind_rows(df_fromFULL, df_fromLCNE)
print(combined_df)

ggplot(combined_df, aes(x = gene, y = mean_log, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ dataset) +
  labs(title = "Mean Log Expression by Gene and Batch", y = "Mean Log Counts") +
  theme_minimal()

ggplot(combined_df, aes(x = gene, y = mean_raw, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ dataset) +
  labs(title = "Mean Raw Counts by Gene and Batch", y = "Mean Raw Counts") +
  theme_minimal()

##########################################################################################################################################################################################
# Gene distribution QC + targeted correction for Dbh
# CONTEXT:
#   - BARseq data, two brains (brain3, brain4) profiling the same
#     noradrenergic neuron population with the same probe panel.
#   - Genes from sequencing cycles (Th, Slc18a2, ...) show good 
#     inter-batch agreement; brain3 is generally slightly higher.
#   - Genes from hybridization cycles (Dbh, Tacr3) show batch
#     differences — Dbh is anomalously HIGH in brain4, opposite 
#     to the general trend. Tacr3 has a modest difference.
#   - Without Dbh/Tacr3, no batch correction is needed.
#
# CORRECTION DECISION:
#   - Dbh: CORRECT. Massive batch effect (logFC=0.53, t=29.4,
#     #1 DE gene). Hybridization artifact clearly dominant.
#   - Tacr3: DO NOT CORRECT. Modest effect (logFC=0.10, t=7.9).
#     Correction in count space is too aggressive because Tacr3
#     has very low counts (median expressing = 1-2), so 
#     round(count * k) with k=0.5 zeroes out many cells,
#     collapsing detection rate from 0.294 to 0.184 and 
#     overshooting the target. Better left uncorrected.
#
# NORMALIZATION PIPELINE (custom BARseq):
#   1. convert_to_cpm(counts, total_counts = 10)
#      -> divides each cell's counts by (colSum / 10)
#   2. logcounts = log1p(cpm) / log(2) = log2(1 + cpm)
#
# CORRECTION STRATEGY:
#   Work directly from raw COUNTS for Dbh only, then 
#   re-normalize using the same pipeline to ensure consistency.
#   Scale factor estimated from median of EXPRESSING cells 
#   (counts > 0) to avoid zero-inflation bias.
##########################################################################################################################################################################################

genes     <- c("Dbh", "Th", "Tacr3", "Slc18a2")
genes_adj <- c("Dbh")  # Dbh only; Tacr3 excluded (see rationale above)

stopifnot(all(genes %in% rownames(fromLCNE)))
stopifnot("batch" %in% colnames(colData(fromLCNE)))
stopifnot("counts" %in% assayNames(fromLCNE))
stopifnot("logcounts" %in% assayNames(fromLCNE))

############################################################
# 1) Helper: make a long dataframe for plotting
############################################################
make_long_expr_df <- function(sce, genes, batch_col = "batch", log_assay_name = "logcounts") {
  stopifnot(batch_col %in% colnames(colData(sce)))
  stopifnot(all(genes %in% rownames(sce)))
  stopifnot("counts" %in% assayNames(sce))
  stopifnot(log_assay_name %in% assayNames(sce))
  
  meta <- as.data.frame(colData(sce))
  meta$cell_id <- colnames(sce)
  
  # counts
  counts_mat <- assay(sce, "counts")[genes, , drop = FALSE]
  df_counts <- as.data.frame(t(as.matrix(counts_mat)))
  df_counts$cell_id <- rownames(df_counts)
  
  counts_long <- df_counts %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "value") %>%
    left_join(meta[, c("cell_id", batch_col)], by = "cell_id") %>%
    rename(batch = all_of(batch_col)) %>%
    mutate(assay = "counts", value = as.numeric(value))
  
  # chosen log assay
  log_mat <- assay(sce, log_assay_name)[genes, , drop = FALSE]
  df_log <- as.data.frame(t(as.matrix(log_mat)))
  df_log$cell_id <- rownames(df_log)
  
  log_long <- df_log %>%
    pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "value") %>%
    left_join(meta[, c("cell_id", batch_col)], by = "cell_id") %>%
    rename(batch = all_of(batch_col)) %>%
    mutate(assay = log_assay_name, value = as.numeric(value))
  
  bind_rows(counts_long, log_long)
}

############################################################
# 2) QC plots BEFORE correction
############################################################
df_before <- make_long_expr_df(fromLCNE, genes, batch_col = "batch", log_assay_name = "logcounts")

# 2a) logcounts histograms (absolute count)
ggplot(df_before %>% filter(assay == "logcounts"),
       aes(x = value, fill = batch)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 50) +
  facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  labs(title = "fromLCNE BEFORE: logcounts distributions by batch",
       x = "logcounts", y = "Number of cells") +
  theme_minimal()

# 2b) logcounts — density (properly normalized per batch)
ggplot(df_before %>% filter(assay == "logcounts"),
       aes(x = value, fill = batch)) +
  geom_histogram(aes(y = after_stat(density)),
                 position = "identity", alpha = 0.45, bins = 50) +
  facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  labs(title = "fromLCNE BEFORE: logcounts distributions by batch (density)",
       x = "logcounts", y = "Density") +
  theme_minimal()

# 2c) log1p(raw counts) histograms
df_counts_before <- df_before %>%
  filter(assay == "counts") %>%
  mutate(value_log1p = log1p(value))

ggplot(df_counts_before, aes(x = value_log1p, fill = batch)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 50) +
  facet_wrap(~ gene, scales = "free_x", ncol = 2) +
  labs(title = "fromLCNE BEFORE: log1p(raw counts) distributions by batch",
       x = "log1p(counts)", y = "Number of cells") +
  theme_minimal()

# 2d) Tacr3 expressing cells only (logcounts > 0)
ggplot(df_before %>% filter(assay == "logcounts", gene == "Tacr3", value > 0),
       aes(x = value, fill = batch)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 50) +
  labs(title = "fromLCNE BEFORE: Tacr3 logcounts (expressing cells only)",
       x = "Tacr3 logcounts (value > 0)", y = "Number of cells") +
  theme_minimal()

# 2e) Per-cell total counts (depth/efficiency proxy)
tot <- colSums(assay(fromLCNE, "counts"))
df_tot <- data.frame(batch = colData(fromLCNE)$batch, total = tot)

ggplot(df_tot, aes(x = log1p(total), fill = batch)) +
  geom_histogram(position = "identity", alpha = 0.45, bins = 60) +
  labs(title = "fromLCNE: per-cell total counts (log1p)",
       x = "log1p(total counts per cell)", y = "Number of cells") +
  theme_minimal()

############################################################
# 3) Targeted correction: Dbh only, in RAW COUNT SPACE
#
#    Scale factor from median of EXPRESSING cells (counts > 0)
#    to avoid zero-inflation bias.
#    Direction: brain4 Dbh is anomalously high → scale DOWN.
#    k = median(brain3_expressing) / median(brain4_expressing)
#    Applied to expressing brain4 cells only; zeros stay zero.
############################################################

counts_mat <- assay(fromLCNE, "counts")
batch      <- colData(fromLCNE)$batch
idx3       <- batch == "brain3"
idx4       <- batch == "brain4"

counts_adj <- counts_mat   # copy; we only modify Dbh

min_pos <- 30  # minimum expressing cells per batch for reliable estimate

# --- Dbh correction ---
g <- "Dbh"
c3     <- counts_mat[g, idx3]
c4     <- counts_mat[g, idx4]
c3_pos <- c3[c3 > 0]
c4_pos <- c4[c4 > 0]

cat("\n--- Dbh diagnostics ---\n")
cat("brain3: n_total =", length(c3), ", n_expressing =", length(c3_pos),
    ", detect_rate =", round(length(c3_pos)/length(c3), 3),
    ", median_expressing =", median(c3_pos),
    ", mean_expressing =", round(mean(c3_pos), 2), "\n")
cat("brain4: n_total =", length(c4), ", n_expressing =", length(c4_pos),
    ", detect_rate =", round(length(c4_pos)/length(c4), 3),
    ", median_expressing =", median(c4_pos),
    ", mean_expressing =", round(mean(c4_pos), 2), "\n")

stopifnot(length(c3_pos) >= min_pos && length(c4_pos) >= min_pos)

k_dbh <- median(c3_pos) / median(c4_pos)
cat("Dbh scale factor (brain4 expressing -> brain3 expressing):", signif(k_dbh, 4), "\n")
cat("Expected effect: brain4 Dbh expressing counts multiplied by", signif(k_dbh, 4), "then rounded\n")

# Apply to expressing brain4 cells only (zeros stay zero)
pos4_dbh <- idx4 & (counts_mat[g, ] > 0)
counts_adj[g, pos4_dbh] <- round(counts_mat[g, pos4_dbh] * k_dbh)

# Sanity checks
stopifnot(all(counts_adj[g, idx4 & (counts_mat[g, ] == 0)] == 0))  # zeros preserved
stopifnot(all(counts_adj[g, idx3] == counts_mat[g, idx3]))          # brain3 unchanged

# --- Confirm all other genes are untouched ---
for (g_check in setdiff(rownames(counts_mat), "Dbh")) {
  stopifnot(identical(as.numeric(counts_adj[g_check, ]),
                      as.numeric(counts_mat[g_check, ])))
}
cat("Sanity checks passed: Dbh zeros preserved, brain3 unchanged, all other genes identical.\n")

# --- Log Tacr3 stats for the record (not corrected, but documented) ---
cat("\n--- Tacr3 (NOT corrected — documented for reference) ---\n")
g <- "Tacr3"
c3t <- counts_mat[g, idx3]; c3t_pos <- c3t[c3t > 0]
c4t <- counts_mat[g, idx4]; c4t_pos <- c4t[c4t > 0]
cat("brain3: n_expressing =", length(c3t_pos), ", detect_rate =", round(length(c3t_pos)/sum(idx3), 3),
    ", median_expressing =", median(c3t_pos), "\n")
cat("brain4: n_expressing =", length(c4t_pos), ", detect_rate =", round(length(c4t_pos)/sum(idx4), 3),
    ", median_expressing =", median(c4t_pos), "\n")
cat("Hypothetical k_tacr3 would be:", signif(median(c3t_pos)/median(c4t_pos), 4),
    "— NOT applied (too aggressive for low-count gene)\n")

# Store corrected counts
assay(fromLCNE, "counts_adj_Dbh") <- counts_adj

############################################################
# 4) Re-normalize corrected counts using the SAME pipeline
#    as the original BARseq processing:
#
#    convert_to_cpm(counts, total_counts = 10):
#      cpm_ij = counts_ij / (colSum_j / 10)
#    logcounts = log1p(cpm) / log(2) = log2(1 + cpm)
#
#    Using ORIGINAL per-cell totals as denominator so that
#    correcting one gene doesn't ripple into normalization
#    of all other genes.
############################################################

cat("\n--- Re-normalizing corrected counts ---\n")

# Original per-cell totals
totals_orig <- Matrix::colSums(counts_mat)
normalization_factor <- 10

# Replicate convert_to_cpm: divide by (colSum / total_counts)
if (is(counts_adj, "dgCMatrix")) {
  cpm_adj <- counts_adj
  cpm_adj@x <- cpm_adj@x / rep.int(totals_orig / normalization_factor, diff(cpm_adj@p))
} else {
  cpm_adj <- sweep(counts_adj, 2, totals_orig / normalization_factor, "/")
}

# log2(1 + cpm), matching: logcounts = log1p(cpm) / log(2)
logcounts_adj <- log1p(cpm_adj) / log(2)

# Store adjusted assays
assay(fromLCNE, "cpm_adj_Dbh")       <- cpm_adj
assay(fromLCNE, "logcounts_adj_Dbh")  <- logcounts_adj
cat("Stored assays: counts_adj_Dbh, cpm_adj_Dbh, logcounts_adj_Dbh\n")

# --- Verification: all genes except Dbh should have identical logcounts ---
cat("\nVerification — max logcounts difference per gene (should be 0 for non-Dbh):\n")
for (g_check in genes) {
  orig_log <- as.numeric(assay(fromLCNE, "logcounts")[g_check, ])
  adj_log  <- as.numeric(assay(fromLCNE, "logcounts_adj_Dbh")[g_check, ])
  max_diff <- max(abs(adj_log - orig_log))
  cat(" ", g_check, ":", signif(max_diff, 4),
      ifelse(g_check == "Dbh", " (corrected — expected difference)", ""), "\n")
}

# Spot-check: does re-normalization reproduce original for untouched gene?
g_test <- "Th"
orig  <- as.numeric(assay(fromLCNE, "logcounts")[g_test, 1:5])
recalc <- as.numeric(logcounts_adj[g_test, 1:5])
cat("\nSpot-check Th (first 5 cells):\n")
cat("  Original:     ", round(orig, 6), "\n")
cat("  Recalculated: ", round(recalc, 6), "\n")
cat("  Match:", all.equal(orig, recalc), "\n")

############################################################
# 5) QC plots AFTER correction
############################################################
df_after <- make_long_expr_df(fromLCNE, genes, batch_col = "batch",
                              log_assay_name = "logcounts_adj_Dbh")

# 5a) Density histograms: Dbh before vs after
plot_density_hist <- function(df, title_prefix) {
  df_log <- df %>% filter(assay != "counts")
  ggplot(df_log, aes(x = value, fill = batch)) +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity", alpha = 0.45, bins = 50) +
    facet_wrap(~ gene, scales = "free_x", ncol = 2) +
    labs(title = title_prefix, x = unique(df_log$assay), y = "Density") +
    theme_minimal()
}

plot_density_hist(df_before %>% filter(gene == "Dbh"),
                  "fromLCNE BEFORE: Dbh (density)")

plot_density_hist(df_after %>% filter(gene == "Dbh"),
                  "fromLCNE AFTER: Dbh (density)")

# 5b) All four genes AFTER — confirm Th/Slc18a2/Tacr3 unchanged
plot_density_hist(df_after,
                  "fromLCNE AFTER correction: all 4 genes (density)")

# 5c) Dbh expressing-only before vs after
plot_gene_pos <- function(df, gene_name, title) {
  df_log <- df %>% filter(assay != "counts", gene == gene_name, value > 0)
  ggplot(df_log, aes(x = value, fill = batch)) +
    geom_histogram(aes(y = after_stat(density)),
                   position = "identity", alpha = 0.45, bins = 50) +
    labs(title = title, x = unique(df_log$assay), y = "Density") +
    theme_minimal()
}

plot_gene_pos(df_before, "Dbh", "Dbh expressing-only BEFORE (logcounts)")
plot_gene_pos(df_after,  "Dbh", "Dbh expressing-only AFTER (logcounts_adj_Dbh)")

############################################################
# 6) Summary stats before vs after
############################################################
get_gene_stats2 <- function(sce, genes, counts_assay = "counts", log_assay = "logcounts") {
  stats_list <- list()
  for (gene in genes) {
    raw_expr <- assay(sce, counts_assay)[gene, ]
    raw_by_batch <- colData(sce) %>%
      as.data.frame() %>%
      mutate(expr = as.numeric(raw_expr)) %>%
      group_by(batch) %>%
      summarise(
        n_cells = n(),
        n_expressing = sum(expr > 0),
        detect_rate = round(mean(expr > 0), 3),
        mean_raw = round(mean(expr), 3),
        median_raw = median(expr),
        median_raw_expressing = median(expr[expr > 0]),
        sd_raw = round(sd(expr), 3),
        .groups = "drop"
      )
    
    log_expr <- assay(sce, log_assay)[gene, ]
    log_by_batch <- colData(sce) %>%
      as.data.frame() %>%
      mutate(expr = as.numeric(log_expr)) %>%
      group_by(batch) %>%
      summarise(
        mean_log = round(mean(expr), 4),
        median_log = round(median(expr), 4),
        sd_log = round(sd(expr), 4),
        .groups = "drop"
      )
    
    stats_list[[gene]] <- left_join(raw_by_batch, log_by_batch, by = "batch")
  }
  stats_list
}

cat("BEFORE correction\n")
stats_before <- get_gene_stats2(fromLCNE, genes,
                                counts_assay = "counts",
                                log_assay = "logcounts")
print(bind_rows(stats_before, .id = "gene"))


cat("AFTER correction (Dbh only)\n")
stats_after <- get_gene_stats2(fromLCNE, genes,
                               counts_assay = "counts_adj_Dbh",
                               log_assay = "logcounts_adj_Dbh")
print(bind_rows(stats_after, .id = "gene"))

# Side-by-side comparison
df_compare <- bind_rows(
  bind_rows(stats_before, .id = "gene") %>% mutate(stage = "before"),
  bind_rows(stats_after,  .id = "gene") %>% mutate(stage = "after")
)
print(df_compare)

############################################################
# 7) Re-run DE on corrected assay — KEY VALIDATION
#    Expectation: Dbh should drop from #1 DE gene to non-significant.
#    Tacr3/Th/Slc18a2 should be essentially unchanged.
############################################################
de_fromLCNE_adj <- perform_de(fromLCNE, assay_name = "logcounts_adj_Dbh")
de_fromLCNE_adj$gene <- rownames(de_fromLCNE_adj)

cat("DE results BEFORE correction (top 10):\n")
print(head(de_fromLCNE, 10))


cat("DE results AFTER correction — Dbh only (top 10):\n")
print(head(de_fromLCNE_adj, 10))

cat("\nTarget genes before vs after:\n")
cat("--- BEFORE ---\n")
print(de_fromLCNE[de_fromLCNE$gene %in% c("Dbh", "Tacr3", "Th", "Slc18a2"), ])
cat("--- AFTER ---\n")
print(de_fromLCNE_adj[de_fromLCNE_adj$gene %in% c("Dbh", "Tacr3", "Th", "Slc18a2"), ])

# MA plot for corrected data
significant_adj <- de_fromLCNE_adj[de_fromLCNE_adj$P.Value < 0.05, ]
ggplot(de_fromLCNE_adj, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = significant_adj, aes(x = AveExpr, y = logFC), color = "red") +
  geom_text_repel(data = significant_adj, aes(label = gene), size = 3, max.overlaps = 10) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "fromLCNE AFTER Dbh correction",
       x = "Average Expression", y = "logFC") +
  theme_minimal()

####################################################################################################################################################################################
# UMAP + clustering on Dbh-corrected fromLCNE object
#
# Before running this, ensure you have:
#   - fromLCNE with assays: counts_adj_Dbh, cpm_adj_Dbh, logcounts_adj_Dbh
#   - The analyze_barseq function available in your environment
#
# Key decision: analyze_barseq likely uses logcounts internally
# for PCA/UMAP. We need to temporarily set the "logcounts" assay
# to our corrected version so that the function picks it up,
# OR pass the corrected assay name if analyze_barseq supports it.
#
# Strategy below: create a working copy, swap in the corrected
# logcounts as the default "logcounts" assay, run the analysis,
# then visualize.
####################################################################################################################################################################################

# Create working copy from the Dbh-corrected fromLCNE
Dbh_corrected_test <- fromLCNE

# Replace the logcounts assay with the Dbh-corrected version
# so that analyze_barseq (which reads "logcounts") uses the corrected values
assay(Dbh_corrected_test, "logcounts_original") <- assay(Dbh_corrected_test, "logcounts")
assay(Dbh_corrected_test, "logcounts") <- assay(Dbh_corrected_test, "logcounts_adj_Dbh")

# Sanity check: confirm the swap
cat("Logcounts assay now points to Dbh-corrected values.\n")
cat("Dbh mean_log brain3:", round(mean(assay(Dbh_corrected_test, "logcounts")["Dbh", colData(Dbh_corrected_test)$batch == "brain3"]), 4), "\n")
cat("Dbh mean_log brain4:", round(mean(assay(Dbh_corrected_test, "logcounts")["Dbh", colData(Dbh_corrected_test)$batch == "brain4"]), 4), "\n")
cat("Th  mean_log brain3:", round(mean(assay(Dbh_corrected_test, "logcounts")["Th", colData(Dbh_corrected_test)$batch == "brain3"]), 4), "(should be unchanged)\n")

# Run analysis (PCA + UMAP + clustering)
v <- analyze_barseq(Dbh_corrected_test, "Dbh_corrected_test")
new_barseq <- v[[1]]
clusters <- v[[2]]

# Set up color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)

# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))

# Create plotting data frames
plot_data_cluster <- data.frame(
  UMAP1 = umap_data[,1],
  UMAP2 = umap_data[,2],
  cluster = factor(clusters[["label"]])
)
plot_data_genes <- data.frame(
  UMAP1 = umap_data[,1],
  UMAP2 = umap_data[,2],
  TotalGenes = total_genes
)
plot_data_batch <- data.frame(
  UMAP1 = umap_data[,1],
  UMAP2 = umap_data[,2],
  batch = colData(Dbh_corrected_test)$batch
)

# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
  dplyr::arrange(TotalGenes)

# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))

# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = color_palette) +
  geom_text(data = centroid_data, aes(x = x, y = y, label = cluster),
            colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +
  theme_minimal() +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")

# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
  geom_point(size = 0.05) +
  scale_color_gradient(low = "grey", high = "magenta") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")

# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.05, alpha = 0.1) +
  scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")

# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNE_cells_clusters_genes_batch_Dbh_corrected.pdf",
       plot = p_combined1, device = "pdf", width = 14, height = 8)

# Gene expression on UMAP — now using corrected logcounts
# Include Dbh and Tacr3 alongside the other markers to visualize
# how the correction affects their spatial distribution
genes <- c("Dbh", "Th", "Tacr3", "Slc18a2", "Ddc", "Dlk1")
plots_genes <- list()

for (gene in genes) {
  # Use the corrected logcounts (which is now the default "logcounts" assay)
  plot_data_cluster[[gene]] <- logcounts(Dbh_corrected_test)[gene, ]
  
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}

p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 3, nrow = 2)
print(p_combined3)
ggsave("LCNE_cells_gene_expression_Dbh_corrected.pdf",
       plot = p_combined3, device = "pdf", width = 10, height = 14)

# Include cluster assignment data in the SCE object
table(clusters$label)
all(clusters[['sample']] == colnames(Dbh_corrected_test))  # should be TRUE
colData(Dbh_corrected_test)$louvain_cluster <- as.factor(clusters[["label"]])

# Save — note: this object has "logcounts" = corrected, "logcounts_original" = uncorrected
saveRDS(Dbh_corrected_test, "Dbh_corrected_fromLCNE_CCFv2_uid_cpm_log_clust.rds")


# Export corrected logcounts for Dbh, Th, Slc18a2
# Uses the corrected assay (logcounts_adj_Dbh) from fromLCNE

genes_export <- c("Dbh", "Th", "Slc18a2")

# Extract corrected logcounts for selected genes
log_mat <- assay(fromLCNE, "logcounts_adj_Dbh")[genes_export, , drop = FALSE]

# Transpose so cells are rows, genes are columns
df_export <- as.data.frame(t(as.matrix(log_mat)))

# Add batch info
df_export$batch <- colData(fromLCNE)$batch

# Cell IDs are already the rownames from colnames(fromLCNE)
cat("First few rows:\n")
head(df_export)

# Write to CSV with row.names = TRUE to retain cell IDs
write.csv(df_export, file = "Dbh_Th_Slc18a2_logcounts_adj_Dbh.csv", row.names = TRUE)
cat("Exported", nrow(df_export), "cells x", ncol(df_export), "columns to Dbh_Th_Slc18a2_logcounts_adj_Dbh.csv\n")











############################################################# CCA integrate brain3 and 4 gene expression to batch correct ###############################################################
# Convert SCE to Seurat
brain3_seurat <- as.Seurat(brain3_subset, counts = "counts", data = NULL)  # Use appropriate assay slot
brain4_seurat <- as.Seurat(brain4, counts = "counts", data = NULL)

# Add batch metadata
brain3_seurat$batch <- "brain3"
brain4_seurat$batch <- "brain4"

# Prepare individual objects (normalize and find features)
brain3_seurat <- NormalizeData(brain3_seurat)
brain3_seurat <- FindVariableFeatures(brain3_seurat, selection.method = "vst", nfeatures = 100)  # Reduce to <=102

brain4_seurat <- NormalizeData(brain4_seurat)
brain4_seurat <- FindVariableFeatures(brain4_seurat, selection.method = "vst", nfeatures = 100)

brain3_seurat$batch <- "brain3"
brain4_seurat$batch <- "brain4"

# Find integration anchors (reduce dims to avoid over-fitting)
anchors <- FindIntegrationAnchors(object.list = list(brain3_seurat, brain4_seurat), 
                                  dims = 1:10, reduction = "cca")  # Or 1:20 max

# Integrate
integrated_seurat <- IntegrateData(anchorset = anchors, dims = 1:10)
# Set default assay
DefaultAssay(integrated_seurat) <- "integrated"
integrated_seurat <- ScaleData(integrated_seurat)
integrated_seurat <- RunPCA(integrated_seurat, npcs = 10)  # Max ~102
integrated_seurat <- RunUMAP(integrated_seurat, reduction = "pca", dims = 1:10)
integrated_seurat <- FindNeighbors(integrated_seurat, reduction = "pca", dims = 1:10)
integrated_seurat <- FindClusters(integrated_seurat, resolution = 0.1)

p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "seurat_clusters") +
  ggtitle("Integrated Clusters")
p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "batch") +
  ggtitle("Integrated Batch")
plot_grid(p1, p2, ncol = 2)

# Plot clustree
resolutions <- seq(0.1, 1.0, by = 0.1)
for (res in resolutions) {
  integrated_seurat  <- FindClusters(integrated_seurat , resolution = res)
  integrated_seurat [[paste0("res.", res)]] <- integrated_seurat $seurat_clusters
}
clustree(integrated_seurat, prefix = "res.") + ggtitle("Cluster Stability Across Resolutions")
# Plot UMAPs at different resolutions
p1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "res.0.1") + ggtitle("Resolution 0.1")
p2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "res.0.2") + ggtitle("Resolution 0.2")
p3 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "res.0.3") + ggtitle("Resolution 0.3")
p4 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "res.0.3") + ggtitle("Resolution 0.4")
plot_grid(p1, p2, p3, p4, ncol = 2)

# Export integrated expression matrix (genes x cells)
expr <- GetAssayData(integrated_seurat, assay = "integrated", layer = "data")
head(expr)
dim(expr)
write.csv(expr, "integrated_gene_expression.csv")

# Export metadata (includes cell info like barcodes, clusters, batch, projections if added)
meta <- integrated_seurat@meta.data
head(meta)
dim(meta)
write.csv(meta, "integrated_metadata.csv")

# Alternatively, save the full Seurat object for external use
saveRDS(integrated_seurat, "integrated_seurat.rds")

############################################################# drop Dbh and Tacr3 from analyses ###############################################################
# # Load separately processed datasets and concatenate them
# brain3 <- readRDS("/results/BARseq_780345/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
# dim(brain3)
# brain4 <- readRDS("/results/BARseq_780346/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
# dim(brain4)
# # Examine why gene vectors are different length, only keep shared genes
# # Extract gene names from both objects
# genes_brain3 <- rownames(brain3)
# genes_brain4 <- rownames(brain4)
# # Find genes unique to brain3 (not in brain4)
# unique_to_brain3 <- setdiff(genes_brain3, genes_brain4)
# # Find genes unique to brain4 (not in brain3)
# unique_to_brain4 <- setdiff(genes_brain4, genes_brain3)
# # Find common genes
# common_genes <- intersect(genes_brain3, genes_brain4)
# # Print lengths for summary
# cat("Genes unique to brain3:", length(unique_to_brain3), "\n")
# cat("Genes unique to brain4:", length(unique_to_brain4), "\n")
# cat("Common genes:", length(common_genes), "\n")
# # View first few unique genes
# head(unique_to_brain3)
# head(unique_to_brain4)
# # Subset brain3 to only include genes (rows) shared with brain4
# brain3_subset <- brain3[rownames(brain3) %in% rownames(brain4), ]
# # Verify the new dimensions
# dim(brain3_subset)
# # Concatenate along columns (cells)
# combined_sce <- cbind(brain3_subset, brain4)
# # Verify dimensions
# dim(combined_sce)
# # Add a batch column
# colData(combined_sce)$batch <- c(rep("brain3", ncol(brain3_subset)), rep("brain4", ncol(brain4)))

head(combined_sce)
gene_names <- rownames(combined_sce)
genes_to_remove <- c("Dbh", "Tacr3")
combined_sce <- combined_sce[!gene_names %in% genes_to_remove, ]
dim(combined_sce)

Dbh_Tacr3_dropped_test <- combined_sce
v<-analyze_barseq(Dbh_Tacr3_dropped_test, "Dbh_Tacr3_dropped_test")
new_barseq <- v[[1]]
clusters <- v[[2]]
#visualize UMAP of samples - set color palette
n_clusters <- length(unique(clusters[["label"]]))
color_palette <- get_cluster_colors(n_clusters)
# Extract UMAP coordinates and total gene counts
umap_data <- reducedDim(new_barseq, "UMAP")
total_genes <- colSums(counts(new_barseq))
# Create plotting data frames
plot_data_cluster <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], cluster = factor(clusters[["label"]]))
plot_data_genes <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], TotalGenes = total_genes)
plot_data_batch <- data.frame(UMAP1 = umap_data[,1], UMAP2 = umap_data[,2], batch = colData(Dbh_Tacr3_dropped_test)$batch)
# Sort data for TotalGenes to plot higher values on top
plot_data_genes <- plot_data_genes %>%
dplyr::arrange(TotalGenes)
# Calculate cluster centroids
centroid_data <- plot_data_cluster %>%
dplyr::group_by(cluster) %>%
dplyr::summarise(x = mean(UMAP1), y = mean(UMAP2))
# Plot UMAP with clusters
p1 <- ggplot(plot_data_cluster, aes(x = UMAP1, y = UMAP2, color = cluster)) +
geom_point(size = 0.05) +  # Smaller point size
scale_color_manual(values = color_palette) +
geom_text(data = centroid_data, aes(x = x, y = y, label = cluster),
colour = "black", vjust = 1.6, hjust = 0.5, size = 3.5) +  # Add cluster numbers
theme_minimal() +
theme(legend.position = "none", panel.grid = element_blank()) +  # Remove grids
labs(title = "UMAP plot", x = "UMAP1", y = "UMAP2", color = "Cluster")
# Plot UMAP with total gene counts
p2 <- ggplot(plot_data_genes, aes(x = UMAP1, y = UMAP2, color = TotalGenes)) +
geom_point(size = 0.05) +  # Smaller point size
scale_color_gradient(low = "grey", high = "magenta") +
theme_minimal() +
theme(panel.grid = element_blank()) +  # Remove grids
labs(title = "Total counts", x = "UMAP1", y = "UMAP2", color = "Gene Count")
# Plot UMAP with batch
p3 <- ggplot(plot_data_batch, aes(x = UMAP1, y = UMAP2, color = batch)) +
geom_point(size = 0.01, alpha = 0.1) +  # Smaller point size, alpha for density
scale_color_manual(values = c("brain3" = "blue", "brain4" = "green")) +  # Custom colors
theme_minimal() +
theme(panel.grid = element_blank()) +  # Remove grids
labs(title = "Batch", x = "UMAP1", y = "UMAP2", color = "Batch")
# 3 panels (clusters, total genes, batch)
p_combined1 <- grid.arrange(grobs = list(p1, p2, p3), ncol = 3)
print(p_combined1)
ggsave("LCNE_cells_clusters_genes_batch_Dbh-Tacr3_dropped.pdf", plot = p_combined1, device = "pdf", width = 14, height = 8)


genes <- c("Th", "Ddc", "Slc18a2", "Dlk1", "Eya2","Pdyn")
plots_genes <- list()
for (gene in genes) {
  plot_data_cluster[[gene]] <- logcounts(Dbh_Tacr3_dropped_test)[gene, ]
  # Sort data for each gene to plot higher expression values on top
  plot_data_sorted <- plot_data_cluster %>%
    dplyr::arrange(!!sym(gene))
  
  p <- ggplot(plot_data_sorted, aes(x = UMAP1, y = UMAP2, color = !!sym(gene))) +
    geom_point(size = 0.02) +  # Smaller point size
    scale_color_gradient(low = "cyan", high = "red") +
    theme_minimal() +
    theme(panel.grid = element_blank()) +  # Remove grids
    labs(title = paste(gene, "Expression"), x = "UMAP1", y = "UMAP2", color = "Logcounts")
  
  plots_genes[[length(plots_genes) + 1]] <- p
}
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 3, nrow = 2)
print(p_combined3)
ggsave("LCNE_cells_gene_expression_Dbh-Tacr3_dropped.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

#include cluster assignment data to the original SingleCellExperiment object
table(clusters$label)
all(clusters[['sample']]==colnames(Dbh_Tacr3_dropped_test)) # should be true
colData(Dbh_Tacr3_dropped_test)$louvain_cluster <- as.factor(clusters[["label"]])

saveRDS(Dbh_Tacr3_dropped_test, "Dbh_Tacr3_dropped_test_CCFv2_uid_cpm_log_clust.rds")
