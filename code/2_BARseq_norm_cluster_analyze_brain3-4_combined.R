################# BARseq Analysis Pipeline ##################
# Normalization, clustering, and spatial coherence analyses to isolate LC-NE neurons 
# have been performed on brain#3 and brain#4 datasets independently
# now bringing the data together to perform transcriptomic group identification, harmonize normalized matrices if necessary

# Load functions
source("~/capsule/code/1_BARseq_analyses_functions_brain3-4_combined.R")

# Set working directory
setwd(BARSEQ_OUTPUT_DIR)

############################################################################################################################################################################################################
# load LC-NE isolated, normalized gene expression data, concatenate the matrices
brain3 <- readRDS("/scratch/BARseq_780345/combined_neurons_clust_CCFv2_uid_cpm_log.rds") 

brain4 <- readRDS("/scratch/BARseq_780346/combined_neurons_clust_CCFv2_uid_cpm_log.rds") 
dim(brain4)

# Examine why gene vectors are different length, only keep shared genes
# Extract gene names from both objects
genes_brain3 <- rownames(brain3)
genes_brain4 <- rownames(brain4)
# Find genes unique to brain3 (not in brain4)
unique_to_brain3 <- setdiff(genes_brain3, genes_brain4)
# Find genes unique to brain4 (not in brain3)
unique_to_brain4 <- setdiff(genes_brain4, genes_brain3)
# Find common genes
common_genes <- intersect(genes_brain3, genes_brain4)
# Print lengths for summary
cat("Genes unique to brain3:", length(unique_to_brain3), "\n")
cat("Genes unique to brain4:", length(unique_to_brain4), "\n")
cat("Common genes:", length(common_genes), "\n")
# View first few unique genes 
head(unique_to_brain3)
head(unique_to_brain4)
# Subset brain3 to only include genes (rows) shared with brain4
brain3_subset <- brain3[rownames(brain3) %in% rownames(brain4), ]
# Verify the new dimensions
dim(brain3_subset)

# Concatenate along columns (cells)
combined_sce <- cbind(brain3_subset, brain4)
# Verify dimensions
dim(combined_sce) 
# Add a batch column
colData(combined_sce)$batch <- c(rep("brain3", ncol(brain3_subset)), rep("brain4", ncol(brain4)))

# # Extract colData from combined_sce and convert to data.frame
# metadata_df <- as.data.frame(colData(combined_sce))
# # Save to CSV in the working directory
# write.csv(metadata_df, file = "/scratch/BARseq_780345-780346_combined/LCNE_combined_metadata.csv", row.names = TRUE)

####################################### perform initial clustering and check for batch effects ###############################################################################################
# Check if clustering analysis already exists
output_dir <- "analysis/barseq_all_QCed_cells"
umap_file <- file.path(output_dir, "umap.csv")
cluster_file <- file.path(output_dir, "cluster.csv")
annot_file <- file.path(output_dir, "cluster_annotation.csv")

if (dir.exists(output_dir) && file.exists(umap_file) && file.exists(cluster_file) && file.exists(annot_file)) {
  cat("All analysis results already exist for barseq_all_QCed_cells. Skipping.\n")
} else if (dir.exists(output_dir) && file.exists(umap_file)) {
  cat("UMAP exists but clustering incomplete. Loading UMAP and running PCA-based clustering...\n")
  # Load UMAP
  umap_mat <- as.matrix(read_csv(umap_file))
  rownames(umap_mat) <- colnames(combined_sce)
  reducedDim(combined_sce, "UMAP") <- umap_mat
  # Run PCA and clustering
  v <- analyze_barseq_pca_only(combined_sce, "barseq_all_QCed_cells")
} else {
  cat("Running full clustering analysis for barseq_all_QCed_cells...\n")
  v <- analyze_barseq(combined_sce, "barseq_all_QCed_cells")  # Full pipeline if nothing exists
}

# plot the results from clustering on the whole dataset
# load saved umap and cluster information for data including all the neurons
file_path <- file.path(BARSEQ_OUTPUT_DIR, "analysis/barseq_all_QCed_cells/umap.csv")
umap_data <- read.csv(file_path)
file_path <- file.path(BARSEQ_OUTPUT_DIR, "analysis/barseq_all_QCed_cells/cluster.csv")
clusters <- read.csv(file_path)
x<-umap_data[['UMAP1']]
y<-umap_data[['UMAP2']]

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
  labs(title = "Total Genes", x = "UMAP1", y = "UMAP2", color = "Gene Count")

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
write.csv(gene_expression_summary, file = "gene_expression_summary_by_cluster.csv", row.names = FALSE)

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
write.csv(combination_counts, file = "gene_coexpression_summary_LC_cluster.csv", row.names = FALSE)

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
logcounts(LC) = log1p(cpm(LC))/log(2)
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
  labs(title = "Total genes", x = "UMAP1", y = "UMAP2", color = "Gene Count")

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
ggsave("LCcluster_gene_expression.pdf", plot = p_combined3, device = "pdf", width = 10, height = 8)  

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
LCNE_barseq <- LC_barseq[, colData(LC_barseq)$louvain_cluster %in% c(1,2,3)]
dim(LCNE_barseq)
saveRDS(LCNE_barseq, "LCNE_cluster_neurons_CCFv2_uid.rds")

# Normalize and save LC-NE cluster object
LCNE_barseq <- load_barseq(filename="LCNE_cluster_neurons_CCFv2_uid.rds", from_output = TRUE)
colData(LCNE_barseq)
logcounts(LCNE_barseq) = log1p(cpm(LCNE_barseq))/log(2)
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
  labs(title = "Total genes", x = "UMAP1", y = "UMAP2", color = "Gene Count")

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
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
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

umap_path <- "/scratch/BARseq_780345-780346_combined/analysis/barseq_LC_NE_cells/umap.csv"
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
logcounts(LCNE_barseq_clusters_filtered_coherence_filtered) = log1p(cpm(LCNE_barseq_clusters_filtered_coherence_filtered))/log(2)
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
cols_to_drop <- c("louvain_cluster", "spatial_coherence_3D_weighted", "spatial_density_3D")
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
  labs(title = "Total genes", x = "UMAP1", y = "UMAP2", color = "Gene Count")

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
p_combined3 <- grid.arrange(grobs = plots_genes, ncol = 2, nrow = 2)
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

######################################################## compare the LC-NE final overlap between combined and separate processing ##############################################################################
fromFULL <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")

# Load separately processed datasets and concatenate them
brain3 <- readRDS("/scratch/BARseq_780345/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(brain3)
brain4 <- readRDS("/scratch/BARseq_780346/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(brain4)
# Examine why gene vectors are different length, only keep shared genes
# Extract gene names from both objects
genes_brain3 <- rownames(brain3)
genes_brain4 <- rownames(brain4)
# Find genes unique to brain3 (not in brain4)
unique_to_brain3 <- setdiff(genes_brain3, genes_brain4)
# Find genes unique to brain4 (not in brain3)
unique_to_brain4 <- setdiff(genes_brain4, genes_brain3)
# Find common genes
common_genes <- intersect(genes_brain3, genes_brain4)
# Print lengths for summary
cat("Genes unique to brain3:", length(unique_to_brain3), "\n")
cat("Genes unique to brain4:", length(unique_to_brain4), "\n")
cat("Common genes:", length(common_genes), "\n")
# View first few unique genes
head(unique_to_brain3)
head(unique_to_brain4)
# Subset brain3 to only include genes (rows) shared with brain4
brain3_subset <- brain3[rownames(brain3) %in% rownames(brain4), ]
# Verify the new dimensions
dim(brain3_subset)

# # Data export for Shuonan
# cnt3 <- assay(brain3_subset, "counts")
# cnt4 <- assay(brain4, "counts")
# colnames(cnt3) <- make.unique(colnames(cnt3), sep="__dup")
# colnames(cnt4) <- make.unique(colnames(cnt4), sep="__dup")
# colnames(cnt3) <- paste0("brain3|", colnames(cnt3))
# colnames(cnt4) <- paste0("brain4|", colnames(cnt4))
# preview_genes <- rownames(cnt3)[1:10]
# preview_cells3 <- colnames(cnt3)[1:5]
# preview_cells4 <- colnames(cnt4)[1:5]
# preview3 <- as.data.frame(as.matrix(t(cnt3[preview_genes, preview_cells3])))
# preview4 <- as.data.frame(as.matrix(t(cnt4[preview_genes, preview_cells4])))
# preview3
# preview4
# dt3 <- as.data.table(as.matrix(t(cnt3)), keep.rownames = "cell_id")
# library(data.table)
# dt3 <- as.data.table(as.matrix(t(cnt3)), keep.rownames = "cell_id")
# dt4 <- as.data.table(as.matrix(t(cnt4)), keep.rownames = "cell_id")
# fwrite(dt3, "brain3_counts_cells_by_genes.csv.gz")
# fwrite(dt4, "brain4_counts_cells_by_genes.csv.gz")
# cd3 <- as.data.table(as.data.frame(colData(brain3)))
# cd3[, cell_id := colnames(cnt3)]
# fwrite(cd3, "brain3_colData.csv.gz")
# cd4 <- as.data.table(as.data.frame(colData(brain4)))
# cd4[, cell_id := colnames(cnt4)]
# fwrite(cd4, "brain4_colData.csv.gz")
# View(dt4)
# View(dt3)
# View(cd3)
# View(cd4)


# Concatenate along columns (cells)
combined_sce <- cbind(brain3_subset, brain4)
# Verify dimensions
dim(combined_sce)
# Add a batch column
colData(combined_sce)$batch <- c(rep("brain3", ncol(brain3_subset)), rep("brain4", ncol(brain4)))

fromLCNE <- combined_sce

dim(fromFULL)
dim(fromLCNE)

# Extract cell IDs
cells_fromFULL <- colnames(fromFULL)
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
cells_fromFULL <- colnames(fromFULL)
cells_fromLCNE <- colnames(fromLCNE)

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

# Now, get batches for unique to FULL
if (length(mapseq_unique_FULL) > 0) {
  batches_unique_FULL <- colData(fromFULL)[mapseq_unique_FULL, ]$batch
  cat("Batches for MAPseq cells unique to FULL:\n")
  print(table(batches_unique_FULL))
} else {
  cat("No MAPseq cells unique to FULL.\n")
}

# Get batches for unique to LCNE
if (length(mapseq_unique_LCNE) > 0) {
  batches_unique_LCNE <- colData(fromLCNE)[mapseq_unique_LCNE, ]$batch
  cat("Batches for MAPseq cells unique to LCNE:\n")
  print(table(batches_unique_LCNE))
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

############################################################# check DE genes between brain3 and 4 for combined and separate processing ###############################################################
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
head(de_fromFULL)  # View top results
de_fromFULL$gene <- rownames(de_fromFULL)
significant <- de_fromFULL[de_fromFULL$P.Value < 0.05, ]  # Adjust threshold as needed
ggplot(de_fromFULL, aes(x = AveExpr, y = logFC)) +
  geom_point(alpha = 0.5, color = "grey") +
  geom_point(data = significant, aes(x = AveExpr, y = logFC), color = "red") +  # Highlight significant
  geom_text_repel(data = significant, aes(label = gene), size = 3, max.overlaps = 10) +  # Label significant
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
    # Raw counts
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
    # Logcounts (if available)
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

# For fromFULL
stats_fromFULL <- get_gene_stats(fromFULL, genes)
print(stats_fromFULL)

# For fromLCNE
stats_fromLCNE <- get_gene_stats(fromLCNE, genes)
print(stats_fromLCNE)

# Combine stats into a single dataframe
combine_stats <- function(stats_list, dataset_name) {
  bind_rows(stats_list, .id = "gene") %>%
    mutate(dataset = dataset_name)
}

df_fromFULL <- combine_stats(stats_fromFULL, "fromFULL")
df_fromLCNE <- combine_stats(stats_fromLCNE, "fromLCNE")
combined_df <- bind_rows(df_fromFULL, df_fromLCNE)
print(combined_df)  # View as dataframe

# Visualization: Boxplot of logcounts per gene per batch (across datasets)
ggplot(combined_df, aes(x = gene, y = mean_log, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ dataset) +
  labs(title = "Mean Log Expression by Gene and Batch", y = "Mean Log Counts") +
  theme_minimal()

# For raw counts: Replace mean_log with mean_raw
ggplot(combined_df, aes(x = gene, y = mean_raw, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ dataset) +
  labs(title = "Mean Raw Counts by Gene and Batch", y = "Mean Raw Counts") +
  theme_minimal()

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
# brain3 <- readRDS("/scratch/BARseq_780345/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
# dim(brain3)
# brain4 <- readRDS("/scratch/BARseq_780346/LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
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
labs(title = "Total genes", x = "UMAP1", y = "UMAP2", color = "Gene Count")
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
