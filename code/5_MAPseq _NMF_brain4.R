# Load functions which handle pre-processing or organizing of the data
source("~/capsule/code/01_loaders_brain4.R")

#set working directory 
setwd('/results/BARseq_780346/')

dat <- load_data()
names(dat)
mat_ordered <- dat$mat_log_ordered
inRH_lookup <- dat$inRH_lookup
metadata <- dat$metadata

#############################################################################################################################
# # Subset data to only include cells which pass visual QC for segmentation accuracy
# good_cells <- read.csv("LC_visualQC_barcoded_cells.csv")
# # Get the uids of good barcoded cells
# good_uids <- good_cells$uid[good_cells$good_barcoded == 1]
# # Extract the base uid from rownames of mat_ordered (remove .1, .2, etc.)
# base_rownames <- sub("\\.[0-9]+$", "", rownames(mat_ordered))
# # Subset mat_ordered to only include rows with base uid in good_uids
# mat_ordered <- mat_ordered[base_rownames %in% good_uids, ]
# # Subset inRH_lookup in the same way
# inRH_lookup <- inRH_lookup[base_rownames %in% good_uids, , drop = FALSE]
#############################################################################################################################

# Plot heatmap of log-norm projections
heatmap.2(mat_ordered, 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "none", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.6, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,
          main = "Log-transformed projections map")
dev.copy(pdf, "heatmap_lognorm_roisort_clust_RH+LH.pdf", width = 16, height = 12)
dev.off()

# Plot heatmap of log-norm projections, sort rows according to soma location
# Sort rows so inRH==1 comes first
order_idx <- order(-inRH_lookup$inRH)
mat_ordered_sorted <- mat_ordered[order_idx, ]
inRH_lookup_sorted <- inRH_lookup[order_idx, ]
# Create a color vector for inRH (red for 1, blue for 0)
inRH_colors <- ifelse(inRH_lookup_sorted$inRH == 1, "red", "blue")
# Plot the heatmap with row annotation
heatmap.2(mat_ordered_sorted,
          distfun = function(x) proxy::dist(x, method = "cosine"),
          hclustfun = function(x) hclust(x, method = "ward.D"),
          scale = "none",
          Colv = NA,
          Rowv = NA,
          dendrogram = "none",
          trace = "none",
          density.info = "histogram",
          cexRow = 0.3,
          cexCol = 0.6,
          keysize = 1,
          margins = c(10, 8),
          col = viridis,
          RowSideColors = inRH_colors,
          main = "Log-transformed projections map")
# Add a legend for inRH annotation
legend("topright",
       legend = c("inRH = 1", "inRH = 0"),
       fill = c("red", "blue"),
       border = NA,
       bty = "n",
       cex = 0.8)
dev.copy(pdf, "heatmap_lognorm_roisort_clust_RH+LH_soma.pdf", width = 16, height = 12)
dev.off()


#############################################################################################################################
#based on cell belonging to either RH or LH soma location, plot its projections patterns to RH and LH
mat_ordered_1 <- mat_ordered[inRH_lookup$inRH == 1, ]
mat_ordered_0 <- mat_ordered[inRH_lookup$inRH == 0, ]

# Plot the heatmap for RH_soma = 1
heatmap.2(mat_ordered_1, 
          distfun = function(x) proxy::dist(x, method = "cosine"), 
          hclustfun = function(x) hclust(x, method = "ward.D"), 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "row", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.7, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,  
          main = "Log-transformed projections map for RH_soma = 1")
dev.copy(pdf, "heatmap_lognorm_roisort_clust_RH+LH_somasRH.pdf", width = 12, height = 10)
dev.off()

# Plot the heatmap for RH_soma = 0
heatmap.2(mat_ordered_0, 
          distfun = function(x) proxy::dist(x, method = "cosine"), 
          hclustfun = function(x) hclust(x, method = "ward.D"), 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "row", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.7, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,  
          main = "Log-transformed projections map for RH_soma = 0")
dev.copy(pdf, "heatmap_lognorm_roisort_clust_RH+LH_somasLH.pdf", width = 12, height = 10)
dev.off()

#################################################################################################################
#function to perform NMF on projections data and save relevant heatmap plots 
perform_nmf <- function(mat_ordered, rank, nrun, seed = NULL) {
  require(NMF)
  require(gplots)
  require(viridisLite)
  
  cat("Running NMF with rank =", rank, ", nrun =", nrun, "\n")
  nmf_result <- nmf(mat_ordered, rank = rank, method = "lee", nrun = nrun, seed = seed)
  
  W <- coef(nmf_result)
  H <- basis(nmf_result)
  Hsquared <- H^2
  
  ## ---- W heatmap: show + save ----
  heatmap.2(W,
            scale = "none",
            Colv = TRUE,
            Rowv = NA,
            dendrogram = "none",
            trace = "none",
            density.info = "histogram",
            cexRow = 0.7,
            cexCol = 0.15,
            keysize = 0.8,
            margins = c(5, 5),
            col = viridisLite::viridis,
            main = "Cells by factors")
  
  pdf(paste0("NMF_W_cells_rank_", rank, "_nrun", nrun, ".pdf"),
      width = 10, height = 10)
  heatmap.2(W,
            scale = "none",
            Colv = TRUE,
            Rowv = NA,
            dendrogram = "none",
            trace = "none",
            density.info = "histogram",
            cexRow = 0.7,
            cexCol = 0.15,
            keysize = 0.8,
            margins = c(5, 5),
            col = viridisLite::viridis,
            main = "Cells by factors")
  dev.off()
  
  ## ---- H heatmap: show + save ----
  heatmap.2(H,
            scale = "none",
            Colv = NA,
            Rowv = NA,
            dendrogram = "none",
            trace = "none",
            density.info = "histogram",
            cexRow = 0.25,
            cexCol = 0.7,
            keysize = 0.8,
            margins = c(5, 5),
            col = viridisLite::viridis,
            main = "Factors by ROI")
  
  pdf(paste0("NMF_H_regions_rank_", rank, "_nrun", nrun, ".pdf"),
      width = 10, height = 10)
  heatmap.2(H,
            scale = "none",
            Colv = NA,
            Rowv = NA,
            dendrogram = "none",
            trace = "none",
            density.info = "histogram",
            cexRow = 0.4,
            cexCol = 0.7,
            keysize = 0.8,
            margins = c(5, 5),
            col = viridisLite::viridis,
            main = "Factors by ROI")
  dev.off()
  
  ## ---- Consensus: show + save ----
  consensus_matrix <- consensus(nmf_result)
  
  consensusmap(nmf_result)  # show on screen
  
  pdf(paste0("Consensus_rank", rank, "_nrun", nrun, ".pdf"),
      width = 10, height = 10)
  consensusmap(nmf_result)  # draw directly into PDF
  dev.off()
  
  return(list(W = W, H = H, Hsquared = Hsquared,
              consensus_matrix = consensus_matrix,
              nmf_result = nmf_result))
}
############################################ ipsi and contralateral NMF ##################################################################
# Create ipsi-contra matrix using the updated loader function, use log norm data as input
ipsi_contra <- create_ipsi_contra_from_log(mat_ordered, inRH_lookup)

#plot the heatmap to visualize ipsi and contralateral patterns
numeric_matrix <- data.matrix(ipsi_contra)
heatmap.2(numeric_matrix, 
          distfun = function(x) proxy::dist(x, method = "cosine"), 
          hclustfun = function(x) hclust(x, method = "ward.D"), 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "none", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.6, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,  
          main = "Log-transformed projections map")
dev.copy(pdf, "heatmap_lognorm_roisort_clust_ipsi-contra.pdf", width = 16, height = 12)
dev.off()

# Transpose for NMF  as NMF seems to internally assume that columns are samples
t_ipsi_contra <- t(ipsi_contra)
t_ipsi_contra_clean <- t_ipsi_contra[
  rowSums(t_ipsi_contra) != 0,
  colSums(t_ipsi_contra) != 0,
  drop = FALSE
]

#perform NMF on projections data by calling the function
result_rank2 <- perform_nmf(t_ipsi_contra_clean, rank = 2, nrun = 20, seed = 12345)
result_rank3 <- perform_nmf(t_ipsi_contra_clean, rank = 3, nrun = 20, seed = 12345)
result_rank4 <- perform_nmf(t_ipsi_contra_clean, rank = 4, nrun = 20, seed = 12345)
result_rank5 <- perform_nmf(t_ipsi_contra_clean, rank = 5, nrun = 20, seed = 12345)
result_rank6 <- perform_nmf(t_ipsi_contra_clean, rank = 6, nrun = 20, seed = 12345)
result_rank7 <- perform_nmf(t_ipsi_contra_clean, rank = 7, nrun = 20, seed = 12345)

##########################################################################################################################################
# W: factors x cells (from NMF result)
W <- result_rank4$W
n_factors <- nrow(W)
n_cells   <- ncol(W)

# Hard assignment of each cell to its dominant factor - extract cell IDs and their corresponding factor assignments
# For each cell (column), find the factor with the largest loading
factor_per_cell <- apply(W, 2, which.max)

cell_factors <- data.frame(
  cellID = colnames(W),
  factor = factor_per_cell,
  stringsAsFactors = FALSE
)
cell_factors <- cell_factors[order(cell_factors$factor), ]
table(cell_factors$factor)

# Add anatomical interpretation of the factors
# collapse factors into anterior and posterior projecting, and annotate major projeciton groups
cell_factors <- cell_factors %>%
  mutate(
    proj_pattern = ifelse(factor %in% c(4), "Posterior", "Anterior"),
    proj_target  = case_when(
      factor == 4 ~ "medulla-SP",
      factor == 2 ~ "olf-ant_ctx",
      factor == 1 ~ "midbrain-hindbrain",
      factor == 3 ~ "ctx-posterior_hippocampus",
      TRUE        ~ NA_character_
    )
  )

# Strip suffix (.1, .2, etc.) to recover the base cell ID
cell_factors$split_cellID <- sub("\\.[0-9]+.*$", "", as.character(cell_factors$cellID))
write.csv(cell_factors, "NMF_ids_factors.csv", row.names = FALSE)

# examine Wnorm and Wmax to assess fidelity of factor assignment for a given cell
#  Cell-centered normalization:
#  For each cell, factors sum to 1.
#  Interpretable as "factor membership of each cell".
W_cellnorm <- sweep(W, 2, colSums(W), FUN = "/")
W_cellnorm[!is.finite(W_cellnorm)] <- 0  # protect against zero-sum columns

# Factor-centered normalization:
#  For each factor, cell weights sum to 1.
#  Interpretable as "distribution of this factor across cells".
W_factornorm <- sweep(W, 1, rowSums(W), FUN = "/")
W_factornorm[!is.finite(W_factornorm)] <- 0

# Factor-centered relative-to-max:
# For each factor, divide by the strongest cell for that factor.
# Interpretable as "relative strength compared to the strongest cell in that factor".
row_max <- apply(W_factornorm, 1, max)
W_factor_relmax <- sweep(W_factornorm, 1, row_max, FUN = "/")
W_factor_relmax[!is.finite(W_factor_relmax)] <- 0


# Heatmap of cell-centered membership:
# How cleanly does each cell belong to a single factor
membership_mat <- t(W_cellnorm)   # dimensions: cells x factors
heatmap.2(membership_mat,
          scale = "none",
          Colv = NA,              # do not cluster factors (only 5)
          Rowv = TRUE,            # cluster cells if you like
          dendrogram = "row",
          trace = "none",
          main = "Cell membership in NMF factors",
          cexRow = 0.2,
          cexCol = 1)

dev.copy(pdf, "heatmap_W_cellnorm.pdf", width = 8, height = 8)
dev.off()

# Heatmap of factor-centered relative-to-max:
# Is each factor carried by a few cells or many cells
heatmap.2(W_factor_relmax,
          scale = "none",
          Colv = TRUE,              # you can set TRUE if you want to cluster cells here too
          Rowv = NA,              # only 5 factors, no clustering needed
          dendrogram = "col",
          trace = "none",
          main = "Factor strength across cells\nnormalized per factor",
          cexCol = 0.2,
          cexRow = 1)

dev.copy(pdf, "heatmap_W_factor_relmax.pdf", width = 8, height = 8)
dev.off()

# check how factor loading relates to transcriptomic identity
# ensure that only data from properly segmented cells is being utilized here
good_cells <- read.csv("LC_visualQC_barcoded_cells.csv")
good_uids <- good_cells$uid[good_cells$good_barcoded == 1]
# Remove suffixes from cell IDs in W
base_cell_ids <- sub("\\..*$", "", colnames(W))  # remove everything after first dot  sub("\\.[0-9]+$", "", colnames(W))
# Logical vector: TRUE if cell is in good_uids
keep_cols <- base_cell_ids %in% good_uids
# Subset W to only good cells
W_qc <- W[, keep_cols, drop = FALSE]
# Update colnames to match original
colnames(W_qc) <- colnames(W)[keep_cols]

cp <- t(W_qc)                      # get cells x patterns
cp_normalized <- as.data.frame(cp/rowSums(cp)) # normalize per cell to interpret rows as membership of that cell to each pattern.
#find per-transcriptomic type avegrage membership to patterns
cp_normalized$BarSeq_cluster <- sub(".*\\.", "", rownames(cp_normalized))
cp_normalized <- cp_normalized %>% arrange(BarSeq_cluster)
aggregated_df <- cp_normalized %>%
  group_by(BarSeq_cluster) %>%
  summarise(across(starts_with("V"), ~mean(., na.rm = TRUE)), .groups = 'drop')
# Define the color palette
color_palette <- c("#dcbeff", "#aaffc3", "#fabed4", "#e6194B", "#ffe119", 
                   "#469990", "#fffac8", "#f032e6", "#ffd8b1", "#bfef45", 
                   "#3cb44b", "#42d4f4", "#911eb4", "#4363d8", "#f58231",
                   "#800000")
# Reshape the data to long format
df_long <- reshape2::melt(aggregated_df, id.vars = "BarSeq_cluster")
# Create the plot
p <- ggplot(df_long, aes(x = variable, y = value, fill = BarSeq_cluster)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major.x = element_blank(),  
        panel.grid.minor.x = element_blank()) +  
  labs(x = "NMF factor", y = "Mean factor membership", fill = "BarSeq class") +
  scale_fill_manual(values = color_palette)
print(p)
ggsave("BarSeqclass_vs_NMFgroup.pdf", plot = p, width = 6, height = 6, dpi = 600)


##########################################################################################################################################
# Tests association between factor loadings and transcriptomic clusters
# Extract cluster information
clusters <- sub(".*\\.", "", colnames(W_qc))
cluster_counts <- table(clusters)
print(cluster_counts)

# Check if we have enough samples per cluster
if (any(cluster_counts < 10)) {
  warning("Some clusters have very few cells - results may be unreliable")
}

# Initialize results storage
p_values <- numeric(nrow(W_qc))
test_types <- character(nrow(W_qc))

print("Running statistical tests for each factor...")
# Test each factor
for (i in 1:nrow(W_qc)) {
  factor_loadings <- W_qc[i, ]
  
  # Check for sufficient variation
  if (var(factor_loadings) < 1e-10) {
    print(paste("Factor", i, "has no variation - skipping"))
    p_values[i] <- NA
    test_types[i] <- "No variation"
    next
  }
  
  # Check if we have multiple clusters
  if (length(unique(clusters)) < 2) {
    print("Only one cluster found - cannot perform statistical test")
    p_values[i] <- NA
    test_types[i] <- "Single cluster"
    next
  }
  
  # Test normality (for smaller samples)
  if (length(factor_loadings) <= 5000) {
    shapiro_p <- shapiro.test(factor_loadings)$p.value
    is_normal <- shapiro_p > 0.05
  } else {
    # For large samples, default to non-parametric
    is_normal <- FALSE
  }
  
  if (is_normal) {
    # Use ANOVA for normally distributed data
    aov_result <- aov(factor_loadings ~ clusters)
    p_values[i] <- summary(aov_result)[[1]][["Pr(>F)"]][1]
    test_types[i] <- "ANOVA"
    
    # Post-hoc test if significant
    if (p_values[i] < 0.05) {
      print(paste("Factor", i, "- ANOVA p-value:", format(p_values[i], scientific = TRUE)))
      print("TukeyHSD post-hoc results:")
      posthoc <- TukeyHSD(aov_result)
      print(posthoc)
      cat("\n")
    }
    
  } else {
    # Use Kruskal-Wallis for non-normal data
    kw_result <- kruskal.test(factor_loadings ~ clusters)
    p_values[i] <- kw_result$p.value
    test_types[i] <- "Kruskal-Wallis"
    
    # Post-hoc test if significant
    if (p_values[i] < 0.05) {
      print(paste("Factor", i, "- Kruskal-Wallis p-value:", format(p_values[i], scientific = TRUE)))
      print("Pairwise Wilcoxon post-hoc results:")
      
      # Pairwise Wilcoxon with Bonferroni correction
      pairwise_result <- pairwise.wilcox.test(factor_loadings, clusters, 
                                              p.adjust.method = "bonferroni",
                                              paired = FALSE)
      print(pairwise_result)
      cat("\n")
    }
  }
}
# Apply multiple testing correction across all factors
p_adjusted <- p.adjust(p_values, method = "bonferroni")
# Create summary results
results_df <- data.frame(
  Factor = 1:nrow(W_qc),
  Test_Type = test_types,
  P_value = p_values,
  P_adjusted = p_adjusted,
  Significant_raw = p_values < 0.05,
  Significant_corrected = p_adjusted < 0.05
)
# Remove rows with NA p-values for cleaner output
results_summary <- results_df[!is.na(results_df$P_value), ]
print(results_summary)
n_sig_corrected <- sum(results_summary$Significant_corrected)
# Save results
write.csv(results_df, "factor_cluster_associations.csv", row.names = FALSE)


# Additional summary by cluster
print("CLUSTER SUMMARY:")
for (cluster in names(cluster_counts)) {
  cluster_cells <- clusters == cluster
  if (sum(cluster_cells) > 0) {
    print(paste("Cluster", cluster, "- n =", sum(cluster_cells)))
    for (i in 1:nrow(W_qc)) {
      cluster_mean <- mean(W_qc[i, cluster_cells])
      cluster_sd <- sd(W_qc[i, cluster_cells])
      print(paste("  Factor", i, "mean:", round(cluster_mean, 4), 
                  "± SD:", round(cluster_sd, 4)))
    }
    cat("\n")
  }
}

# Create visualization
print("Creating visualization...")
df <- data.frame(
  Category = sub(".*\\.", "", colnames(W_qc)),
  t(W_qc)
)
# Convert to long format
df_long <- reshape2::melt(df, id.vars = "Category")
# Create improved boxplot
p <- ggplot(df_long, aes(x = Category, y = value, fill = Category)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.3, size = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c(
    "1" = "#dcbeff",
    "2" = "#aaffc3",
    "3" = "#fabed4",
    "4" = "#e6194B"
  )) +
  facet_wrap(~ variable, scales = "free",
             labeller = labeller(variable = function(x) paste("Factor", gsub("V", "", x)))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_rect(fill = "lightgray", color = "black"),
        legend.position = "none") +
  labs(
    x = "Transcriptomic Cluster",
    y = "Factor Loading",
    title = "Factor Loadings by Transcriptomic Cluster",
    subtitle = paste("Statistical tests completed -", n_sig_corrected,
                     "significant associations after correction")
  )
print(p)
ggsave("BarSeqclust_vs_NMF_loadings.pdf", plot = p, width = 12, height = 8)

##########################################################################################################################################
# HEATMAP ANALYSIS
# Use existing clusters variable and create mean matrix
cluster_names <- sort(unique(clusters))
n_factors <- nrow(W_qc)
# Calculate mean factor loadings (reusing cluster extraction from above)
mean_matrix <- matrix(0, nrow = n_factors, ncol = length(cluster_names))
rownames(mean_matrix) <- paste("Factor", 1:n_factors)
colnames(mean_matrix) <- paste("Cluster", cluster_names)

for (i in 1:n_factors) {
  for (j in 1:length(cluster_names)) {
    cluster_cells <- clusters == cluster_names[j]
    mean_matrix[i, j] <- mean(W_qc[i, cluster_cells])
  }
}

# Create normalized versions
# Row-normalize: for each factor, distribution across clusters
row_sums <- rowSums(mean_matrix)
row_sums[row_sums == 0] <- NA  # avoid division by zero
mean_matrix_row_norm <- sweep(mean_matrix, 1, row_sums, "/") * 100
mean_matrix_row_norm[!is.finite(mean_matrix_row_norm)] <- 0

# Column-normalize: for each cluster, distribution across factors
col_sums <- colSums(mean_matrix)
col_sums[col_sums == 0] <- NA  # avoid division by zero
mean_matrix_col_norm <- sweep(mean_matrix, 2, col_sums, "/") * 100
mean_matrix_col_norm[!is.finite(mean_matrix_col_norm)] <- 0

# Heatmap function
create_heatmap <- function(data_matrix, title, subtitle, value_name, add_percent = FALSE) {
  heatmap_data <- reshape2::melt(data_matrix)
  colnames(heatmap_data) <- c("Factor", "Cluster", "Value")
  
  label_text <- if(add_percent) {
    paste0(round(heatmap_data$Value, 1), "%")
  } else {
    round(heatmap_data$Value, 1)
  }
  
  ggplot(heatmap_data, aes(x = Cluster, y = Factor, fill = Value)) +
    geom_tile(color = "white", size = 0.5) +
    geom_text(aes(label = label_text), color = "black", size = 3) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = median(heatmap_data$Value),
                         name = value_name) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title, subtitle = subtitle,
         x = "Transcriptomic Cluster", y = "NMF Factor")
}

# Generate heatmaps
h1 <- create_heatmap(mean_matrix, "Raw Factor Loadings: connectivity × transcriptomics association", 
                     "Absolute magnitude", "Mean\nLoading")
h2 <- create_heatmap(mean_matrix_row_norm, "Factor Composition by Cluster", 
                     "Which clusters dominate each projection pattern?", 
                     "% of\nFactor", add_percent = TRUE)
h3 <- create_heatmap(mean_matrix_col_norm, "Cluster Projection Profiles", 
                     "How does each cluster distribute its projections?", 
                     "% of\nCluster", add_percent = TRUE)

print(h1); 
print(h2);
print(h3)

ggsave("heatmap_raw.pdf", h1, width = 8, height = 6)
ggsave("heatmap_factor_composition.pdf", h2, width = 8, height = 6)
ggsave("heatmap_cluster_profiles.pdf", h3, width = 8, height = 6)


# SUMMARY TABLE
# Identify which FACTOR IDs are significant after correction
sig_factor_ids <- results_summary$Factor[results_summary$Significant_corrected]
# Map adjusted p-values to all factors (some may be NA if skipped earlier)
p_adj_all <- results_summary$P_adjusted[match(1:n_factors, results_summary$Factor)]
# Precompute strongest / weakest clusters and loadings
max_vals <- apply(mean_matrix, 1, max)
min_vals <- apply(mean_matrix, 1, min)

strongest_cluster <- apply(mean_matrix, 1, function(x) colnames(mean_matrix)[which.max(x)])
weakest_cluster   <- apply(mean_matrix, 1, function(x) colnames(mean_matrix)[which.min(x)])

# Fold-difference with protection against division by zero
fold_diff <- ifelse(min_vals > 0, round(max_vals / min_vals, 2), NA)
# Construct summary table
summary_comprehensive <- data.frame(
  Factor = 1:n_factors,
  Projection_Target = c(
    "midbrain-hindbrain",
    "olf-ant_ctx",
    "ctx-posterior_hippocampus",
    "medulla-SP"
  ),
  Statistically_Significant = ifelse(1:n_factors %in% sig_factor_ids, "YES", "NO"),
  Test_P_Value = p_adj_all,
  Strongest_Cluster = strongest_cluster,
  Strongest_Loading = round(max_vals, 1),
  Weakest_Cluster = weakest_cluster,
  Weakest_Loading = round(min_vals, 1),
  Fold_Difference = fold_diff,
  stringsAsFactors = FALSE
)

print(summary_comprehensive)
write.csv(summary_comprehensive, "comprehensive_factor_cluster_summary.csv", row.names = FALSE)

##########################################################################################################################################
########plot gene expression based on NMF group assignment to identify transcriptomic trends################
# --------- GENE EXPRESSION ANALYSIS BASED ON NMF GROUPS (QC-passing cells only) ---------
#  Filter cell_factors to only QC-passing cells
good_cells <- read.csv("LC_visualQC_barcoded_cells.csv")
good_uids <- good_cells$uid[good_cells$good_barcoded == 1]
cell_factors_qc <- cell_factors[cell_factors$split_cellID %in% good_uids, ]

# Load LCNEneurons
LCNEneurons <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")

# Subset to NMF cells that passed QC
subset <- match(cell_factors_qc$split_cellID, colnames(LCNEneurons))
subset <- subset[!is.na(subset)]
LCNEneurons_subset <- LCNEneurons[, subset]

# Ensure correct order for NMF info
cell_ids <- colnames(LCNEneurons_subset)

# Assign NMF info by matching cell_ids to cell_factors_qc$split_cellID
colData(LCNEneurons_subset)$NMF_factor <- cell_factors_qc$factor[match(cell_ids, cell_factors_qc$split_cellID)]
colData(LCNEneurons_subset)$NMF_group  <- cell_factors_qc$proj_pattern[match(cell_ids, cell_factors_qc$split_cellID)]
colData(LCNEneurons_subset)$NMF_target <- cell_factors_qc$proj_target[match(cell_ids, cell_factors_qc$split_cellID)]

# Extract log-normalized counts
expr_matrix <- logcounts(LCNEneurons_subset)

# Extract NMF factor labels
nmf_labels <- colData(LCNEneurons_subset)$NMF_factor

# Sort the expression matrix by NMF factor labels
sorted_indices <- order(nmf_labels)
expr_matrix_sorted <- expr_matrix[, sorted_indices]

# Add gene names as row names (if not already present)
rownames(expr_matrix_sorted) <- rownames(LCNEneurons_subset)

# Create a data frame for column annotations
annotation_col <- data.frame(NMF_Factor = factor(nmf_labels[sorted_indices]))
rownames(annotation_col) <- paste0("Cell_", seq_len(ncol(expr_matrix_sorted)))

# Update column names of expr_matrix_sorted to match annotation_col row names
colnames(expr_matrix_sorted) <- rownames(annotation_col)

# Map the custom colors to the NMF factor levels
color_palette <- c("#469990", "#f032e6", "#3cb44b", "#42d4f4", "#911eb4")
nmf_levels <- levels(annotation_col$NMF_Factor)
cluster_colors <- setNames(color_palette[1:length(nmf_levels)], nmf_levels)
annotation_colors <- list(NMF_Factor = cluster_colors)

# Plot heatmap using pheatmap with annotations
pheatmap::pheatmap(expr_matrix_sorted,
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   scale = "none",
                   color = viridis::viridis(50),
                   show_rownames = TRUE,
                   show_colnames = FALSE,
                   annotation_col = annotation_col,
                   annotation_colors = annotation_colors,
                   main = "Log-count Gene Expression by NMF Factor",
                   fontsize_row = 6,
                   legend = TRUE)
dev.copy(pdf, "geneexpr_heatmap_logcounts_ipsi_contra_nmf.pdf", width = 10, height = 12)
dev.off()

########plot gene expression based on NMF group assignment to identify transcriptomic trends using most variable genes################
# Calculate mean expression for each factor
mean_expr_by_factor <- matrix(0, nrow = nrow(expr_matrix_sorted), ncol = 4)
rownames(mean_expr_by_factor) <- rownames(expr_matrix_sorted)
colnames(mean_expr_by_factor) <- paste("Factor", 1:4)

# Also calculate sample sizes for each factor
factor_sizes <- numeric(4)
for(f in 1:4) {
  factor_cells <- which(nmf_labels[sorted_indices] == f)
  factor_sizes[f] <- length(factor_cells)
  if(length(factor_cells) > 0) {
    mean_expr_by_factor[, f] <- rowMeans(expr_matrix_sorted[, factor_cells, drop = FALSE], na.rm = TRUE)
  }
}
print(setNames(factor_sizes, paste("Factor", 1:4)))

# Debug data ranges
print("Data range diagnostics:")
print(paste("Individual cell expression range:", round(min(expr_matrix_sorted, na.rm = TRUE), 2), "to", round(max(expr_matrix_sorted, na.rm = TRUE), 2)))

# Find most variable genes across factors
gene_vars_across_factors <- apply(mean_expr_by_factor, 1, var, na.rm = TRUE)
gene_means_across_factors <- rowMeans(mean_expr_by_factor, na.rm = TRUE)

print(paste("Gene mean range (across factors):", round(min(gene_means_across_factors, na.rm = TRUE), 2), "to", round(max(gene_means_across_factors, na.rm = TRUE), 2)))
print(paste("Gene variance range:", round(min(gene_vars_across_factors, na.rm = TRUE), 4), "to", round(max(gene_vars_across_factors, na.rm = TRUE), 4)))

# Adaptive filtering - use data-driven thresholds
mean_threshold <- quantile(gene_means_across_factors, 0.25, na.rm = TRUE)  # 25th percentile # reasonable expression level
var_threshold <- quantile(gene_vars_across_factors, 0.7, na.rm = TRUE)     # 70th percentile # high variability genes

print(paste("Using adaptive mean threshold (25th percentile):", round(mean_threshold, 3)))
print(paste("Using adaptive variance threshold (70th percentile):", round(var_threshold, 4)))

# Filter for expressed and variable genes
expressed_variable_genes <- names(gene_vars_across_factors)[
  gene_vars_across_factors > var_threshold & 
    gene_means_across_factors > mean_threshold
]

print(paste("Number of genes meeting criteria:", length(expressed_variable_genes)))

# Adaptive gene selection
min_genes <- 20
max_genes <- min(75, nrow(mean_expr_by_factor))

if(length(expressed_variable_genes) < min_genes) {
  print(paste("Too few genes, using top", max_genes, "most variable..."))
  expressed_variable_genes <- names(sort(gene_vars_across_factors, decreasing = TRUE)[1:max_genes])
} else if(length(expressed_variable_genes) > max_genes) {
  print(paste("Too many genes, selecting top", max_genes, "by variance..."))
  gene_vars_subset <- gene_vars_across_factors[expressed_variable_genes]
  expressed_variable_genes <- names(sort(gene_vars_subset, decreasing = TRUE)[1:max_genes])
}

print(paste("Final gene count for analysis:", length(expressed_variable_genes)))

# Take final gene set
top_genes_final <- expressed_variable_genes

# Create enhanced annotation
mean_expr_subset <- mean_expr_by_factor[top_genes_final, ]
mean_expr_zscore <- t(scale(t(mean_expr_subset)))

# Add projection target info to column names
proj_targets <- c("midbrain-hindbrain",
                  "olf-ant_ctx",
                  "ctx-posterior_hippocampus",
                  "medulla-SP")
colnames(mean_expr_zscore) <- paste0("Factor ", 1:4, "\n(", proj_targets, ")")

# Enhanced plot
pheatmap::pheatmap(mean_expr_zscore,
                   cluster_rows = TRUE,
                   cluster_cols = FALSE,
                   scale = "none",
                   color = colorRampPalette(c("navy", "blue", "white", "red", "darkred"))(100),
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   main = "Gene Expression Signatures by NMF Projection Factor",
                   subtitle = paste("Top", length(top_genes_final), "genes (Z-scored mean expression)"),
                   fontsize_row = 7,
                   fontsize_col = 9,
                   breaks = seq(-2.5, 2.5, length.out = 101),
                   legend = TRUE,
                   border_color = "grey90")
dev.copy(pdf, "geneexpr_heatmap_logcounts_ipsi_contra_nmf_HVG.pdf", width = 8, height = 8)
dev.off()

# Save the results
write.csv(mean_expr_zscore, "gene_expression_by_nmf_factor_zscore.csv")

# Create summary of top genes per factor
print("\nTop upregulated genes per factor:")
for(f in 1:4) {
  factor_scores <- mean_expr_zscore[, f]
  top_genes_factor <- names(sort(factor_scores, decreasing = TRUE)[1:10])
  cat(sprintf("\nFactor %d (%s):\n", f, proj_targets[f]))
  cat(paste(top_genes_factor, collapse = ", "))
}
# Optional: Look for factor-specific gene sets
# Find genes most specific to each factor
factor_specific_genes <- list()

for(f in 1:4) {
  # Genes highly expressed in this factor vs others
  factor_specific <- mean_expr_zscore[, f] > 1.5 & 
    apply(mean_expr_zscore[, -f], 1, max) < 0.5
  factor_specific_genes[[f]] <- rownames(mean_expr_zscore)[factor_specific]
  
  cat(sprintf("\nFactor %d specific genes (%s): %d genes\n", 
              f, proj_targets[f], length(factor_specific_genes[[f]])))
  if(length(factor_specific_genes[[f]]) > 0) {
    cat(paste(factor_specific_genes[[f]], collapse = ", "))
  }
}

# Save factor-specific genes
factor_specific_df <- data.frame(
  Factor = rep(1:4, sapply(factor_specific_genes, length)),
  Projection_Target = rep(proj_targets, sapply(factor_specific_genes, length)),
  Gene = unlist(factor_specific_genes)
)
write.csv(factor_specific_df, "factor_specific_genes.csv", row.names = FALSE)

# Save top 10 genes per factor
top_genes_df <- data.frame(
  Factor = rep(1:4, each = 10),
  Projection_Target = rep(proj_targets, each = 10),
  Gene = unlist(lapply(1:4, function(f) {
    names(sort(mean_expr_zscore[, f], decreasing = TRUE)[1:10])
  }))
)
write.csv(top_genes_df, "top_genes_per_factor.csv", row.names = FALSE)


# Check how NMF relates to BarSeq_subclass (QC-passing cells only)
df <- as.data.frame(colData(LCNEneurons_subset))
contingency_table <- table(df$louvain_cluster, df$NMF_factor)
color_palette <- c("#dcbeff","#aaffc3", "#fabed4",  "#e6194B")
barplot(contingency_table, main="BarSeq_subclass vs NMF_factor", xlab="NMF_factor", ylab="number of cells", col=color_palette, 
        legend = rownames(contingency_table),args.legend = list(x = "topright"))
dev.copy(pdf, "BarSeq_subclass_vs_NMF_factor.pdf",
         width = 7, height = 7)
dev.off()

contingency_table <- table(df$louvain_cluster, df$NMF_group)
barplot(contingency_table, main="BarSeq_subclass vs NMF_group", xlab="NMF_group", ylab="number of cells", col=color_palette, 
        legend = rownames(contingency_table),args.legend = list(x = "topright"))
dev.copy(pdf, "BarSeq_subclass_vs_NMF_group.pdf",
         width = 7, height = 7)
dev.off()

contingency_table <- table(df$louvain_cluster, df$NMF_target)
barplot(contingency_table, main="BarSeq_subclass vs NMF_target", xlab="NMF_target", ylab="number of cells", col=color_palette, 
        legend = rownames(contingency_table),args.legend = list(x = "topright"))
dev.copy(pdf, "BarSeq_subclass_vs_NMF_target.pdf",
         width = 7, height = 7)
dev.off()

######################### implement plotting where bars indicate NMF belonging to A and P groups ################################
################################################### pheatmap sorted implementation #####################################################
library(pheatmap)
library(viridis)

# Use the already filtered and annotated objects:
# LCNEneurons_subset: SingleCellExperiment with only QC-passing, NMF-assigned cells
# cell_factors_qc: data.frame with NMF assignments for those cells

# Extract logcounts matrix for only these cells
expr_df <- logcounts(LCNEneurons_subset)

# Get the NMF group assignment for each cell (column)
# Ensure the order matches the columns of expr_df
cell_ids <- colnames(expr_df)
nmf_group <- colData(LCNEneurons_subset)$NMF_group

# Create column annotation data frame
col_anno <- data.frame(factor = as.factor(nmf_group))
rownames(col_anno) <- cell_ids

# Check group sizes
print(table(col_anno$factor))  # Should show counts for "Anterior" and "Posterior"

# Split by factor
expr_df_1 <- expr_df[, col_anno$factor == "Anterior", drop = FALSE]
expr_df_2 <- expr_df[, col_anno$factor == "Posterior", drop = FALSE]

# Cluster only if enough cells
get_hclust_order <- function(df, type) {
  if (type == "row") {
    d <- dist(as.matrix(df), method = "euclidean")
  }
  if (type == "col") {
    d <- dist(t(as.matrix(df)), method = "euclidean")
  }
  hclust_ob <- hclust(d, method = "ward.D2")
  ord <- hclust_ob$order
  return(ord)
}

if (ncol(expr_df_1) >= 2) {
  col_order_1 <- get_hclust_order(expr_df_1, type = "col")
  expr_df_1 <- expr_df_1[, col_order_1, drop = FALSE]
}
if (ncol(expr_df_2) >= 2) {
  col_order_2 <- get_hclust_order(expr_df_2, type = "col")
  expr_df_2 <- expr_df_2[, col_order_2, drop = FALSE]
}
row_order <- get_hclust_order(expr_df, type = "row")
expr_df_1 <- expr_df_1[row_order, , drop = FALSE]
expr_df_2 <- expr_df_2[row_order, , drop = FALSE]

# Combine
expr_df_reordered <- cbind(expr_df_1, expr_df_2)
colnames(expr_df_reordered) <- make.unique(colnames(expr_df_reordered))

# New column annotation for reordered matrix
factor <- c(rep("Anterior", ncol(expr_df_1)), rep("Posterior", ncol(expr_df_2)))
col_anno = data.frame(factor = as.factor(factor))
rownames(col_anno) <- colnames(expr_df_reordered)

# Remove genes with <10% nonzero
rows_ <- rowSums(expr_df_reordered > 0) > 0.1 * ncol(expr_df_reordered)

# Plot
pheatmap(
  expr_df_reordered[rows_, ],
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  annotation_col = col_anno,
  color = viridis(100),
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 3
)
dev.copy(pdf, "logcounts expression by projection pattern pheatmap.pdf", width = 13, height = 12)
dev.off()
