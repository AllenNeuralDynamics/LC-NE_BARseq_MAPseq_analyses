##########################################################################################################################################################################################################
#################### GLOBAL FDR correction cross-correlogram analyses #################
# Load functions which handle pre-processing or organizing of the data
source("~/capsule/code/01_loaders_brain3.R")

#set working directory 
setwd('/results/BARseq_780345/')

dat <- load_data()
names(dat)
raw_matrix <- dat$proj_matrix_raw
log_matrix <- dat$mat_log_ordered      # For visualization
rownorm_matrix <- dat$mat_rownorm_ordered  # For abundance analysis
stopifnot(identical(rownames(rownorm_matrix), dat$inRH_lookup$row_id))

# Improved assess_data_quality function with CSV export and fixed sparsity
assess_data_quality <- function(matrix_data, matrix_name = "Matrix", save_csv = TRUE) {
  cat("\n=== Data Quality Assessment for", matrix_name, "===\n")
  
  # Region coverage (how many cells project to each region)
  region_coverage <- colSums(matrix_data > 0)
  cat("Region coverage statistics:\n")
  cat("  - Total regions:", ncol(matrix_data), "\n")
  cat("  - Mean cells per region:", round(mean(region_coverage), 1), "\n")
  cat("  - Median cells per region:", median(region_coverage), "\n")
  cat("  - Min cells per region:", min(region_coverage), "\n")
  cat("  - Max cells per region:", max(region_coverage), "\n")
  
  # Show regions with sparse coverage
  sparse_regions <- names(region_coverage[region_coverage < 20])
  if(length(sparse_regions) > 0) {
    cat("  - Regions with <20 projecting cells (", length(sparse_regions), "):\n")
    cat("   ", paste(sparse_regions, collapse = ", "), "\n")
  }
  
  # Cell projection breadth (how many regions each cell projects to)
  cell_breadth <- rowSums(matrix_data > 0)
  cat("\nCell projection breadth statistics:\n")
  cat("  - Total cells:", nrow(matrix_data), "\n")
  cat("  - Mean regions per cell:", round(mean(cell_breadth), 1), "\n")
  cat("  - Median regions per cell:", median(cell_breadth), "\n")
  cat("  - Min regions per cell:", min(cell_breadth), "\n")
  cat("  - Max regions per cell:", max(cell_breadth), "\n")
  
  # Show cells with narrow projection patterns
  narrow_cells <- sum(cell_breadth < 3)
  broad_cells <- sum(cell_breadth > 20)
  cat("  - Cells projecting to <3 regions:", narrow_cells, "(", round(100*narrow_cells/nrow(matrix_data), 1), "%)\n")
  cat("  - Cells projecting to >20 regions:", broad_cells, "(", round(100*broad_cells/nrow(matrix_data), 1), "%)\n")
  
  # Distribution of projection values (non-zero)
  non_zero_values <- matrix_data[matrix_data > 0]
  
  # FIXED: Correct sparsity calculation
  total_entries <- nrow(matrix_data) * ncol(matrix_data)
  zero_entries <- sum(matrix_data == 0)
  sparsity_percent <- round(100 * zero_entries / total_entries, 1)
  
  cat("\nProjection value distribution (non-zero only):\n")
  cat("  - Mean projection value:", round(mean(non_zero_values), 4), "\n")
  cat("  - Median projection value:", round(median(non_zero_values), 4), "\n")
  cat("  - Sparsity (% zeros):", sparsity_percent, "%\n")  # FIXED
  
  # Create structured results for CSV export
  quality_summary <- data.frame(
    Matrix_Name = matrix_name,
    Total_Regions = ncol(matrix_data),
    Total_Cells = nrow(matrix_data),
    Mean_Cells_Per_Region = round(mean(region_coverage), 1),
    Median_Cells_Per_Region = median(region_coverage),
    Min_Cells_Per_Region = min(region_coverage),
    Max_Cells_Per_Region = max(region_coverage),
    Regions_With_Sparse_Coverage = length(sparse_regions),
    Mean_Regions_Per_Cell = round(mean(cell_breadth), 1),
    Median_Regions_Per_Cell = median(cell_breadth),
    Min_Regions_Per_Cell = min(cell_breadth),
    Max_Regions_Per_Cell = max(cell_breadth),
    Cells_Narrow_Projection_Count = narrow_cells,
    Cells_Narrow_Projection_Percent = round(100*narrow_cells/nrow(matrix_data), 1),
    Cells_Broad_Projection_Count = broad_cells,
    Cells_Broad_Projection_Percent = round(100*broad_cells/nrow(matrix_data), 1),
    Mean_Projection_Value = round(mean(non_zero_values), 4),
    Median_Projection_Value = round(median(non_zero_values), 4),
    Sparsity_Percent = sparsity_percent,
    Total_Entries = total_entries,
    Zero_Entries = zero_entries,
    NonZero_Entries = total_entries - zero_entries
  )
  
  # Save to CSV if requested
  if(save_csv) {
    # Clean matrix name for filename
    clean_name <- gsub("[^A-Za-z0-9_-]", "_", matrix_name)
    
    # Summary statistics
    summary_filename <- paste0("DataQuality_Summary_", clean_name, ".csv")
    write.csv(quality_summary, summary_filename, row.names = FALSE)
    cat("  - Summary saved to:", summary_filename, "\n")
    
    # Region coverage details
    region_details <- data.frame(
      Matrix_Name = matrix_name,
      Region = names(region_coverage),
      Cells_Projecting = as.numeric(region_coverage),
      Is_Sparse = region_coverage < 20
    )
    region_filename <- paste0("DataQuality_Regions_", clean_name, ".csv")
    write.csv(region_details, region_filename, row.names = FALSE)
    cat("  - Region details saved to:", region_filename, "\n")
  }
  
  return(list(
    summary = quality_summary,
    region_coverage = region_coverage,
    cell_breadth = cell_breadth,
    sparse_regions = sparse_regions,
    narrow_cells = narrow_cells,
    broad_cells = broad_cells,
    sparsity_percent = sparsity_percent
  ))
}

# Function to combine all assessments into one comparison file
save_combined_quality_assessment <- function(quality_list, filename = "DataQuality_All_Matrices_Comparison.csv") {
  combined_summary <- do.call(rbind, lapply(quality_list, function(x) x$summary))
  write.csv(combined_summary, filename, row.names = FALSE)
  cat("Combined quality assessment saved to:", filename, "\n")
  return(combined_summary)
}

# Assess data quality
quality_raw <- assess_data_quality(raw_matrix, "Raw counts")
quality_log <- assess_data_quality(log_matrix, "Log-transformed matrix")
quality_rownorm <- assess_data_quality(rownorm_matrix, "Row-normalized matrix")

# Extract hemisphere data using inclusive pattern to capture all variants
RH_cols <- grep("_RH", colnames(rownorm_matrix))
LH_cols <- grep("_LH", colnames(rownorm_matrix))

RH <- rownorm_matrix[, RH_cols]
LH <- rownorm_matrix[, LH_cols]

# Also assess hemisphere-specific quality
quality_RH <- assess_data_quality(RH, "RH hemisphere")
quality_LH <- assess_data_quality(LH, "LH hemisphere")

# Function to create balanced matrices for cross-hemisphere analysis
create_balanced_matrices <- function(RH, LH) {
  RH_base <- gsub("_RH(\\.\\d+)?$", "", colnames(RH))
  LH_base <- gsub("_LH(\\.\\d+)?$", "", colnames(LH))
  
  # Keep only regions present in both hemispheres
  common_bases <- intersect(RH_base, LH_base)
  
  # For each common base, select variant with highest total projections
  RH_keep <- c()
  LH_keep <- c()
  
  for (base in common_bases) {
    # Select variant with highest total projections
    RH_candidates <- which(RH_base == base)
    LH_candidates <- which(LH_base == base)
    
    RH_sums <- colSums(RH[, RH_candidates, drop = FALSE])
    LH_sums <- colSums(LH[, LH_candidates, drop = FALSE])
    
    RH_best <- RH_candidates[which.max(RH_sums)]
    LH_best <- LH_candidates[which.max(LH_sums)]
    
    RH_keep <- c(RH_keep, RH_best)
    LH_keep <- c(LH_keep, LH_best)
  }
  
  RH_balanced <- RH[, RH_keep]
  LH_balanced <- LH[, LH_keep]
  
  cat("Balanced matrices created - RH:", ncol(RH_balanced), "columns, LH:", ncol(LH_balanced), "columns\n")
  cat("Selected highest-projection variants for", length(common_bases), "common regions\n")
  
  return(list(RH = RH_balanced, LH = LH_balanced))
}

# Calculate cross-correlation function (NO FDR correction here)
calculate_cross_correlation <- function(hemisphere1, hemisphere2) {
  # Initialize matrices to store the correlation values and p-values
  correlation_matrix <- matrix(nrow = ncol(hemisphere1), ncol = ncol(hemisphere2), dimnames = list(colnames(hemisphere1), colnames(hemisphere2)))
  pvalue_matrix <- matrix(nrow = ncol(hemisphere1), ncol = ncol(hemisphere2), dimnames = list(colnames(hemisphere1), colnames(hemisphere2)))
  
  # Compute the correlation and p-value between each pair of corresponding regions
  for (i in 1:ncol(hemisphere1)) {
    for (j in 1:ncol(hemisphere2)) {
      # Remove rows with missing values and add a small constant to avoid division by zero
      h1_clean <- na.omit(hemisphere1[, i]) + 1e-10
      h2_clean <- na.omit(hemisphere2[, j]) + 1e-10
      
      test_result <- cor.test(h1_clean, h2_clean, method = "spearman")
      correlation_matrix[i, j] <- test_result$estimate
      pvalue_matrix[i, j] <- test_result$p.value
    }
  }
  
  # DON'T adjust p-values here - return raw p-values
  # The adjustment will be done globally across all analyses
  
  # Return the correlation matrix and p-value matrix
  list(correlation_matrix = correlation_matrix, pvalue_matrix = pvalue_matrix)
}

# Function for global FDR correction
apply_global_fdr <- function(pvalue_list) {
  # Combine all p-values from all analyses
  all_pvalues <- unlist(lapply(pvalue_list, as.vector))
  
  # Apply global FDR correction
  adjusted_pvalues_global <- p.adjust(all_pvalues, method = "BH")
  
  # Split back into original matrix structures
  start_idx <- 1
  adjusted_list <- list()
  
  for(i in 1:length(pvalue_list)) {
    n_values <- length(as.vector(pvalue_list[[i]]))
    end_idx <- start_idx + n_values - 1
    
    adjusted_matrix <- matrix(adjusted_pvalues_global[start_idx:end_idx], 
                              nrow = nrow(pvalue_list[[i]]), 
                              ncol = ncol(pvalue_list[[i]]))
    rownames(adjusted_matrix) <- rownames(pvalue_list[[i]])
    colnames(adjusted_matrix) <- colnames(pvalue_list[[i]])
    
    adjusted_list[[i]] <- adjusted_matrix
    start_idx <- end_idx + 1
  }
  
  return(adjusted_list)
}

# Generate heatmaps function with midline exclusion for all plots
generate_heatmaps <- function(correlation_matrix, pvalue_matrix, filename_prefix) {
  # Create a color palette function
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  
  # Define midline regions by name pattern (regions that cross midline)
  midline_patterns <- c("olf bulb", "orb ctx",  "AON", "ctx 3", "cc", "septum", "NAc", "BNST", "hippocampus", "thalamus", "hypothalamus", "midbrain", "hindbrain", "cerebellum", "medulla")
  
  # Get indices of midline regions in the current matrix
  row_names <- rownames(correlation_matrix)
  col_names <- colnames(correlation_matrix)
  
  midline_row_indices <- which(sapply(row_names, function(x) any(sapply(midline_patterns, function(pattern) grepl(pattern, x, ignore.case = TRUE)))))
  midline_col_indices <- which(sapply(col_names, function(x) any(sapply(midline_patterns, function(pattern) grepl(pattern, x, ignore.case = TRUE)))))
  
  # Function to apply midline exclusion to diagonal elements (for same-hemisphere comparisons)
  apply_midline_exclusion <- function(matrix_data) {
    if (nrow(matrix_data) == ncol(matrix_data)) {
      for (i in midline_row_indices) {
        if (i %in% midline_col_indices) {
          matrix_data[i, i] <- NA
        }
      }
    }
    return(matrix_data)
  }
  
  # Apply midline exclusion to correlation matrix for all plots
  correlation_matrix_midline <- apply_midline_exclusion(correlation_matrix)
  
  # Visualize the correlation matrix WITH midline exclusion
  heatmap.2(correlation_matrix_midline, Rowv = NA, Colv = NA, cexRow = 0.5, cexCol = 0.5, keysize = 1, col = my_palette(100), 
            trace = "none", margins = c(8, 8), dendrogram = "none", na.color = "lightgrey", main = "Correlation Matrix Heatmap (Midline Excluded)")
  dev.copy(pdf, paste0(filename_prefix, " Row-norm cross-correlation matched ROI pairs.pdf"), width = 14, height = 11)
  dev.off()
  
  # Visualize the p-value matrix
  heatmap.2(pvalue_matrix, Rowv = NA, Colv = NA, cexRow = 0.5, cexCol = 0.5, keysize = 1, col =  my_palette(100), 
            trace = "none", margins = c(8, 8), dendrogram = "none", main = "P-value Matrix Heatmap")
  dev.copy(pdf, paste0(filename_prefix, " Row-norm p-values matched ROI pairs.pdf"), width = 14, height = 11)
  dev.off()
  
  # Create a new matrix where all values corresponding to p-values below the threshold are set to NA
  threshold <- 0.05
  correlation_matrix_threshold_1 <- correlation_matrix_midline  # Start with midline-excluded matrix
  correlation_matrix_threshold_1[pvalue_matrix > threshold] <- NA
  
  heatmap.2(correlation_matrix_threshold_1, Rowv = NA, Colv = NA, cexRow = 0.5, cexCol = 0.5, keysize = 1, col =  my_palette(100), 
            trace = "none", margins = c(8, 8), dendrogram = "none", na.color = "lightgrey", main = "Thresholded Correlation Matrix Heatmap (Midline Excluded)")
  dev.copy(pdf, paste0(filename_prefix, " Row-norm cross-correlation matched ROI pairs with threshold.pdf"), width = 14, height = 11)
  dev.off()
  
  #Return the thresholded correlation matrix (with midline exclusion applied)
  return(correlation_matrix_threshold_1)
}

# Run analyses and collect raw results (no individual FDR correction):
# For cross-hemisphere analysis, create balanced matrices
balanced <- create_balanced_matrices(RH, LH)
result_RH_LH <- calculate_cross_correlation(balanced$RH, balanced$LH)

# Within-hemisphere analyses (use all available columns)
result_RH_RH <- calculate_cross_correlation(RH, RH)
result_LH_LH <- calculate_cross_correlation(LH, LH)

# Apply global FDR correction across all analyses
pvalue_list <- list(result_RH_LH$pvalue_matrix, result_RH_RH$pvalue_matrix, result_LH_LH$pvalue_matrix)
adjusted_pvalue_list <- apply_global_fdr(pvalue_list)

cat("Global FDR correction applied across", sum(sapply(pvalue_list, length)), "total tests\n")

# Generate heatmaps with globally corrected p-values
thresholded_correlation_matrix_RH_LH <- generate_heatmaps(result_RH_LH$correlation_matrix, adjusted_pvalue_list[[1]], "RH-LH")
thresholded_correlation_matrix_RH_RH <- generate_heatmaps(result_RH_RH$correlation_matrix, adjusted_pvalue_list[[2]], "RH-RH")
thresholded_correlation_matrix_LH_LH <- generate_heatmaps(result_LH_LH$correlation_matrix, adjusted_pvalue_list[[3]], "LH-LH")


##########################################################################################################################################################################################################
#################### clustering analyses for regional similarity #################
# Function to perform clustering analysis using projection patterns (not correlation matrices)
perform_clustering_analysis <- function(projection_matrix, correlation_matrix, analysis_name, k_clusters = 4) {
  cat("\n=== Clustering analysis for", analysis_name, "===\n")
  
  # Use the original projection data (regions as columns, cells as rows)
  # Transpose so regions are rows (what we want to cluster)
  region_data <- t(projection_matrix)
  
  # Remove regions with too many NAs/zeros if needed
  valid_regions <- rowSums(!is.na(region_data)) > ncol(region_data) * 0.1  # At least 10% non-NA
  region_data_clean <- region_data[valid_regions, ]
  
  cat("Clustering", nrow(region_data_clean), "regions based on projection patterns across", ncol(region_data_clean), "cells\n")
  
  # Compare different distance metrics
  distance_methods <- list(
    "Correlation" = function(x) as.dist(1 - cor(t(x), use = "pairwise.complete.obs")),
    "Cosine" = function(x) proxy::dist(x, method = "cosine"),
    "Euclidean" = function(x) dist(x)  # FIXED: Don't scale, just use raw distances between regions
  )
  
  # Store results for comparison
  clustering_results <- list()
  
  for (method_name in names(distance_methods)) {
    cat("Computing", method_name, "distance...\n")
    
    # Compute distance matrix
    dist_matrix <- distance_methods[[method_name]](region_data_clean)
    
    # Handle any remaining NAs
    if(any(is.na(dist_matrix))) {
      max_dist <- max(dist_matrix, na.rm = TRUE)
      dist_matrix[is.na(dist_matrix)] <- max_dist * 1.5
    }
    
    # Hierarchical clustering
    hc <- hclust(dist_matrix, method = "ward.D2")
    
    # For k-means, use TRANSPOSED data so we're clustering regions (not cells)
    k.max <- min(15, nrow(region_data_clean) - 1)
    wss <- sapply(1:k.max, function(k) {
      if(k < nrow(region_data_clean) && k > 0) {
        tryCatch({
          # FIXED: Use region_data_clean directly (regions as rows)
          kmeans(region_data_clean, centers = k, nstart = 50, iter.max = 100)$tot.withinss
        }, error = function(e) NA)
      } else NA
    })
    
    # Store results
    clustering_results[[method_name]] <- list(
      hclust = hc,
      distance = dist_matrix,
      wss = wss,
      clusters = cutree(hc, k = k_clusters)  
    )
  }
  
  # Create comprehensive comparison plots
  # First display on screen
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  
  # Plot dendrograms for each method
  for (method_name in names(distance_methods)) {
    hc <- clustering_results[[method_name]]$hclust
    clusters <- clustering_results[[method_name]]$clusters
    
    plot(hc, hang = -1, main = paste(method_name, "Distance -", analysis_name), 
         xlab = "", ylab = "Height", cex.main = 0.9, cex.axis = 0.7, cex = 0.6)  
    rect.hclust(hc, k = k_clusters, border = 2:(k_clusters+1))
  }
  
  # Plot elbow plots for each method
  for (method_name in names(distance_methods)) {
    wss <- clustering_results[[method_name]]$wss
    valid_k <- which(!is.na(wss))
    
    if(length(valid_k) > 1) {
      plot(valid_k, wss[valid_k], type="b", pch = 19, frame = FALSE, 
           xlab="Number of clusters K", ylab="WSS",
           main = paste("Elbow -", method_name), cex.main = 0.9, cex.axis = 0.7)
    }
  }
  
  # Reset to single plot layout
  par(mfrow = c(1, 1))
  
  # Save detailed comparison to PDF
  pdf(paste0("Cluster_analysis_", gsub("-", "_", analysis_name), "_comparison.pdf"), width = 18, height = 12)
  par(mfrow = c(2, 3), mar = c(6, 4, 3, 1))  # FIXED: More space for labels
  
  # Plot dendrograms
  for (method_name in names(distance_methods)) {
    hc <- clustering_results[[method_name]]$hclust
    clusters <- clustering_results[[method_name]]$clusters
    
    plot(hc, hang = -1, main = paste(method_name, "Distance -", analysis_name), 
         xlab = "", ylab = "Height", cex.main = 1.1, cex.axis = 0.8, cex = 0.7)  
    rect.hclust(hc, k = k_clusters, border = 2:(k_clusters+1))
  }
  
  # Plot elbow plots
  for (method_name in names(distance_methods)) {
    wss <- clustering_results[[method_name]]$wss
    valid_k <- which(!is.na(wss))
    
    if(length(valid_k) > 1) {
      plot(valid_k, wss[valid_k], type="b", pch = 19, frame = FALSE, 
           xlab="Number of clusters K", ylab="Within-cluster sum of squares",
           main = paste("Elbow Plot -", method_name), cex.main = 1.1, cex.axis = 0.8)
    }
  }
  
  dev.off()
  
  # Reset par settings
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
  
  # Print cluster comparison - FIXED: Now all vectors have same length
  cat("\nCluster assignments comparison (k =", k_clusters, "):\n")
  cluster_comparison <- data.frame(
    Region = rownames(region_data_clean),
    Correlation = clustering_results[["Correlation"]]$clusters,
    Cosine = clustering_results[["Cosine"]]$clusters,
    Euclidean = clustering_results[["Euclidean"]]$clusters
  )
  print(cluster_comparison)
  
  # Manual Adjusted Rand Index calculation (no mclust needed)
  adjusted_rand_index <- function(x, y) {
    # Contingency table
    tab <- table(x, y)
    n <- sum(tab)
    
    # Row and column sums
    a <- rowSums(tab)
    b <- colSums(tab)
    
    # Calculate components
    sum_comb_a <- sum(choose(a, 2))
    sum_comb_b <- sum(choose(b, 2))
    sum_comb_tab <- sum(choose(tab, 2))
    
    # Expected index
    expected_index <- sum_comb_a * sum_comb_b / choose(n, 2)
    
    # Max index
    max_index <- (sum_comb_a + sum_comb_b) / 2
    
    # Adjusted Rand Index
    if (max_index - expected_index == 0) {
      return(0)
    } else {
      return((sum_comb_tab - expected_index) / (max_index - expected_index))
    }
  }
  
  # Calculate agreement between methods
  cat("\nCluster agreement between methods:\n")
  methods <- names(distance_methods)
  for(i in 1:(length(methods)-1)) {
    for(j in (i+1):length(methods)) {
      # Calculate adjusted rand index for cluster agreement
      ari <- adjusted_rand_index(clustering_results[[methods[i]]]$clusters, 
                                 clustering_results[[methods[j]]]$clusters)
      cat(paste(methods[i], "vs", methods[j], "- ARI:", round(ari, 3), "\n"))
    }
  }
  
  # Return comprehensive results (default to correlation distance for main results)
  return(list(
    clusters = clustering_results[["Correlation"]]$clusters,
    hclust = clustering_results[["Correlation"]]$hclust,
    wss = clustering_results[["Correlation"]]$wss,
    all_methods = clustering_results,
    cluster_comparison = cluster_comparison
  ))
}

# Modified function calls that pass both projection and correlation matrices:
# Clustering analysis for all three analyses
clustering_RH_RH <- perform_clustering_analysis(RH, thresholded_correlation_matrix_RH_RH, "RH-RH", k_clusters = 4)
clustering_LH_LH <- perform_clustering_analysis(LH, thresholded_correlation_matrix_LH_LH, "LH-LH", k_clusters = 4)

# For RH-LH, use the balanced matrices
balanced <- create_balanced_matrices(RH, LH)
clustering_RH_LH <- perform_clustering_analysis(balanced$RH, thresholded_correlation_matrix_RH_LH, "RH-LH", k_clusters = 4)

##########################################################################################################################################################################################################
#################### correlational and clustering analyses for regional similarity with ipsi-contra matrix#################
# Function to create ipsi-contra matrices using your approach
create_ipsi_contra_matrices_clean <- function(projection_matrix, inRH_lookup) {
  
  cat("\n=== Creating Ipsi/Contra matrices ===\n")
  
  df <- as.data.frame(projection_matrix)
  
  # --- NEW: safer way to attach inRH using what the loader provides ---
  # If the loader has already aligned inRH_lookup to these rows, use it directly.
  if (all(rownames(df) == inRH_lookup$row_id)) {
    df$inRH <- inRH_lookup$inRH
  } else {
    # Fallback to old behavior: match by row names to row_id
    idx <- match(rownames(df), inRH_lookup$row_id)
    if (any(is.na(idx))) {
      stop("Some rows of projection_matrix not found in inRH_lookup$row_id – check alignment.")
    }
    df$inRH <- inRH_lookup$inRH[idx]
  }
  # --- end new bit ---
  
  cat("Initial data:\n")
  cat("- Total cells:", nrow(df), "\n")
  cat("- LH soma cells:", sum(df$inRH == 0, na.rm = TRUE), "\n")
  cat("- RH soma cells:", sum(df$inRH == 1, na.rm = TRUE), "\n")
  cat("- Total projection columns:", ncol(df) - 1, "\n")
  
  # Split by soma location
  L_soma <- df[df$inRH == 0, , drop = FALSE]
  R_soma <- df[df$inRH == 1, , drop = FALSE]
  
  # Relabel for L_soma and R_soma
  # For L_soma (LH soma: ipsi = LH, contra = RH)
  colnames(L_soma) <- gsub("_LH(\\.[0-9]+)?$", "-ipsi\\1",   colnames(L_soma), perl = TRUE)
  colnames(L_soma) <- gsub("_RH(\\.[0-9]+)?$", "-contra\\1", colnames(L_soma), perl = TRUE)
  
  # For R_soma (RH soma: ipsi = RH, contra = LH)
  colnames(R_soma) <- gsub("_RH(\\.[0-9]+)?$", "-ipsi\\1",   colnames(R_soma), perl = TRUE)
  colnames(R_soma) <- gsub("_LH(\\.[0-9]+)?$", "-contra\\1", colnames(R_soma), perl = TRUE)
  
  # Remove inRH column
  L_soma$inRH <- NULL
  R_soma$inRH <- NULL
  
  cat("After relabeling:\n")
  cat("- L_soma columns:", ncol(L_soma), "\n")
  cat("- R_soma columns:", ncol(R_soma), "\n")
  
  # Align columns, only keep shared ROIs
  common_cols <- intersect(colnames(L_soma), colnames(R_soma))
  L_soma <- L_soma[, common_cols, drop = FALSE]
  R_soma <- R_soma[, common_cols, drop = FALSE]
  
  cat("After matching:\n")
  cat("- Common projection columns:", length(common_cols), "\n")
  cat("- Regions dropped from L_soma:", ncol(df) - 1 - length(common_cols), "\n")
  cat("- Regions dropped from R_soma:", ncol(df) - 1 - length(common_cols), "\n")
  
  # Combine all cells
  ipsi_contra_combined <- rbind(L_soma, R_soma)
  
  # Split into separate ipsi and contra matrices
  ipsi_cols   <- grep("-ipsi",   colnames(ipsi_contra_combined))
  contra_cols <- grep("-contra", colnames(ipsi_contra_combined))
  
  ipsi_matrix   <- ipsi_contra_combined[, ipsi_cols,   drop = FALSE]
  contra_matrix <- ipsi_contra_combined[, contra_cols, drop = FALSE]
  
  # Clean up column names (remove -ipsi/-contra suffixes for easier analysis)
  colnames(ipsi_matrix)   <- gsub("-ipsi.*",   "", colnames(ipsi_matrix))
  colnames(contra_matrix) <- gsub("-contra.*", "", colnames(contra_matrix))
  
  cat("Final matrices:\n")
  cat("- Ipsi matrix:",   ncol(ipsi_matrix),   "regions x", nrow(ipsi_matrix),   "cells\n")
  cat("- Contra matrix:", ncol(contra_matrix), "regions x", nrow(contra_matrix), "cells\n")
  cat("- Total cells used:", nrow(ipsi_matrix), "(doubled from", nrow(projection_matrix), ")\n")
  
  return(list(
    ipsi   = ipsi_matrix,
    contra = contra_matrix,
    combined = ipsi_contra_combined,
    summary = list(
      n_cells_total     = nrow(ipsi_matrix),
      n_regions_matched = ncol(ipsi_matrix),
      n_lh_soma         = sum(df$inRH == 0, na.rm = TRUE),
      n_rh_soma         = sum(df$inRH == 1, na.rm = TRUE)
    )
  ))
}


# Improved heatmap function with clear axis labels for laterality analysis
generate_heatmaps_laterality <- function(correlation_matrix, pvalue_matrix, analysis_type) {
  # Create a color palette function
  my_palette <- colorRampPalette(c("blue", "white", "red"))
  
  # Define axis labels based on analysis type
  if(analysis_type == "Ipsi-Ipsi") {
    x_label <- "Ipsilateral Regions"
    y_label <- "Ipsilateral Regions"
    main_title <- "Ipsi-Ipsi Correlations: How ipsilateral projections co-vary"
    filename_prefix <- "Ipsi-Ipsi"
  } else if(analysis_type == "Contra-Contra") {
    x_label <- "Contralateral Regions"
    y_label <- "Contralateral Regions"
    main_title <- "Contra-Contra Correlations: How contralateral projections co-vary"
    filename_prefix <- "Contra-Contra"
  } else if(analysis_type == "Ipsi-Contra") {
    x_label <- "Ipsilateral Regions"
    y_label <- "Contralateral Regions"
    main_title <- "Ipsi-Contra Correlations: Ipsi vs Contra to same regions"
    filename_prefix <- "Ipsi-Contra"
  }
  
  # Apply midline exclusion (same as before)
  midline_patterns <- c("olf bulb", "orb ctx", "AON", "ctx 3", "cc", "septum", "NAc", "BNST", "hippocampus", "thalamus", "hypothalamus", "midbrain", "hindbrain", "cerebellum", "medulla")
  
  row_names <- rownames(correlation_matrix)
  col_names <- colnames(correlation_matrix)
  
  midline_row_indices <- which(sapply(row_names, function(x) any(sapply(midline_patterns, function(pattern) grepl(pattern, x, ignore.case = TRUE)))))
  midline_col_indices <- which(sapply(col_names, function(x) any(sapply(midline_patterns, function(pattern) grepl(pattern, x, ignore.case = TRUE)))))
  
  apply_midline_exclusion <- function(matrix_data) {
    if (nrow(matrix_data) == ncol(matrix_data)) {
      for (i in midline_row_indices) {
        if (i %in% midline_col_indices) {
          matrix_data[i, i] <- NA
        }
      }
    }
    return(matrix_data)
  }
  
  correlation_matrix_midline <- apply_midline_exclusion(correlation_matrix)
  
  # Generate heatmaps with clear labels
  # 1. Basic correlation matrix
  heatmap.2(correlation_matrix_midline, Rowv = NA, Colv = NA, 
            cexRow = 0.5, cexCol = 0.5, keysize = 1, col = my_palette(100), 
            trace = "none", margins = c(10, 10), dendrogram = "none", 
            na.color = "lightgrey", 
            main = main_title,
            xlab = x_label, ylab = y_label)
  dev.copy(pdf, paste0(filename_prefix, " Row-norm cross-correlation.pdf"), width = 14, height = 11)
  dev.off()
  
  # 2. P-value matrix
  heatmap.2(pvalue_matrix, Rowv = NA, Colv = NA, 
            cexRow = 0.5, cexCol = 0.5, keysize = 1, col = my_palette(100), 
            trace = "none", margins = c(10, 10), dendrogram = "none",
            main = paste("P-values:", analysis_type),
            xlab = x_label, ylab = y_label)
  dev.copy(pdf, paste0(filename_prefix, " Row-norm p-values.pdf"), width = 14, height = 11)
  dev.off()
  
  # 3. Thresholded matrix
  threshold <- 0.05
  correlation_matrix_threshold <- correlation_matrix_midline
  correlation_matrix_threshold[pvalue_matrix > threshold] <- NA
  
  heatmap.2(correlation_matrix_threshold, Rowv = NA, Colv = NA, 
            cexRow = 0.5, cexCol = 0.5, keysize = 1, col = my_palette(100), 
            trace = "none", margins = c(10, 10), dendrogram = "none", 
            na.color = "lightgrey",
            main = paste("Significant Correlations:", analysis_type, "(p < 0.05)"),
            xlab = x_label, ylab = y_label)
  dev.copy(pdf, paste0(filename_prefix, " Row-norm significant-correlations.pdf"), width = 14, height = 11)
  dev.off()
  
  return(correlation_matrix_threshold)
}
# Create ipsi/contra matrices using row-normalized data
ipsi_contra <- create_ipsi_contra_matrices_clean(rownorm_matrix, dat$inRH_lookup)

# Assess data quality for ipsi/contra (should be much better with doubled sample size)
quality_ipsi <- assess_data_quality(ipsi_contra$ipsi, "Ipsilateral projections")
quality_contra <- assess_data_quality(ipsi_contra$contra, "Contralateral projections")

# Combine all assessments at the end
all_quality <- list(quality_ipsi, quality_contra, quality_RH, quality_LH)
combined_quality <- save_combined_quality_assessment(all_quality)

# Run correlation analyses with raw p-values first
result_ipsi_ipsi <- calculate_cross_correlation(ipsi_contra$ipsi, ipsi_contra$ipsi)
result_contra_contra <- calculate_cross_correlation(ipsi_contra$contra, ipsi_contra$contra)
result_ipsi_contra <- calculate_cross_correlation(ipsi_contra$ipsi, ipsi_contra$contra)

cat("Correlation analyses completed for:\n")
cat("- Ipsi-Ipsi:", nrow(result_ipsi_ipsi$correlation_matrix), "x", ncol(result_ipsi_ipsi$correlation_matrix), "correlations\n")
cat("- Contra-Contra:", nrow(result_contra_contra$correlation_matrix), "x", ncol(result_contra_contra$correlation_matrix), "correlations\n")
cat("- Ipsi-Contra:", nrow(result_ipsi_contra$correlation_matrix), "x", ncol(result_ipsi_contra$correlation_matrix), "correlations\n")

# Apply global FDR correction across all laterality analyses
pvalue_list_laterality <- list(
  result_ipsi_ipsi$pvalue_matrix, 
  result_contra_contra$pvalue_matrix, 
  result_ipsi_contra$pvalue_matrix
)
adjusted_pvalue_list_laterality <- apply_global_fdr(pvalue_list_laterality)

total_tests_laterality <- sum(sapply(pvalue_list_laterality, length))
cat("Global FDR correction applied across", total_tests_laterality, "laterality tests\n")

thresholded_correlation_matrix_ipsi_ipsi <- generate_heatmaps_laterality(
  result_ipsi_ipsi$correlation_matrix, 
  adjusted_pvalue_list_laterality[[1]], 
  "Ipsi-Ipsi"
)

thresholded_correlation_matrix_contra_contra <- generate_heatmaps_laterality(
  result_contra_contra$correlation_matrix, 
  adjusted_pvalue_list_laterality[[2]], 
  "Contra-Contra"
)

thresholded_correlation_matrix_ipsi_contra <- generate_heatmaps_laterality(
  result_ipsi_contra$correlation_matrix, 
  adjusted_pvalue_list_laterality[[3]], 
  "Ipsi-Contra"
)


# Run clustering analyses on projection patterns
cat("\n=== Running laterality clustering analyses ===\n")
clustering_ipsi_ipsi <- perform_clustering_analysis(
  ipsi_contra$ipsi, 
  thresholded_correlation_matrix_ipsi_ipsi, 
  "Ipsi-Ipsi", 
  k_clusters = 4
)

clustering_contra_contra <- perform_clustering_analysis(
  ipsi_contra$contra, 
  thresholded_correlation_matrix_contra_contra, 
  "Contra-Contra", 
  k_clusters = 4
)

