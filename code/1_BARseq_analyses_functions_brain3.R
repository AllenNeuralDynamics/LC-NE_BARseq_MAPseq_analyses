################# BARseq Analysis Functions #################

# Set up paths
BARSEQ_INPUT_DIR <- '/data/barseq_780345_2025-02-24_12-00-00_processed_MAT2RDS_2026-05-28_17-49-02/BARseq/'  # read-only
BARSEQ_OUTPUT_DIR <- '/results/BARseq_780345/'  # writable
BARSEQ_CLUSTERING_DIR <- '/code/cached_clustering/BARseq_780345/'  # pre-computed clustering, read-only — see clustering_freeze.md for provenance

# Ensure output directory exists
if (!dir.exists(BARSEQ_OUTPUT_DIR)) {
  dir.create(BARSEQ_OUTPUT_DIR, recursive = TRUE)
}

################## Load and process sce file for normalization and QC cutoffs ##################
load_barseq <- function(filename = "combined_neurons_clust_CCFv2_uid.rds", slice = NULL, normalization_factor = 10, min_genes = 5, min_counts = 20, from_output = FALSE) {
  base_dir <- if(from_output) BARSEQ_OUTPUT_DIR else BARSEQ_INPUT_DIR
  sce <- readRDS(file.path(base_dir, filename))
  sce <- sce[, colSums(counts(sce)) >= min_counts & colSums(counts(sce) > 0) >= min_genes]
  cpm(sce) <- convert_to_cpm(counts(sce), normalization_factor)
  if (!is.null(slice)) {
    sce <- sce[, sce$slice %in% slice]
  }
  return(sce)
}

################## Convert sparse matrix to CPM ##################
convert_to_cpm <- function(M, total_counts = 1000000) {
  normalization_factor <- Matrix::colSums(M) / total_counts
  if (is(M, "dgCMatrix")) {
    M@x <- M@x / rep.int(normalization_factor, diff(M@p))
    return(M)
  } else {
    return(scale(M, center = FALSE, scale = normalization_factor))
  }
}

################## Perform BARseq clustering ##################
analyze_barseq <- function(barseq, output_name, n_pca = 30, k_umap = 100, k_cluster = 50, seed = 42) {
  output_dir <- file.path("analysis", make.names(output_name))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  barseq <- scater::runPCA(barseq, ncomponents = n_pca, exprs_values = "logcounts")
  barseq <- scater::runUMAP(barseq, dimred = "PCA", n_neighbors = k_umap)
  write_csv(as.data.frame(reducedDim(barseq, "UMAP")), file.path(output_dir, "umap.csv"))
  
  clusters <- cluster_data(barseq, k_cluster, "jaccard", seed=seed)
  write_csv(clusters, file.path(output_dir, "cluster.csv"))
  
  # create default cluster annotation
  labels <- as.factor(clusters$label)
  cluster_annotation <- tibble(cluster_id = levels(labels)) %>%
    mutate(cluster_name = paste0(output_name, "_", cluster_id))
  write_csv(cluster_annotation, file.path(output_dir, "cluster_annotation.csv"))
  return(list(barseq, clusters))
}

################## Cluster data using SNN graph ##################
cluster_data <- function(barseq, k = 10, type = "rank", seed=42) {
  # type -> "rank" (preservation of neighbor rank),
  # "number" (number of shared neighbors), "jaccard"
  snn <- scran::buildSNNGraph(barseq, k = k, use.dimred = "PCA", type = type,
                              BNPARAM = BiocNeighbors::AnnoyParam(),
                              BPPARAM = BiocParallel::MulticoreParam(12, RNGseed=seed))
  label <- igraph::cluster_louvain(snn)$membership
  result <- tibble(sample = colnames(barseq), label)
  return(result)
}

################## Get color palette for clusters ##################
get_cluster_colors <- function(n_clusters) {
  base_palette <- c("#dcbeff", "#aaffc3", "#fabed4", "#e6194B", "#ffe119", 
                    "#469990", "#fffac8", "#f032e6", "#ffd8b1", "#bfef45", 
                    "#3cb44b", "#42d4f4", "#911eb4", "#4363d8", "#f58231",
                    "#800000")
  
  if (n_clusters <= length(base_palette)) {
    return(base_palette[1:n_clusters])
  } else {
    return(c(base_palette, rainbow(n_clusters - length(base_palette))))
  }
}

################## Calculate 3D spatial coherence ##################
#Computes kNN-based cluster purity: proportion of same-cluster neighbors.
# Distance weighted spatial coherence score calculation
calculate_spatial_coherence_3D <- function(sce, cluster_col = "louvain_cluster", k = 10, slice_thickness = 20) {
  metadata <- as.data.frame(colData(sce))
  
  metadata <- metadata %>%
    mutate(
      CCF_ML = as.numeric(CCF_ML) * 25,
      CCF_DV = as.numeric(CCF_DV) * 25,
      slice = as.numeric(slice),
      depth = slice * slice_thickness
    )
  
  coords_3d <- as.matrix(metadata[, c("CCF_ML", "CCF_DV", "depth")])
  clusters <- metadata[[cluster_col]]
  
  knn_result <- FNN::get.knn(data = coords_3d, k = k)
  neighbor_indices <- knn_result$nn.index
  neighbor_distances <- knn_result$nn.dist
  
  epsilon <- 1e-6
  
  coherence_scores <- sapply(seq_len(nrow(coords_3d)), function(i) {
    own_cluster <- clusters[i]
    neighbors <- neighbor_indices[i, ]
    distances <- neighbor_distances[i, ]
    weights <- 1 / (distances + epsilon)
    
    same_cluster <- as.numeric(clusters[neighbors] == own_cluster)
    weighted_same <- sum(weights * same_cluster)
    total_weight <- sum(weights)
    
    weighted_same / total_weight
  })
  
  colData(sce)$spatial_coherence_3D_weighted <- coherence_scores
  return(sce)
}

################## Calculate 3D spatial density ##################
# Calculates average kNN distance in 3D space, and inverts this distance → spatial density score (larger = denser).
calculate_spatial_density <- function(sce, k = 10, slice_thickness = 20) {
  metadata <- as.data.frame(colData(sce)) %>%
    mutate(
      CCF_ML = as.numeric(CCF_ML) * 25,
      CCF_DV = as.numeric(CCF_DV) * 25,
      slice = as.numeric(slice),
      depth = slice * slice_thickness
    )
  
  coords_3d <- as.matrix(metadata[, c("CCF_ML", "CCF_DV", "depth")])
  
  knn_result <- FNN::get.knn(data = coords_3d, k = k)
  mean_distances <- rowMeans(knn_result$nn.dist)
  
  density_score <- 1 / (mean_distances + 1e-6)
  
  colData(sce)$spatial_density_3D <- density_score
  return(sce)
}

################## Clean up work space between processing steps ##################
clear_objects_except_functions <- function() {
  obj_list <- ls(envir = .GlobalEnv)

  # Keep functions and the clearing function itself
  to_keep <- obj_list[sapply(obj_list, function(x) {
    obj <- get(x, envir = .GlobalEnv)
    is.function(obj)
  })]

  # Also keep directory paths and the cache-control flag
  to_keep <- c(to_keep, "BARSEQ_INPUT_DIR", "BARSEQ_OUTPUT_DIR", "BARSEQ_CLUSTERING_DIR", "fig_dir", "RECOMPUTE_CLUSTERING")

  to_remove <- setdiff(obj_list, to_keep)

  if(length(to_remove) > 0) {
    cat("Removing:", paste(to_remove, collapse = ", "), "\n")
    rm(list = to_remove, envir = .GlobalEnv)
  }
  
  # Force garbage collection
  gc()
}
