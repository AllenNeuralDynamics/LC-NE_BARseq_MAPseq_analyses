####### Prepare data for analyses following MAT-RDS conversion by performing cpm normalization and thresholding #######
# script order 02
# barseq-r4 environment to run


############################################################################################################################################################################################################
#specify working directory
BARSEQ_DIR <- '/data/BARseq_MATtoRDSfiles_brain3_brain4/BARseq_780345/'

#function to load and process sce file for normalization and QC cut offs
load_barseq = function(filename="combined_neurons_clust_CCFv2_uid.rds", slice = NULL, normalization_factor = 10, min_genes = 5, min_counts = 20) {
  sce = readRDS(file.path(BARSEQ_DIR, filename))
  sce = sce[, colSums(counts(sce)) >= min_counts & colSums(counts(sce)>0) >= min_genes]
  cpm(sce) = convert_to_cpm(counts(sce), normalization_factor)
  if (!is.null(slice)) {
    sce = sce[, sce$slice %in% slice]
  }
  return(sce)
}

#function to convert sparse matrix of counts into cpm
convert_to_cpm = function(M, total_counts = 1000000) {
  normalization_factor = Matrix::colSums(M) / total_counts
  if (is(M, "dgCMatrix")) {
    M@x = M@x / rep.int(normalization_factor, diff(M@p))
    return(M)
  } else {
    return(scale(M, center = FALSE, scale = normalization_factor))
  }
}

############################################################################################################################################################################################################
#load sce file, convert to cpm using SingleCellExperiment package, then log normalize and save data
# All cells
all = load_barseq(filename="combined_neurons_clust_CCFv2_uid.rds")
colData(all)
assayNames(all)
logcounts(all) = log1p(cpm(all))/log(2)
assayNames(all)
# Check if the directory exists, create it if not, then save the RDS file
dir_path <- "/scratch/BARseq_780345"
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}
saveRDS(all, file.path(dir_path, "combined_neurons_clust_CCFv2_uid_cpm_log.rds"))

############################################################################################################################################################################################################
rm(list = ls()[sapply(mget(ls(), .GlobalEnv), function(x) !is.function(x))])
#specify working directory
BARSEQ_DIR <- '/scratch/BARseq_780345/'

# Normalize and save LC cluster object
LC = load_barseq(filename="LCcluster_neurons_CCFv2_uid.rds")
colData(LC)
logcounts(LC) = log1p(cpm(LC))/log(2)
saveRDS(LC, "LCcluster_neurons_CCFv2_uid_cpm_log.rds")

############################################################################################################################################################################################################
rm(list = ls()[sapply(mget(ls(), .GlobalEnv), function(x) !is.function(x))])
#specify working directory
BARSEQ_DIR <- '/scratch/BARseq_780345/'

# Load and Normalize data for  LC-NE cells specifically
LCNE = load_barseq(filename="LCNE_cluster_neurons_CCFv2_uid.rds")
logcounts(LCNE) = log1p(cpm(LCNE))/log(2)
saveRDS(LCNE, "LCNE_cluster_neurons_CCFv2_uid_cpm_log.rds")

############################################################################################################################################################################################################
rm(list = ls()[sapply(mget(ls(), .GlobalEnv), function(x) !is.function(x))])
#specify working directory
BARSEQ_DIR <- '/scratch/BARseq_780345/'

# Load and Normalize cleaned up LC-NE cells
LCNE = load_barseq(filename="LCNE_clusters_filtered.rds")
logcounts(LCNE) = log1p(cpm(LCNE))/log(2)
saveRDS(LCNE, "LCNE_clusters_filtered_cpm_log.rds")
