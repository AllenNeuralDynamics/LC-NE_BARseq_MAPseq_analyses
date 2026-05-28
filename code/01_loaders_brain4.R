MAPSEQ_INPUT_DIR <- '/data/mapseq_780346_2025-07-23_12-00-00/MAPseq/M305_20251030_USEthis/'

# Function to pre-process and load the data in a format directly usable by downstream analyses
load_data <-function() {

  # Set working directory
  setwd('/results/BARseq_780346/')
  
  # Read in raw projection matrix
  projection_matrix <- readr::read_csv("./MapSeq_matched_projections_1_mismatch.csv", col_names = TRUE)
  
  # Read and clean black list (skip header, remove decimal suffix)
  black_list <- readr::read_csv("./duplicate_blacklist_informed_final.csv", col_names = FALSE, skip = 1) %>%
    mutate(X1 = str_remove(X1, "\\..*$"))  # Remove from "." to end
  
  # Remove blacklisted cells
  rows_before <- nrow(projection_matrix)
  projection_matrix <- projection_matrix %>%
    filter(!CellID %in% black_list$X1)
  duplicates_removed <- rows_before - nrow(projection_matrix)
  
  # Report results
  cat("Duplicate handling summary:\n")
  cat("  Original rows:", rows_before, "\n")
  cat("  Duplicates removed:", duplicates_removed, "\n")
  cat("  Final rows:", nrow(projection_matrix), "\n")
  
  if(duplicates_removed > 0) {
    cat("✓ Removed", duplicates_removed, "duplicate CellIDs (via blacklist)\n")
  } else {
    cat("✓ No duplicates found\n")
  }
  
  # Load cell metadata
  LCNEneurons <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
  
  # Extract complete metadata
  cell_metadata <- as.data.frame(colData(LCNEneurons))
  
  # Load sample information
  proj_index <- read_tsv(file.path(MAPSEQ_INPUT_DIR, "M305sampleinfo.tsv"))
  head(proj_index)
  sampleinfo <- read_excel("/data/barseq_780346_2025-06-13_12-00-00_processed_MAT2RDS_2026-05-28_17-50-51/MAPseq/sampleinfo_780346.xlsx", sheet = "Sample information", skip = 1,range = "A2:J122")
  head(sampleinfo)
  colnames(sampleinfo) <- c("usertube", "ourtube", "samplename", "siteinfo", "QC_qPCR", "rtprimer", "brain", "hemisphere", "ROI", "notes")

  # Need to create a vector to index samples based on slide they came from to circumvent erroneous indexing due lost samples
  # Keep only the 117 "target" samples
  design_df <- sampleinfo %>%
    filter(siteinfo == "target")
  nrow(design_df)
  slide_counts <- c(2, 2, 6, 6, 14, 14, 14, 14, 16, 12, 6, 6, 2, 3)
  stopifnot(sum(slide_counts) == nrow(design_df))
  slide_vector <- rep(seq_along(slide_counts), times = slide_counts)
  length(slide_vector)  # 117
  design_df$slide <- slide_vector
  proj_index <- proj_index %>%
    left_join(
      design_df %>% select(rtprimer, ROI, hemisphere, slide),
      by = "rtprimer"
    )
  
  # Format the data by putting it together
  # Get the BC columns and match them to brain region names and samplename
  bc_cols <- grep("^BC\\d+$", colnames(projection_matrix), value = TRUE)
  rtprimer_numbers <- as.integer(sub("BC", "", bc_cols))
  # For each BC column, find the corresponding region and samplename
  hemi   <- proj_index$hemisphere[match(rtprimer_numbers, proj_index$rtprimer)]
  roi    <- proj_index$ROI[match(rtprimer_numbers, proj_index$rtprimer)]
  slide  <- proj_index$slide[match(rtprimer_numbers, proj_index$rtprimer)]
  # Combine region and samplename for new column names
  region_samplename <- paste0(roi, "_", hemi, ".", slide)
  # Rename BC columns to region_samplename
  colnames(projection_matrix)[match(bc_cols, colnames(projection_matrix))] <- region_samplename
  head(projection_matrix)
  
  # Merge cell metadata
  cell_metadata_subset <- cell_metadata %>%
    select(uid, slice, barcode, CCF_DV, CCF_ML, CCF_AP, louvain_cluster)
  df <- projection_matrix %>%
    left_join(cell_metadata_subset, by = c("CellID" = "uid")) %>%
    as.data.frame()
  
  # Add hemisphere classification (ML total is 456, midline around 228)
  df$inRH <- ifelse(df$CCF_ML > 228, 1, 0)
  df$CellID <- as.character(df$CellID)
  df$louvain_cluster <- as.character(df$louvain_cluster)
  df$row_id <- paste(df$CellID, df$louvain_cluster, sep=".")
  rownames(df) <- make.unique(df$row_id)
  
  # Create a lookup table for row_id and inRH
  inRH_lookup <- data.frame(
    row_id = rownames(df),
    inRH = df$inRH,
    stringsAsFactors = FALSE
  )
  
  # Extract ROI/hemisphere info
  col_names <- colnames(df)
  extract_info <- function(name) {
    hemi <- ifelse(str_detect(name, "_RH"), "RH",
                   ifelse(str_detect(name, "_LH"), "LH",
                          ifelse(str_detect(name, "_SP"), "SP", NA)))
    base_roi <- str_replace(name, "(_RH|_LH|_SP).*", "")
    return(data.frame(colname = name, base_roi = base_roi, hemisphere = hemi, stringsAsFactors = FALSE))
  }
  
  meta <- do.call(rbind, lapply(col_names, extract_info))
  
  # Extract projection columns only
  proj_cols <- meta$colname[!is.na(meta$hemisphere)]
  proj_matrix_raw <- df[, proj_cols, drop = FALSE]
  
  # Store metadata columns separately
  metadata_cols <- df[, !(colnames(df) %in% proj_cols), drop = FALSE]
  
  # Remove rows with all zero projections
  non_zero_rows <- rowSums(proj_matrix_raw, na.rm = TRUE) > 0
  proj_matrix_raw <- proj_matrix_raw[non_zero_rows, ]
  metadata_final <- metadata_cols[non_zero_rows, ]
  
  # Update inRH_lookup to match filtered data
  inRH_lookup_filtered <- inRH_lookup[non_zero_rows, ]
  
  # Order columns: RH, LH, SP
  ordered_cols <- c(
    meta$colname[meta$hemisphere == "RH"],
    meta$colname[meta$hemisphere == "LH"],
    meta$colname[meta$hemisphere == "SP"]
  )
  ordered_cols <- ordered_cols[!is.na(ordered_cols) & ordered_cols %in% colnames(proj_matrix_raw)]
  proj_matrix_raw <- proj_matrix_raw[, ordered_cols, drop = FALSE]
  
  
  # Cluster rows for ordering
  hc <- hclust(dist(proj_matrix_raw))
  row_order <- hc$order
  proj_matrix_raw <- proj_matrix_raw[row_order, ]
  metadata_final <- metadata_final[row_order, ]
  inRH_lookup_filtered <- inRH_lookup_filtered[row_order, ]
  
  # Log-transform all projection columns
  proj_matrix_log <- log10(1 + proj_matrix_raw * 100)
  proj_matrix_log <- as.matrix(proj_matrix_log)
  proj_matrix_log[is.na(proj_matrix_log) | is.nan(proj_matrix_log) | is.infinite(proj_matrix_log)] <- 0
  
  # Row-normalize all projection columns
  proj_matrix_rownorm <- proj_matrix_raw / rowSums(proj_matrix_raw)
  proj_matrix_rownorm[is.na(proj_matrix_rownorm)] <- 0
  proj_matrix_rownorm <- as.matrix(proj_matrix_rownorm)
  
  # Ensure numeric matrices and preserve row/col names
  rnames <- rownames(proj_matrix_raw)
  cnames <- colnames(proj_matrix_raw)
  
  if (!is.matrix(proj_matrix_raw) || !is.numeric(proj_matrix_raw)) {
    proj_matrix_raw <- as.matrix(sapply(proj_matrix_raw, as.numeric))
  }
  rownames(proj_matrix_raw) <- rnames
  colnames(proj_matrix_raw) <- cnames
  
  if (!is.matrix(proj_matrix_log) || !is.numeric(proj_matrix_log)) {
    proj_matrix_log <- as.matrix(sapply(proj_matrix_log, as.numeric))
  }
  rownames(proj_matrix_log) <- rnames
  colnames(proj_matrix_log) <- cnames
  
  if (!is.matrix(proj_matrix_rownorm) || !is.numeric(proj_matrix_rownorm)) {
    proj_matrix_rownorm <- as.matrix(sapply(proj_matrix_rownorm, as.numeric))
  }
  rownames(proj_matrix_rownorm) <- rnames
  colnames(proj_matrix_rownorm) <- cnames
  
  return(list(
    proj_matrix_raw = proj_matrix_raw,
    mat_log_ordered = proj_matrix_log,
    mat_rownorm_ordered = proj_matrix_rownorm,
    inRH_lookup = inRH_lookup_filtered,
    metadata = metadata_final
  ))
}


################################################################################################################################################################################################################################

# Function to create ipsi-contra matrices from any projection matrix
create_ipsi_contra_from_matrix <- function(proj_matrix, inRH_lookup, matrix_type = "projection") {
  
  cat("\n=== Creating Ipsi/Contra matrices from", matrix_type, "data ===\n")
  
  # Ensure proj_matrix is a data frame for easier manipulation
  if (is.matrix(proj_matrix)) {
    proj_matrix <- as.data.frame(proj_matrix, stringsAsFactors = FALSE)
  }
  
  # Create combined dataframe - ensure we get a data frame
  combined_df <- data.frame(proj_matrix, inRH = inRH_lookup$inRH, stringsAsFactors = FALSE)
  
  cat("Initial data:\n")
  cat("- Total cells:", nrow(combined_df), "\n")
  cat("- LH soma cells:", sum(combined_df$inRH == 0), "\n")
  cat("- RH soma cells:", sum(combined_df$inRH == 1), "\n")
  cat("- Total projection columns:", ncol(combined_df) - 1, "\n")
  
  # Split by soma location
  L_soma <- combined_df[combined_df$inRH == 0, , drop = FALSE]
  R_soma <- combined_df[combined_df$inRH == 1, , drop = FALSE]
  
  # Relabel for L_soma and R_soma
  # For L_soma (LH soma: ipsi = LH, contra = RH)
  colnames(L_soma) <- gsub("_LH(\\.[0-9]+)?$", "-ipsi\\1", colnames(L_soma), perl=TRUE)
  colnames(L_soma) <- gsub("_RH(\\.[0-9]+)?$", "-contra\\1", colnames(L_soma), perl=TRUE)
  
  # For R_soma (RH soma: ipsi = RH, contra = LH)
  colnames(R_soma) <- gsub("_RH(\\.[0-9]+)?$", "-ipsi\\1", colnames(R_soma), perl=TRUE)
  colnames(R_soma) <- gsub("_LH(\\.[0-9]+)?$", "-contra\\1", colnames(R_soma), perl=TRUE)
  
  # Remove inRH column
  L_soma$inRH <- NULL
  R_soma$inRH <- NULL
  
  # Align columns, only keep shared ROIs
  common_cols <- intersect(colnames(L_soma), colnames(R_soma))
  L_soma <- L_soma[, common_cols, drop = FALSE]
  R_soma <- R_soma[, common_cols, drop = FALSE]
  
  cat("After matching:\n")
  cat("- Common projection columns:", length(common_cols), "\n")
  
  # Combine all cells
  ipsi_contra <- rbind(L_soma, R_soma)
  
  cat("Final ipsi-contra matrix (", matrix_type, "):\n")
  cat("- Total cells:", nrow(ipsi_contra), "\n")
  cat("- Total regions:", ncol(ipsi_contra), "\n")
  
  return(ipsi_contra)
}

# Convenience wrapper functions for backward compatibility and clarity
create_ipsi_contra_from_raw <- function(proj_matrix_raw, inRH_lookup) {
  return(create_ipsi_contra_from_matrix(proj_matrix_raw, inRH_lookup, "raw counts"))
}

create_ipsi_contra_from_log <- function(proj_matrix_log, inRH_lookup) {
  return(create_ipsi_contra_from_matrix(proj_matrix_log, inRH_lookup, "log-transformed"))
}

create_ipsi_contra_from_rownorm <- function(proj_matrix_rownorm, inRH_lookup) {
  return(create_ipsi_contra_from_matrix(proj_matrix_rownorm, inRH_lookup, "row-normalized"))
}
