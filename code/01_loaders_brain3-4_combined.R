# Loader function for brain #3 and brain #4 MAPseq data where each dataset gets pre-processed and normalized, and then combined into one
# ROI order is indexed by slide it was collected from to enhance compatibility across datasets and lost samples

################################################################################################################################################################################################################################
# Function to pre-process and load the data in a format directly usable by downstream analyses - 780345 brain #3
load_data_brain3 <-function() {
  
  # Set working directory
  setwd('/scratch/BARseq_780345/')
  
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
  proj_index <- readr::read_tsv("sampleinfo.tsv")
  head(proj_index)
  sampleinfo <- read_excel("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/M295_20250721.sampleinfo.xlsx", sheet = "Sample information", skip = 1, range = "A2:G110")
  head(sampleinfo)
  colnames(sampleinfo) <- c("usertube", "ourtube", "samplename", "siteinfo", "QC_qPCR", "rtprimer", "brain")
  proj_index$MapSeqV1_tube <- sampleinfo$rtprimer[match(proj_index$usertube, sampleinfo$usertube)]
  head(proj_index)
  
  # Need to create a vector to index samples based on slide they came from to circumvent erroneous indexing due lost samples
  design_df <- proj_index %>%
    filter(brain == 780345)
  nrow(design_df)  # Should be 108, as this includes 7 samples lost during repeated processing
  slide_counts <- c(2, 2, 6, 6, 14, 14, 14, 14, 16, 5, 4, 4, 4, 3) 
  stopifnot(sum(slide_counts) == nrow(design_df))
  slide_vector <- rep(seq_along(slide_counts), times = slide_counts)
  length(slide_vector)
  design_df$slide <- slide_vector
  design_df$hemisphere <- design_df$samplename
  design_df$ROI <- design_df$region
  proj_index <- proj_index %>%
    left_join(
      design_df %>% select(rtprimer, ROI, hemisphere, slide),
      by = "rtprimer"
    )
  
  # Rename BC columns to region_samplename_unique
  bc_cols <- grep("^BC\\d+$", colnames(projection_matrix), value = TRUE)
  rtprimer_numbers <- as.integer(sub("BC", "", bc_cols))
  # Index metadata from proj_index using rtprimer
  hemi  <- proj_index$hemisphere[match(rtprimer_numbers, proj_index$MapSeqV1_tube)]
  roi   <- proj_index$ROI[match(rtprimer_numbers, proj_index$MapSeqV1_tube)]
  slide <- proj_index$slide[match(rtprimer_numbers, proj_index$MapSeqV1_tube)]
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
  
  # Ensure numeric matrices and preserve row/col names
  rnames <- rownames(proj_matrix_raw)
  cnames <- colnames(proj_matrix_raw)
  
  if (!is.matrix(proj_matrix_raw) || !is.numeric(proj_matrix_raw)) {
    proj_matrix_raw <- as.matrix(sapply(proj_matrix_raw, as.numeric))
  }
  rownames(proj_matrix_raw) <- rnames
  colnames(proj_matrix_raw) <- cnames
  
  # Check for olf bulb mislabeling and fix if needed
  RH_soma_mask <- inRH_lookup_filtered$inRH == 1
  LH_soma_mask <- inRH_lookup_filtered$inRH == 0
  
  # Find all olfactory bulb columns (including numbered variants)
  olf_bulb_cols <- grep("^olf bulb", colnames(proj_matrix_raw), value = TRUE)
  olf_bulb_RH_cols <- grep("^olf bulb.*_RH", olf_bulb_cols, value = TRUE)
  olf_bulb_LH_cols <- grep("^olf bulb.*_LH", olf_bulb_cols, value = TRUE)
  
  if(length(olf_bulb_RH_cols) > 0 && length(olf_bulb_LH_cols) > 0) {
    cat("\n=== Olfactory Bulb Analysis ===\n")
    cat("Found olfactory bulb columns:\n")
    cat("RH columns:", paste(olf_bulb_RH_cols, collapse = ", "), "\n")
    cat("LH columns:", paste(olf_bulb_LH_cols, collapse = ", "), "\n")
    
    # Get cells projecting to each olfactory bulb region
    olf_projecting_cells <- list()
    
    cat("\n--- Cells projecting to each olfactory bulb region ---\n")
    for(col in c(olf_bulb_RH_cols, olf_bulb_LH_cols)) {
      cells <- which(proj_matrix_raw[, col] > 0)
      olf_projecting_cells[[col]] <- cells
      cat(col, ": ", length(cells), " cells\n", sep="")
    }
    
    # Calculate overlaps between all pairs of olfactory bulb regions
    cat("\n--- Overlap Analysis Between Olfactory Bulb Projecting Cells ---\n")
    olf_cols <- c(olf_bulb_RH_cols, olf_bulb_LH_cols)
    
    for(i in 1:(length(olf_cols)-1)) {
      for(j in (i+1):length(olf_cols)) {
        col1 <- olf_cols[i]
        col2 <- olf_cols[j]
        
        cells1 <- olf_projecting_cells[[col1]]
        cells2 <- olf_projecting_cells[[col2]]
        
        overlap <- intersect(cells1, cells2)
        unique_to_1 <- setdiff(cells1, cells2)
        unique_to_2 <- setdiff(cells2, cells1)
        
        cat("\n", col1, " vs ", col2, ":\n", sep="")
        cat("  ", col1, " only: ", length(unique_to_1), " cells\n", sep="")
        cat("  ", col2, " only: ", length(unique_to_2), " cells\n", sep="")
        cat("  Shared cells: ", length(overlap), " cells\n", sep="")
        
        if(length(cells1) > 0) {
          pct_shared_of_1 <- length(overlap) / length(cells1) * 100
          cat("  Shared as % of ", col1, " projectors: ", round(pct_shared_of_1, 1), "%\n", sep="")
        }
        
        if(length(cells2) > 0) {
          pct_shared_of_2 <- length(overlap) / length(cells2) * 100
          cat("  Shared as % of ", col2, " projectors: ", round(pct_shared_of_2, 1), "%\n", sep="")
        }
        
        # Calculate Jaccard index (overlap / union)
        union_size <- length(union(cells1, cells2))
        if(union_size > 0) {
          jaccard <- length(overlap) / union_size * 100
          cat("  Jaccard similarity: ", round(jaccard, 1), "%\n", sep="")
        }
      }
    }
    
    # Special focus on potential swapping pairs (same base region)
    cat("\n--- Sample Swap Detection (Corresponding Regions) ---\n")
    rh_bases <- gsub("_RH.*", "", olf_bulb_RH_cols)
    lh_bases <- gsub("_LH.*", "", olf_bulb_LH_cols)
    
    for(i in seq_along(olf_bulb_RH_cols)) {
      rh_col <- olf_bulb_RH_cols[i]
      rh_base <- rh_bases[i]
      
      matching_lh <- olf_bulb_LH_cols[lh_bases == rh_base]
      
      if(length(matching_lh) > 0) {
        for(lh_col in matching_lh) {
          rh_cells <- olf_projecting_cells[[rh_col]]
          lh_cells <- olf_projecting_cells[[lh_col]]
          shared_cells <- intersect(rh_cells, lh_cells)
          
          cat("\n", rh_col, " vs ", lh_col, " (corresponding regions):\n", sep="")
          cat("  Total unique cells: ", length(union(rh_cells, lh_cells)), "\n", sep="")
          cat("  ", rh_col, " projectors: ", length(rh_cells), "\n", sep="")
          cat("  ", lh_col, " projectors: ", length(lh_cells), "\n", sep="")
          cat("  Shared projectors: ", length(shared_cells), "\n", sep="")
          
          if(length(rh_cells) > 0 && length(lh_cells) > 0) {
            pct_shared_rh <- length(shared_cells) / length(rh_cells) * 100
            pct_shared_lh <- length(shared_cells) / length(lh_cells) * 100
            cat("  % of ", rh_col, " cells that also project to ", lh_col, ": ", round(pct_shared_rh, 1), "%\n", sep="")
            cat("  % of ", lh_col, " cells that also project to ", rh_col, ": ", round(pct_shared_lh, 1), "%\n", sep="")
            
            # High overlap suggests potential sample swap
            if(pct_shared_rh > 50 || pct_shared_lh > 50) {
              cat("  *** HIGH OVERLAP - POTENTIAL SAMPLE SWAP? ***\n")
            }
          }
        }
      }
    }
    
    # Updated mislabeling check for indexed columns
    # Extract unique indices from olfactory bulb columns
    olf_indices <- unique(gsub("olf bulb_(LH|RH)\\.(.*)", "\\2", olf_bulb_cols))
    
    for(index in olf_indices) {
      rh_col <- paste0("olf bulb_RH.", index)
      lh_col <- paste0("olf bulb_LH.", index)
      
      if(rh_col %in% colnames(proj_matrix_raw) && lh_col %in% colnames(proj_matrix_raw)) {
        mean_RH_soma_RH <- mean(proj_matrix_raw[RH_soma_mask, rh_col], na.rm = TRUE)
        mean_RH_soma_LH <- mean(proj_matrix_raw[RH_soma_mask, lh_col], na.rm = TRUE)
        mean_LH_soma_RH <- mean(proj_matrix_raw[LH_soma_mask, rh_col], na.rm = TRUE)
        mean_LH_soma_LH <- mean(proj_matrix_raw[LH_soma_mask, lh_col], na.rm = TRUE)
        
        cat("\n=== Mean Projection Analysis for index", index, "===\n")
        cat("RH soma: mean", rh_col, "=", round(mean_RH_soma_RH, 3), "; mean", lh_col, "=", round(mean_RH_soma_LH, 3), "\n")
        cat("LH soma: mean", rh_col, "=", round(mean_LH_soma_RH, 3), "; mean", lh_col, "=", round(mean_LH_soma_LH, 3), "\n")
        
        # If means indicate a switch, swap the columns in all matrices
        if (mean_RH_soma_RH < mean_RH_soma_LH && mean_LH_soma_RH > mean_LH_soma_LH) {
          # Swap in raw matrix
          tmp <- proj_matrix_raw[, rh_col]
          proj_matrix_raw[, rh_col] <- proj_matrix_raw[, lh_col]
          proj_matrix_raw[, lh_col] <- tmp
          
          cat("Swapped", rh_col, "and", lh_col, "columns due to detected switch.\n")
        }
      }
    }
  }
  
  return(list(
    proj_matrix_raw = proj_matrix_raw,
    inRH_lookup = inRH_lookup_filtered,
    metadata = metadata_final
  ))
}

################################################################################################################################################################################################################################
# Function to pre-process and load the data in a format directly usable by downstream analyses - 780346 brain #4
load_data_brain4 <-function() {
  
  # Set working directory
  setwd('/scratch/BARseq_780346/')
  
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
  proj_index <- read_tsv("/data/780346_2025-06-11_00-00-00/MAPseq/M305_20251030_USEthis/M305sampleinfo.tsv")
  head(proj_index)
  sampleinfo <- read_excel("./M305_sample_information.xlsx", sheet = "Sample information", skip = 1,range = "A2:J122")
  head(sampleinfo)
  colnames(sampleinfo) <- c("usertube", "ourtube", "samplename", "siteinfo", "QC_qPCR", "rtprimer", "brain", "hemisphere", "ROI", "notes")
  
  # Need to create a vector to index samples based on slide they came from to circumvent erroneous indexing due lost samples
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
  
  # Ensure numeric matrices and preserve row/col names
  rnames <- rownames(proj_matrix_raw)
  cnames <- colnames(proj_matrix_raw)
  
  if (!is.matrix(proj_matrix_raw) || !is.numeric(proj_matrix_raw)) {
    proj_matrix_raw <- as.matrix(sapply(proj_matrix_raw, as.numeric))
  }
  rownames(proj_matrix_raw) <- rnames
  colnames(proj_matrix_raw) <- cnames
  
  return(list(
    proj_matrix_raw = proj_matrix_raw,
    inRH_lookup = inRH_lookup_filtered,
    metadata = metadata_final
  ))
}

################################################################################################################################################################################################################################
# Function to sum by base ROI + hemisphere and output a mapping list
sum_by_base_roi <- function(proj_matrix) {
  if (is.matrix(proj_matrix)) proj_matrix <- as.data.frame(proj_matrix)
  
  # Extract base ROI + hemisphere (e.g., "olf bulb_RH" from "olf bulb_RH.1")
  base_cols <- stringr::str_replace(colnames(proj_matrix), "\\..*$", "")  # Remove slide (.1, etc.)
  
  # Get unique base ROIs
  unique_bases <- unique(base_cols)
  
  # Create summed matrix
  summed <- matrix(0, nrow = nrow(proj_matrix), ncol = length(unique_bases))
  colnames(summed) <- unique_bases
  
  # Preserve rownames explicitly
  original_rownames <- rownames(proj_matrix)
  if (is.null(original_rownames)) {
    original_rownames <- 1:nrow(proj_matrix)  # Fallback
  }
  rownames(summed) <- original_rownames
  
  for (i in seq_along(unique_bases)) {
    base <- unique_bases[i]
    cols <- which(base_cols == base)
    if (length(cols) > 0) {
      summed[, i] <- rowSums(proj_matrix[, cols, drop = FALSE])
    }
  }
  
  # Ensure the returned data.frame has rownames
  result <- as.data.frame(summed)
  rownames(result) <- original_rownames
  
  # Create mapping list: for each base, list the original columns summed
  mapping <- list()
  for (i in seq_along(unique_bases)) {
    base <- unique_bases[i]
    cols <- which(base_cols == base)
    if (length(cols) > 0) {
      mapping[[base]] <- colnames(proj_matrix)[cols]  # Store original column names
    }
  }
  
  # Return a list with the summed matrix and the mapping
  return(list(summed_matrix = result, mapping = mapping))
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

# Convenience wrapper function for backward compatibility and clarity
create_ipsi_contra_from_raw <- function(proj_matrix_raw, inRH_lookup) {
  return(create_ipsi_contra_from_matrix(proj_matrix_raw, inRH_lookup, "raw counts"))
}

################################################################################################################################################################################################################################
# Function to normalize projection matrix: row norm, log transform, global max norm
normalize_projection_matrix <- function(proj_matrix, matrix_name = "unknown") {
  
  cat("\n=== Normalizing", matrix_name, "===\n")
  cat("Original matrix dimensions:", dim(proj_matrix), "\n")
  cat("Original range:", range(proj_matrix, na.rm = TRUE), "\n")
  
  # Step 1: Row normalization (normalize by total abundance per cell)
  row_sums <- rowSums(proj_matrix, na.rm = TRUE)
  # Avoid division by zero for rows with all zeros
  row_sums[row_sums == 0] <- 1
  proj_matrix <- proj_matrix / row_sums
  cat("After row normalization - range:", range(proj_matrix, na.rm = TRUE), "\n")
  
  # Step 2: Log transform (log1p to handle zeros and small values)
  proj_matrix <- log1p(proj_matrix)  # log(1 + x)
  cat("After log transform - range:", range(proj_matrix, na.rm = TRUE), "\n")
  
  # Step 3: Normalize by global max
  global_max <- max(proj_matrix, na.rm = TRUE)
  if (global_max > 0) {
    proj_matrix <- proj_matrix / global_max
  }
  cat("After global max normalization - range:", range(proj_matrix, na.rm = TRUE), "\n")
  
  cat("Final matrix dimensions:", dim(proj_matrix), "\n")
  return(proj_matrix)
}