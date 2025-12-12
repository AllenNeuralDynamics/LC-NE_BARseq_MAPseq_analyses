# Function to pre-process and load the data in a format directly usable by downstream analyses
load_data <-function() {

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
  proj_index <- readr::read_tsv("sampleinfo.tsv")
  sampleinfo <- read_excel("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/M295_20250721.sampleinfo.xlsx", sheet = "Sample information")
  colnames(sampleinfo) <- c("usertube", "ourtube", "samplename", "siteinfo", "QC_qPCR", "rtprimer", "brain")
  proj_index$MapSeqV1_tube <- sampleinfo$rtprimer[match(proj_index$usertube, sampleinfo$usertube)]
  
  # Rename BC columns to region_samplename_unique
  bc_cols <- grep("^BC\\d+$", colnames(projection_matrix), value = TRUE)
  rtprimer_numbers <- as.integer(sub("BC", "", bc_cols))
  region_names <- proj_index$region[match(rtprimer_numbers, proj_index$MapSeqV1_tube)]
  samplenames  <- proj_index$samplename[match(rtprimer_numbers, proj_index$MapSeqV1_tube)]
  region_samplename <- paste(region_names, samplenames, sep = "_")
  region_samplename_unique <- make.unique(region_samplename)
  colnames(projection_matrix)[match(bc_cols, colnames(projection_matrix))] <- region_samplename_unique
  
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
    
    # Original mislabeling check using first found columns
    if("olf bulb_RH" %in% colnames(proj_matrix_raw) && "olf bulb_LH" %in% colnames(proj_matrix_raw)) {
      mean_RH_soma_RH <- mean(proj_matrix_raw[RH_soma_mask, "olf bulb_RH"], na.rm = TRUE)
      mean_RH_soma_LH <- mean(proj_matrix_raw[RH_soma_mask, "olf bulb_LH"], na.rm = TRUE)
      mean_LH_soma_RH <- mean(proj_matrix_raw[LH_soma_mask, "olf bulb_RH"], na.rm = TRUE)
      mean_LH_soma_LH <- mean(proj_matrix_raw[LH_soma_mask, "olf bulb_LH"], na.rm = TRUE)
      
      cat("\n=== Mean Projection Analysis ===\n")
      cat("RH soma: mean olf bulb_RH =", round(mean_RH_soma_RH, 3), "; mean olf bulb_LH =", round(mean_RH_soma_LH, 3), "\n")
      cat("LH soma: mean olf bulb_RH =", round(mean_LH_soma_RH, 3), "; mean olf bulb_LH =", round(mean_LH_soma_LH, 3), "\n")
      
      # If means indicate a switch, swap the columns in all matrices
      if (mean_RH_soma_RH < mean_RH_soma_LH && mean_LH_soma_RH > mean_LH_soma_LH) {
        # Swap in raw matrix
        tmp <- proj_matrix_raw[, "olf bulb_RH"]
        proj_matrix_raw[, "olf bulb_RH"] <- proj_matrix_raw[, "olf bulb_LH"]
        proj_matrix_raw[, "olf bulb_LH"] <- tmp
        
        # Swap in log matrix
        tmp <- proj_matrix_log[, "olf bulb_RH"]
        proj_matrix_log[, "olf bulb_RH"] <- proj_matrix_log[, "olf bulb_LH"]
        proj_matrix_log[, "olf bulb_LH"] <- tmp
        
        # Swap in row-normalized matrix
        tmp <- proj_matrix_rownorm[, "olf bulb_RH"]
        proj_matrix_rownorm[, "olf bulb_RH"] <- proj_matrix_rownorm[, "olf bulb_LH"]
        proj_matrix_rownorm[, "olf bulb_LH"] <- tmp
        
        cat("Swapped olf bulb_RH and olf bulb_LH columns due to detected switch.\n")
      }
    }  # Close the if("olf bulb_RH" %in% colnames...) block
  }    # Close the if(length(olf_bulb_RH_cols) > 0...) block
  
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
