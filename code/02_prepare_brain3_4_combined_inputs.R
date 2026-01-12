# 02_prepare_brain3_4_combined_inputs.R
# Purpose: produce combined brain3+brain4 ipsi/contra normalized matrix + aligned metadata
# Depends on: 01_loaders_brain3-4_combined.R (defines load_data_brain3/load_data_brain4 + helpers)

prepare_brain3_4_inputs <- function(
    loaders_path = "~/capsule/code/01_loaders_brain3-4_combined.R",
    out_dir = "/scratch/BARseq_780345-780346_combined/",
    do_cerebellum_merge = TRUE,
    cerebellum_groups = c("RH.12", "LH.12"),
    fix_amyg_gpe = TRUE,
    return_intermediates = FALSE,
    verbose = TRUE,
    restore_wd = TRUE          # <-- ADD THIS
) {
  
  if (!file.exists(loaders_path)) stop("Missing loaders_path: ", loaders_path)
  source(loaders_path)
  
  old_wd <- getwd()
  if (restore_wd) {            # <-- CHANGE THIS
    on.exit(setwd(old_wd), add = TRUE)
  }
  
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  setwd(out_dir)
  
  as_numeric_strict <- function(x, label = "value") {
    na0 <- sum(is.na(x))
    y <- suppressWarnings(as.numeric(as.character(x)))
    na1 <- sum(is.na(y))
    if (na1 > na0) {
      stop("Numeric coercion introduced ", (na1 - na0), " new NA(s) in ", label,
           ". Example original values: ",
           paste(head(unique(as.character(x)[is.na(y) & !is.na(x)]), 10), collapse = ", "))
    }
    y
  }
  
  # ---------- Load both brains ----------
  brain3 <- load_data_brain3()
  brain4 <- load_data_brain4()
  setwd(out_dir)  # restore, because loaders change wd internally
  
  # Basic structure checks
  req_names <- c("proj_matrix_raw", "inRH_lookup", "metadata")
  if (!all(req_names %in% names(brain3))) stop("brain3 loader missing: ", paste(setdiff(req_names, names(brain3)), collapse = ", "))
  if (!all(req_names %in% names(brain4))) stop("brain4 loader missing: ", paste(setdiff(req_names, names(brain4)), collapse = ", "))
  
  # Extract components (preserve rownames exactly)
  brain3_raw <- as.data.frame(brain3$proj_matrix_raw, stringsAsFactors = FALSE)
  rownames(brain3_raw) <- rownames(brain3$proj_matrix_raw)
  brain3_lookup <- brain3$inRH_lookup
  brain3_meta <- brain3$metadata
  
  brain4_raw <- as.data.frame(brain4$proj_matrix_raw, stringsAsFactors = FALSE)
  rownames(brain4_raw) <- rownames(brain4$proj_matrix_raw)
  brain4_lookup <- brain4$inRH_lookup
  brain4_meta <- brain4$metadata
  
  req_meta_cols <- c("slice", "barcode", "CCF_DV", "CCF_ML", "CCF_AP", "louvain_cluster", "CellID", "inRH")
  missing3 <- setdiff(req_meta_cols, colnames(brain3_meta))
  missing4 <- setdiff(req_meta_cols, colnames(brain4_meta))
  if (length(missing3) > 0) stop("brain3_meta missing columns: ", paste(missing3, collapse = ", "))
  if (length(missing4) > 0) stop("brain4_meta missing columns: ", paste(missing4, collapse = ", "))
  
  # Fail-fast: rownames must exist and be unique
  if (is.null(rownames(brain3_raw)) || any(rownames(brain3_raw) == "")) stop("brain3_raw has missing/empty rownames (row_id).")
  if (is.null(rownames(brain4_raw)) || any(rownames(brain4_raw) == "")) stop("brain4_raw has missing/empty rownames (row_id).")
  if (anyDuplicated(rownames(brain3_raw))) stop("brain3_raw rownames not unique (row_id).")
  if (anyDuplicated(rownames(brain4_raw))) stop("brain4_raw rownames not unique (row_id).")
  
  # Fail-fast: lookup must be aligned-able
  if (!all(c("row_id", "inRH") %in% colnames(brain3_lookup))) stop("brain3_inRH_lookup must contain row_id,inRH.")
  if (!all(c("row_id", "inRH") %in% colnames(brain4_lookup))) stop("brain4_inRH_lookup must contain row_id,inRH.")
  if (anyDuplicated(brain3_lookup$row_id)) stop("brain3_inRH_lookup$row_id not unique.")
  if (anyDuplicated(brain4_lookup$row_id)) stop("brain4_inRH_lookup$row_id not unique.")
  
  # Fail-fast: metadata rownames should be row_id (recommended)
  if (is.null(rownames(brain3_meta)) || any(rownames(brain3_meta) == "")) {
    stop("brain3_metadata must have rownames=row_id for safe alignment downstream.")
  }
  if (is.null(rownames(brain4_meta)) || any(rownames(brain4_meta) == "")) {
    stop("brain4_metadata must have rownames=row_id for safe alignment downstream.")
  }
  if (anyDuplicated(rownames(brain3_meta))) stop("brain3_metadata rownames not unique (row_id).")
  if (anyDuplicated(rownames(brain4_meta))) stop("brain4_metadata rownames not unique (row_id).")
  
  # ---------- ROI compatibility tweaks on brain4 raw ----------
  if (fix_amyg_gpe) {
    colnames(brain4_raw) <- stringr::str_replace_all(colnames(brain4_raw), "amyg-Gpe", "amyg-GPe")
  }
  
  if (do_cerebellum_merge) {
    for (group in cerebellum_groups) {
      i_col   <- paste0("cerebellum I_", group)
      ii_col  <- paste0("cerebellum II_", group)
      new_col <- paste0("cerebellum_", group)
      
      if (i_col %in% colnames(brain4_raw) && ii_col %in% colnames(brain4_raw)) {
        # Sum
        brain4_raw[[new_col]] <-
          as_numeric_strict(brain4_raw[[i_col]],  i_col) +
          as_numeric_strict(brain4_raw[[ii_col]], ii_col)
        # Drop originals
        brain4_raw <- brain4_raw[, !(colnames(brain4_raw) %in% c(i_col, ii_col)), drop = FALSE]
        if (verbose) cat("Summed", i_col, "+", ii_col, "into", new_col, "\n")
      } else {
        if (verbose) cat("Warning: Columns for", group, "not found (", i_col, ", ", ii_col, ")\n", sep = "")
      }
    }
  }
  
  # ---------- Sum by base ROI (removes slide suffix) ----------
  brain3_sum_res <- sum_by_base_roi(brain3_raw)
  brain4_sum_res <- sum_by_base_roi(brain4_raw)
  brain3_summed  <- brain3_sum_res$summed_matrix
  brain4_summed  <- brain4_sum_res$summed_matrix
  
  # Preserve rownames (sum_by_base_roi already does, but enforce)
  rownames(brain3_summed) <- rownames(brain3_raw)
  rownames(brain4_summed) <- rownames(brain4_raw)
  
  # ---------- Ipsi/contra conversion (aligns inRH by row_id inside your edited function) ----------
  ipsi3 <- create_ipsi_contra_from_matrix(brain3_summed, brain3_lookup, "brain3 summed")
  ipsi4 <- create_ipsi_contra_from_matrix(brain4_summed, brain4_lookup, "brain4 summed")
  
  # ---------- Shared regions ----------
  brain3_cols <- colnames(ipsi3)
  brain4_cols <- colnames(ipsi4)
  
  shared_regions <- intersect(brain3_cols, brain4_cols)
  unique_brain3 <- setdiff(brain3_cols, brain4_cols)
  unique_brain4 <- setdiff(brain4_cols, brain3_cols)
  
  if (length(shared_regions) == 0) {
    stop("No shared regions between brain3 and brain4 after ipsi/contra. Check ROI naming / hemisphere tokens.")
  }
  
  # Subset to shared
  brain3_shared <- ipsi3[, shared_regions, drop = FALSE]
  brain4_shared <- ipsi4[, shared_regions, drop = FALSE]
  
  # ---------- Normalize (your edited normalizer preserves dimnames) ----------
  brain3_norm <- normalize_projection_matrix(brain3_shared, "brain3 shared ipsi-contra")
  brain4_norm <- normalize_projection_matrix(brain4_shared, "brain4 shared ipsi-contra")
  
  if (!all(is.finite(brain3_norm))) stop("brain3_norm contains non-finite values (Inf/NaN).")
  if (!all(is.finite(brain4_norm))) stop("brain4_norm contains non-finite values (Inf/NaN).")
  
  # Combine
  combined_norm <- rbind(brain3_norm, brain4_norm)
  
  # Ensure rownames exist and unique after combine
  if (is.null(rownames(combined_norm)) || any(rownames(combined_norm) == "")) stop("combined_norm missing/empty rownames.")
  if (anyDuplicated(rownames(combined_norm))) stop("combined_norm rownames not unique after rbind (row_id collision across brains).")
  
  # ---------- Metadata alignment (CRITICAL) ----------
  combined_metadata <- rbind(brain3_meta, brain4_meta)
  if (anyDuplicated(rownames(combined_metadata))) {
    dup <- rownames(combined_metadata)[duplicated(rownames(combined_metadata))]
    stop("combined_metadata has duplicated row_id(s) (match() ambiguous). Example: ",
         paste(head(unique(dup), 20), collapse = ", "))
  }
  # Align metadata to combined_norm row order by row_id (rownames)
  idx_meta <- match(rownames(combined_norm), rownames(combined_metadata))
  if (anyNA(idx_meta)) {
    miss <- rownames(combined_norm)[is.na(idx_meta)]
    stop("combined_metadata missing ", sum(is.na(idx_meta)), " row_id(s) present in combined_norm. Example: ",
         paste(head(miss, 20), collapse = ", "))
  }
  combined_metadata_aligned <- combined_metadata[idx_meta, , drop = FALSE]
  stopifnot(identical(rownames(combined_metadata_aligned), rownames(combined_norm)))
  
  # Optional: align inRH lookup too (sometimes useful downstream)
  align_lookup <- function(lookup, ids, label = "lookup") {
    idx <- match(ids, lookup$row_id)
    if (anyNA(idx)) {
      miss <- ids[is.na(idx)]
      stop(label, " missing ", sum(is.na(idx)), " row_id(s). Example: ",
           paste(head(miss, 20), collapse = ", "))
    }
    out <- lookup[idx, , drop = FALSE]
    if (!identical(out$row_id, ids)) stop(label, " alignment produced unexpected order (bug).")
    out
  }
  
  combined_lookup_aligned <- rbind(
    align_lookup(brain3_lookup, rownames(brain3_norm), "brain3_lookup vs brain3_norm"),
    align_lookup(brain4_lookup, rownames(brain4_norm), "brain4_lookup vs brain4_norm")
  )
  stopifnot(identical(combined_lookup_aligned$row_id, rownames(combined_norm)))
  rownames(combined_lookup_aligned) <- combined_lookup_aligned$row_id
  
  if (verbose) {
    cat("\n=== Combined outputs ===\n")
    cat("combined_norm dim:", dim(combined_norm), "\n")
    cat("combined_metadata dim:", dim(combined_metadata_aligned), "\n")
    cat("shared regions:", length(shared_regions), "\n")
    cat("unique brain3 regions:", length(unique_brain3), "\n")
    cat("unique brain4 regions:", length(unique_brain4), "\n")
  }
  
  out <- list(
    combined_norm = combined_norm,
    combined_metadata = combined_metadata_aligned,
    combined_inRH_lookup = combined_lookup_aligned,
    shared_regions = shared_regions,
    unique_brain3 = unique_brain3,
    unique_brain4 = unique_brain4,
    brain3_mapping = brain3_sum_res$mapping,
    brain4_mapping = brain4_sum_res$mapping
  )
  
  if (return_intermediates) {
    out$brain3 <- brain3
    out$brain4 <- brain4
    out$brain3_metadata <- brain3_meta
    out$brain4_metadata <- brain4_meta
    out$brain3_inRH_lookup <- brain3_lookup
    out$brain4_inRH_lookup <- brain4_lookup
    out$brain3_raw <- brain3_raw
    out$brain4_raw <- brain4_raw
    out$brain3_summed <- brain3_summed
    out$brain4_summed <- brain4_summed
    out$ipsi_contra_brain3 <- ipsi3
    out$ipsi_contra_brain4 <- ipsi4
    out$brain3_norm <- brain3_norm
    out$brain4_norm <- brain4_norm
  }
  
  return(out)
}
