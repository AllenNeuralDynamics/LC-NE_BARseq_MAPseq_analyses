################# MAPseq-BARseq Barcode Matching Pipeline ##################
# Performs converts BARseq barcodes to MAPseq format, truncates MAPseq 32nt barcodes to 15nt to match BARseq barcode length
# Identifies uniquely barcoded BARseq cells which have a match in the MAPseq dataset with 0,1,2,or 3 Hamming distance mismatches allowed

#set working directory 
setwd('/results/BARseq_780345/')
############################################################################################################################################################################################################
# Sanity checks for QC of the MAPseq data
UMI_filt <- read_tsv("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/780345.nbcm.tsv")
UMI_raw <- read_tsv("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/780345.rbcm.tsv")
UMI_spikein <- read_tsv("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/780345.sbcm.tsv")

# Process all three datasets
UMI_filt <- UMI_filt %>%
  mutate(umi_count = rowSums(select(., starts_with("BC"))))

UMI_raw <- UMI_raw %>%
  mutate(umi_count = rowSums(select(., starts_with("BC"))))

UMI_spikein <- UMI_spikein %>%
  mutate(umi_count = rowSums(select(., starts_with("BC"))))

# Use consistent binwidth approach for all three plots
p1 <- ggplot(UMI_filt, aes(x = umi_count)) + 
  geom_histogram(bins = 100, fill = "steelblue", color = "white", linewidth = 0.1) +
  scale_x_log10() +
  labs(title = "Filtered UMI Counts (log10 scale)",
       x = "UMI Count",
       y = "Frequency") +
  theme_minimal()

p2 <- ggplot(UMI_raw, aes(x = umi_count)) +
  geom_histogram(bins = 100, fill = "coral", color = "white", linewidth = 0.1) + 
  scale_x_log10() +
  labs(title = "Raw UMI Counts (log10 scale)",
       x = "UMI Count",
       y = "Frequency") +
  theme_minimal()

p3 <- ggplot(UMI_spikein %>% filter(umi_count > 1), aes(x = umi_count)) +
  geom_histogram(bins = 15, fill = "darkgreen", color = "white", linewidth = 0.5) +
  labs(title = paste("Spike-in UMI Counts >1 - n =", nrow(UMI_spikein %>% filter(umi_count > 1))),
       x = "UMI Count", 
       y = "Frequency") +
  theme_minimal()

plot <- p1 + p2 + p3
print(plot)
ggsave("MapSeqV1_UMIcounts_QCchecks.pdf", plot = plot, device = "pdf", width = 18, height = 6)  

rm(UMI_raw, UMI_filt,UMI_spikein)

############################################################################################################################################################################################################
# BarSeq input which includes all putative LC-NE cells
LCNEcluster <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
dim(LCNEcluster)
colData(LCNEcluster)
table(colData(LCNEcluster)$barcode)
table(colData(LCNEcluster)$louvain_cluster)

############################################################################################################################################################################################################
# Process barcode information for all barcoded cells from BARseq experiment and convert them to MAPseq basecalls format
# Read CSV without column names, and name columns manually
barcodes_raw <- read_csv("/data/BARseq_soma_barcodes/barcodes_BC_qc_780345.csv", col_names = FALSE)
# Rename for clarity
colnames(barcodes_raw)[1] <- "CellID"
colnames(barcodes_raw)[2:16] <- paste0("B", 1:15)  # Barcode positions
# Mapping function: 1=G, 2=T, 3=A, 4=C
base_map <- c("G", "T", "A", "C")
# Convert numeric basecalls to nucleotide letters
convert_to_nucleotides <- function(x) base_map[as.integer(x)]
# Apply across barcode columns (B1–B15)
barseq_nuc <- as.data.frame(lapply(barcodes_raw[, 2:16], convert_to_nucleotides))
# Add CellID and VBC
barcodes_raw$VBC <- apply(barseq_nuc, 1, paste0, collapse = "")
# Check a few rows
head(barcodes_raw[, c("CellID", "VBC")])

############################################################################################################################################################################################################
# Match the barcode sequences to cells in LCNE cluster based on uid
# Ensure 'uid' exists in barseq_raw
if (!"CellID" %in% colnames(barcodes_raw)) {
  stop("No 'uid' column found in barseq_raw.")
}
# Extract 'uid' of all barcoded cells from LCNE cluster object
good_uids <- colData(LCNEcluster)$uid[colData(LCNEcluster)$barcode == TRUE]
# Subset barcodes_raw using the good 'uid'
LC_barcoded_cells <- barcodes_raw[barcodes_raw$CellID %in% good_uids, ]
# Check the number of rows in the subsetted object
nrow(barcodes_raw)
nrow(LC_barcoded_cells)
# Save the barcode information for LCNE cells
write.csv(LC_barcoded_cells, file = "LCNE_barcoded_cells.csv", row.names = FALSE)
rm(barcodes_raw, barseq_nuc)

############################################################################################################################################################################################################
barseq <- LC_barcoded_cells

# MapSeq input, truncate to first 15 characters to ensure that vector legnth is equivalent to BarSeq 15 cycles
mapseq <- read_tsv("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/780345.nbcm.tsv")

# Find which columns are logical (shouldn't be)
logical_cols <- sapply(mapseq, is.logical)
# Convert them to numeric
mapseq[logical_cols] <- lapply(mapseq[logical_cols], as.numeric)
# check for naming consistency
if ("vbc_read_col" %in% names(mapseq)) {
  names(mapseq)[names(mapseq) == "vbc_read_col"] <- "vbc_read"
}
mapseq_vbcs <- mapseq$vbc_read

mapseq_vbcs <- substr(mapseq_vbcs, 1, 15)
mapseq_split <- strsplit(mapseq_vbcs, "")

# Hamming distance vectorized
hamming_distance_vec <- function(str1_chars, vec2_split) {
  sapply(vec2_split, function(str2_chars) sum(str1_chars != str2_chars))
}

# Function to process a single batch (no longer parallel internally, but called in parallel)
process_batch <- function(batch_df, mapseq_split, mapseq_vbcs, max_dist) {
  result_list <- list()
  for (i in seq_len(nrow(batch_df))) {
    vbc1_chars <- strsplit(batch_df$VBC[i], "")[[1]]
    dists <- hamming_distance_vec(vbc1_chars, mapseq_split)
    matched_idx <- which(dists <= max_dist)
    if (length(matched_idx) > 0) {
      result_list[[i]] <- data.frame(
        CellID = rep(batch_df$CellID[i], length(matched_idx)),
        vbc_read = mapseq_vbcs[matched_idx],
        dist = dists[matched_idx]
      )
    }
  }
  do.call(rbind, result_list)
}

# Run in batches with safe progress bar (now using mclapply for parallelism)
match_in_batches <- function(barseq_df, mapseq_split, mapseq_vbcs, max_dist = 1, batch_size = 100) {
  total <- nrow(barseq_df)
  num_batches <- ceiling(total / batch_size)
  pb <- txtProgressBar(min = 0, max = num_batches, style = 3)
  
  # Split into batches
  batches <- list()
  for (i in seq_len(num_batches)) {
    batch_idx <- ((i - 1) * batch_size + 1):min(i * batch_size, total)
    batches[[i]] <- barseq_df[batch_idx, ]
  }
  
  # Use mclapply for parallel processing (forking, no cluster needed)
  all_results <- mclapply(batches, function(batch) {
    process_batch(batch, mapseq_split, mapseq_vbcs, max_dist)
  }, mc.cores = detectCores() - 1)
  
  # Update progress bar (approximate, since mclapply doesn't support direct progress)
  setTxtProgressBar(pb, num_batches)
  close(pb)
  
  # Combine results
  valid_results <- all_results[!sapply(all_results, is.null)]
  if (length(valid_results) == 0) {
    return(data.frame(CellID = character(), vbc_read = character(), dist = integer()))
  }
  
  do.call(rbind, valid_results)
}

# --- Resolve matches: keep unique best hit; drop ties (including MAPseq 15-mer collisions) ---
resolve_best_unique <- function(matches, mapseq_df) {
  # matches: data.frame with columns CellID, vbc_read (15-mer), dist
  # mapseq_df: original MAPseq table with full vbc_read (>=15 nt)
  
  # 1) quantify MAPseq 15-mer collisions (how many MAPseq rows share each 15-mer)
  mapseq_short_counts <- mapseq_df %>%
    dplyr::mutate(VBC_short = substr(vbc_read, 1, 15)) %>%
    dplyr::count(VBC_short, name = "n_mapseq_rows")
  
  # 2) annotate each match with MAPseq multiplicity for that 15-mer
  m2 <- matches %>%
    dplyr::left_join(
      mapseq_short_counts,
      by = c("vbc_read" = "VBC_short")
    ) %>%
    dplyr::mutate(
      n_mapseq_rows = dplyr::coalesce(n_mapseq_rows, 0L)
    )
  
  if (nrow(m2) == 0) {
    return(list(
      selected = m2,
      dropped_ties = m2,
      summary = tibble::tibble(
        n_input_rows = 0L,
        n_cells_with_any_match = 0L,
        n_cells_selected = 0L,
        n_cells_dropped_ties = 0L
      )
    ))
  }
  
  # 3) for each CellID, keep only the minimum distance rows
  m_best <- m2 %>%
    dplyr::group_by(CellID) %>%
    dplyr::filter(dist == min(dist, na.rm = TRUE)) %>%
    dplyr::ungroup()
  
  # 4) define "unique winner":
  #    - exactly 1 best row for the cell
  #    - AND that 15-mer is unique in MAPseq (n_mapseq_rows == 1)
  per_cell <- m_best %>%
    dplyr::group_by(CellID) %>%
    dplyr::summarise(
      n_best_rows = dplyr::n(),
      n_mapseq_rows_winner = dplyr::first(n_mapseq_rows),
      dist_winner = dplyr::first(dist),
      .groups = "drop"
    )
  
  selected_cells <- per_cell %>%
    dplyr::filter(n_best_rows == 1, n_mapseq_rows_winner == 1) %>%
    dplyr::pull(CellID)
  
  selected <- m_best %>%
    dplyr::filter(CellID %in% selected_cells) %>%
    dplyr::select(CellID, vbc_read, dist)
  
  dropped_ties <- m_best %>%
    dplyr::filter(!CellID %in% selected_cells) %>%
    dplyr::arrange(CellID, dist, vbc_read)
  
  summary <- tibble::tibble(
    n_input_rows = nrow(matches),
    n_cells_with_any_match = dplyr::n_distinct(matches$CellID),
    n_cells_selected = dplyr::n_distinct(selected$CellID),
    n_cells_dropped_ties = dplyr::n_distinct(dropped_ties$CellID)
  )
  
  list(selected = selected, dropped_ties = dropped_ties, summary = summary)
}

# Run matching
matches_0 <- match_in_batches(barseq, mapseq_split, mapseq_vbcs, max_dist = 0)
matches_1 <- match_in_batches(barseq, mapseq_split, mapseq_vbcs, max_dist = 1)
matches_2 <- match_in_batches(barseq, mapseq_split, mapseq_vbcs, max_dist = 2)
matches_3 <- match_in_batches(barseq, mapseq_split, mapseq_vbcs, max_dist = 3)

# Resolve (keep unique best hit; drop ties/collisions)
res_0 <- resolve_best_unique(matches_0, mapseq)
res_1 <- resolve_best_unique(matches_1, mapseq)
res_2 <- resolve_best_unique(matches_2, mapseq)
res_3 <- resolve_best_unique(matches_3, mapseq)

print(res_0$summary)
print(res_1$summary)
print(res_2$summary)
print(res_3$summary)

# Optional: save dropped ties for auditing
readr::write_csv(res_0$dropped_ties, "matches_0_dropped_ties.csv")
readr::write_csv(res_1$dropped_ties, "matches_1_dropped_ties.csv")
readr::write_csv(res_2$dropped_ties, "matches_2_dropped_ties.csv")
readr::write_csv(res_3$dropped_ties, "matches_3_dropped_ties.csv")

summary_table_resolved <- tibble::tibble(
  Mismatches = 0:3,
  Cells_with_any_match = c(
    dplyr::n_distinct(matches_0$CellID),
    dplyr::n_distinct(matches_1$CellID),
    dplyr::n_distinct(matches_2$CellID),
    dplyr::n_distinct(matches_3$CellID)
  ),
  Cells_selected_unique_best = c(
    dplyr::n_distinct(res_0$selected$CellID),
    dplyr::n_distinct(res_1$selected$CellID),
    dplyr::n_distinct(res_2$selected$CellID),
    dplyr::n_distinct(res_3$selected$CellID)
  ),
  Cells_dropped_as_ties = c(
    dplyr::n_distinct(res_0$dropped_ties$CellID),
    dplyr::n_distinct(res_1$dropped_ties$CellID),
    dplyr::n_distinct(res_2$dropped_ties$CellID),
    dplyr::n_distinct(res_3$dropped_ties$CellID)
  )
)

print(summary_table_resolved)
readr::write_csv(summary_table_resolved, "summary_table_resolved_unique_best.csv")

############################################################################################################################################################################################################
# Subset projection data
subset_mapseq_by_matches <- function(matches, mapseq_df) {
  # Add a 15-base truncated barcode column to MapSeq (matching BarSeq)
  mapseq_df$VBC_short <- substr(mapseq_df$vbc_read, 1, 15)
  
  # Perform merge to retain full MapSeq info and bring in CellID/dist from matches
  merged <- merge(
    x = mapseq_df,
    y = matches,
    by.x = "VBC_short",
    by.y = "vbc_read"
  )
  
  # Optional: remove helper column
  merged$VBC_short <- NULL
  
  # Reorder columns for clarity (CellID, dist, original mapseq columns)
  merged <- merged[, c("CellID", "dist", setdiff(names(merged), c("CellID", "dist")))]
  
  return(merged)
}

check_duplicates <- function(proj) {
  duplicated_cellIDs <- proj$CellID[duplicated(proj$CellID)]
  unique_duplicated_cellIDs <- unique(duplicated_cellIDs)
  n_dups <- length(unique_duplicated_cellIDs)
  
  cat("Number of CellIDs with duplicates:", n_dups, "\n")
  cat("Unique duplicated CellIDs:\n")
  print(unique_duplicated_cellIDs)
  
  duplicated_rows <- proj[proj$CellID %in% unique_duplicated_cellIDs, ]
  cat("Rows with duplicated CellIDs:\n")
  print(head(duplicated_rows))  # print first few rows to avoid flooding the console
  
  table_counts_cellID <- table(proj$CellID)
  counts_more_than_one <- table_counts_cellID[table_counts_cellID > 1]
  cat("Counts of CellIDs appearing more than once:\n")
  print(counts_more_than_one)
  
  invisible(list(
    n_dups = n_dups,
    unique_duplicated_cellIDs = unique_duplicated_cellIDs,
    duplicated_rows = duplicated_rows,
    counts_more_than_one = counts_more_than_one
  ))
}

proj_0 <- subset_mapseq_by_matches(res_0$selected, mapseq)
proj_0_nodup <- distinct(proj_0)
write_csv(proj_0_nodup, "MapSeq_matched_projections_exact.csv")
result_0 <- check_duplicates(proj_0_nodup)

proj_1 <- subset_mapseq_by_matches(res_1$selected, mapseq)
proj_1_nodup <- distinct(proj_1)
write_csv(proj_1_nodup, "MapSeq_matched_projections_1_mismatch.csv")
result_1 <- check_duplicates(proj_1_nodup)

proj_2 <- subset_mapseq_by_matches(res_2$selected, mapseq)
proj_2_nodup <- distinct(proj_2)
write_csv(proj_2_nodup, "MapSeq_matched_projections_2_mismatch.csv")
result_2 <- check_duplicates(proj_2_nodup)

proj_3 <- subset_mapseq_by_matches(res_3$selected, mapseq)
proj_3_nodup <- distinct(proj_3)
write_csv(proj_3_nodup, "MapSeq_matched_projections_3_mismatch.csv")
result_3 <- check_duplicates(proj_3_nodup)

print_and_save_duplicates <- function(result, proj_nodup, filename) {
  dups <- result$unique_duplicated_cellIDs
  for (cell_id in dups) {
    indices <- which(proj_nodup$CellID == cell_id)
    cat("\nCellID:", cell_id, "\n")
    cat("Indices:", indices, "\n")
    cat("Rows:\n")
    print(proj_nodup[indices, , drop = FALSE])
  }
  write.csv(result$duplicated_rows, file = filename, row.names = FALSE)
}

print_and_save_duplicates(result_0, proj_0_nodup, "proj_0_duplicated_rows.csv")
print_and_save_duplicates(result_1, proj_1_nodup, "proj_1_duplicated_rows.csv")
print_and_save_duplicates(result_2, proj_2_nodup, "proj_2_duplicated_rows.csv")

#####################################################################################################################################################################################
# Only retain cells for genes vs projections processing which pass visual QC
# Load visual QC info CSV 
visualQC <- read.csv("/data/BARseq_soma_barcodes/LC_visualQC_barcoded_cells_780345.csv", header = TRUE, stringsAsFactors = FALSE)
# Check that the 'uid' column exists
if (!"uid" %in% colnames(visualQC)) {
  stop("No 'uid' column found in LC_visualQC_barcoded_cells.csv")
}

# Get the list of uids where good_barcoded is TRUE
good_uids <- visualQC$uid[visualQC$good_barcoded == 1]
length(good_uids)
# Subset cells for analyses involving gene expression profiles plus projections
proj_0_GENES <- proj_0_nodup[proj_0_nodup$CellID %in% good_uids, ]
dim(proj_0_GENES)
result_0_GENES <- check_duplicates(proj_0_GENES)

proj_1_GENES <- proj_1_nodup[proj_1_nodup$CellID %in% good_uids, ]
dim(proj_1_GENES)
result_1_GENES <- check_duplicates(proj_1_GENES)

proj_2_GENES <- proj_2_nodup[proj_2_nodup$CellID %in% good_uids, ]
dim(proj_2_GENES)
result_2_GENES <- check_duplicates(proj_2_GENES)

# Save
write_csv(proj_0_GENES, "MapSeq_matched_projections_exact_GENES.csv")
write_csv(proj_1_GENES, "MapSeq_matched_projections_1_mismatch_GENES.csv")
write_csv(proj_2_GENES, "MapSeq_matched_projections_2_mismatch_GENES.csv")


#####################################################################################################################################################################################
# FPR calculation to select optimal Hamming distance for matching (UPDATED for resolve_best_unique)
# Estimates probability that a random (or shuffled) BARseq barcode would be *selected* as a unique best MAPseq match
# under the same resolver you now use (min dist; drop ties; drop MAPseq 15-mer collisions).

# ---------------- SETTINGS ----------------
barcode_length <- 15
n_simulations <- 20  # Number of random simulations
max_mismatches <- 3
# ------------------------------------------

# Step 1: Generate Random Barcodes
generate_random_barcodes <- function(n, length = 15) {
  replicate(n, paste0(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = ""))
}

# Step 2: Run Random Matching Simulation (apply resolver!)
run_random_matching_simulation <- function(n_simulations, n_barseq, mapseq_split, mapseq_vbcs, mapseq_df, max_mismatches) {
  random_results <- vector("list", n_simulations)
  
  for (sim in seq_len(n_simulations)) {
    cat("Running simulation", sim, "of", n_simulations, "\n")
    
    # Generate random barcodes
    random_barcodes <- generate_random_barcodes(n_barseq, barcode_length)
    random_df <- data.frame(CellID = paste0("Random_", seq_len(n_barseq)), VBC = random_barcodes)
    
    # Run matching for each mismatch threshold
    matches_0 <- match_in_batches(random_df, mapseq_split, mapseq_vbcs, max_dist = 0)
    matches_1 <- match_in_batches(random_df, mapseq_split, mapseq_vbcs, max_dist = 1)
    matches_2 <- match_in_batches(random_df, mapseq_split, mapseq_vbcs, max_dist = 2)
    matches_3 <- match_in_batches(random_df, mapseq_split, mapseq_vbcs, max_dist = 3)
    
    # Resolve (keep unique best hit; drop ties/collisions) — SAME RULE AS REAL PIPELINE
    res_0 <- resolve_best_unique(matches_0, mapseq_df)
    res_1 <- resolve_best_unique(matches_1, mapseq_df)
    res_2 <- resolve_best_unique(matches_2, mapseq_df)
    res_3 <- resolve_best_unique(matches_3, mapseq_df)
    
    # Store results as "accepted/selected cells" (and ties for auditing)
    random_results[[sim]] <- list(
      selected_0 = dplyr::n_distinct(res_0$selected$CellID),
      selected_1 = dplyr::n_distinct(res_1$selected$CellID),
      selected_2 = dplyr::n_distinct(res_2$selected$CellID),
      selected_3 = dplyr::n_distinct(res_3$selected$CellID),
      ties_0 = dplyr::n_distinct(res_0$dropped_ties$CellID),
      ties_1 = dplyr::n_distinct(res_1$dropped_ties$CellID),
      ties_2 = dplyr::n_distinct(res_2$dropped_ties$CellID),
      ties_3 = dplyr::n_distinct(res_3$dropped_ties$CellID)
    )
  }
  
  return(random_results)
}

# Step 3: Summarize Random Selected/Tie Counts
summarize_random_results <- function(sim_results) {
  selected_means <- c(
    mean(sapply(sim_results, function(x) x$selected_0)),
    mean(sapply(sim_results, function(x) x$selected_1)),
    mean(sapply(sim_results, function(x) x$selected_2)),
    mean(sapply(sim_results, function(x) x$selected_3))
  )
  ties_means <- c(
    mean(sapply(sim_results, function(x) x$ties_0)),
    mean(sapply(sim_results, function(x) x$ties_1)),
    mean(sapply(sim_results, function(x) x$ties_2)),
    mean(sapply(sim_results, function(x) x$ties_3))
  )
  list(selected_means = selected_means, ties_means = ties_means)
}

# Step 4: Calculate Error Rates (FPR) for selected matches
calculate_error_rates <- function(mean_selected, n_barseq) {
  mean_selected / n_barseq
}

# ---------------- RUN (RANDOM) ----------------

# Real BarSeq Data
n_barseq <- nrow(barseq)

# IMPORTANT: Real counts should reflect what you now *accept* (res_*), not raw matches_*
real_selected_counts <- c(
  dplyr::n_distinct(res_0$selected$CellID),
  dplyr::n_distinct(res_1$selected$CellID),
  dplyr::n_distinct(res_2$selected$CellID),
  dplyr::n_distinct(res_3$selected$CellID)
)

real_ties_counts <- c(
  dplyr::n_distinct(res_0$dropped_ties$CellID),
  dplyr::n_distinct(res_1$dropped_ties$CellID),
  dplyr::n_distinct(res_2$dropped_ties$CellID),
  dplyr::n_distinct(res_3$dropped_ties$CellID)
)

# Random Simulations
random_sim_file <- "random_simulation_results_resolved.rds"
if (file.exists(random_sim_file)) {
  random_results <- readRDS(random_sim_file)
} else {
  random_results <- run_random_matching_simulation(
    n_simulations, n_barseq, mapseq_split, mapseq_vbcs, mapseq, max_mismatches
  )
  saveRDS(random_results, file = random_sim_file)
}

# Summarize Random Results
random_summary <- summarize_random_results(random_results)
random_mean_selected <- random_summary$selected_means
random_mean_ties <- random_summary$ties_means

# Calculate Error Rates (FPR for being selected)
error_rates_selected <- calculate_error_rates(random_mean_selected, n_barseq)

# Construct comparison_df
comparison_df <- data.frame(
  Mismatches = 0:max_mismatches,
  Real_Selected = real_selected_counts,
  Real_Dropped_Ties = real_ties_counts,
  Random_Mean_Selected = round(random_mean_selected),
  Random_Mean_Dropped_Ties = round(random_mean_ties),
  FPR_Selected = round(error_rates_selected, 4)
)
comparison_df$Avg_Random_Selected_Per_BarSeq <- round(random_mean_selected / n_barseq, 4)

# Print and save
print(comparison_df)
write.csv(comparison_df, file = "comparison_df_FPR_error_resolved.csv", row.names = FALSE)

# Visualize Results (RANDOM)
p <- ggplot(comparison_df, aes(x = Mismatches, y = Avg_Random_Selected_Per_BarSeq)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "Random *Selected* Match Rate per BarSeq Cell (post-resolver)",
    x = "Allowed Mismatches",
    y = "Avg Random Selected per BarSeq"
  ) +
  theme_minimal()
print(p)
ggsave("BarSeq-MapSeq_FPR_resolved.pdf", plot = p, device = "pdf", width = 6, height = 6)

#####################################################################################################################################################################################
# Shuffled BarSeq Data Simulation for FPR estimation
# Function to shuffle barcode keeping 9th nucleotide fixed
shuffle_barcode <- function(vbc) {
  chars <- strsplit(vbc, "")[[1]]
  fixed <- chars[9]
  others <- chars[-9]
  shuffled_others <- sample(others)
  new_chars <- chars
  new_chars[-9] <- shuffled_others
  paste0(new_chars, collapse = "")
}

# Run Shuffled Matching Simulation (apply resolver!)
run_shuffled_matching_simulation <- function(n_simulations, barseq_df, mapseq_split, mapseq_vbcs, mapseq_df, max_mismatches) {
  shuffled_results <- vector("list", n_simulations)
  
  for (sim in seq_len(n_simulations)) {
    cat("Running shuffled simulation", sim, "of", n_simulations, "\n")
    
    # Shuffle barcodes
    shuffled_barcodes <- sapply(barseq_df$VBC, shuffle_barcode)
    shuffled_df <- data.frame(CellID = barseq_df$CellID, VBC = shuffled_barcodes)
    
    # Run matching for each mismatch threshold
    matches_0 <- match_in_batches(shuffled_df, mapseq_split, mapseq_vbcs, max_dist = 0)
    matches_1 <- match_in_batches(shuffled_df, mapseq_split, mapseq_vbcs, max_dist = 1)
    matches_2 <- match_in_batches(shuffled_df, mapseq_split, mapseq_vbcs, max_dist = 2)
    matches_3 <- match_in_batches(shuffled_df, mapseq_split, mapseq_vbcs, max_dist = 3)
    
    # Resolve (same rule as real pipeline)
    res_0 <- resolve_best_unique(matches_0, mapseq_df)
    res_1 <- resolve_best_unique(matches_1, mapseq_df)
    res_2 <- resolve_best_unique(matches_2, mapseq_df)
    res_3 <- resolve_best_unique(matches_3, mapseq_df)
    
    shuffled_results[[sim]] <- list(
      selected_0 = dplyr::n_distinct(res_0$selected$CellID),
      selected_1 = dplyr::n_distinct(res_1$selected$CellID),
      selected_2 = dplyr::n_distinct(res_2$selected$CellID),
      selected_3 = dplyr::n_distinct(res_3$selected$CellID),
      ties_0 = dplyr::n_distinct(res_0$dropped_ties$CellID),
      ties_1 = dplyr::n_distinct(res_1$dropped_ties$CellID),
      ties_2 = dplyr::n_distinct(res_2$dropped_ties$CellID),
      ties_3 = dplyr::n_distinct(res_3$dropped_ties$CellID)
    )
  }
  
  return(shuffled_results)
}

# Run shuffled simulations (use a NEW filename; old file was pre-resolver format)
shuf_sim_file <- "shuffled_simulation_results_resolved.rds"
if (file.exists(shuf_sim_file)) {
  shuffled_results <- readRDS(shuf_sim_file)
} else {
  shuffled_results <- run_shuffled_matching_simulation(
    n_simulations, barseq, mapseq_split, mapseq_vbcs, mapseq, max_mismatches
  )
  saveRDS(shuffled_results, file = shuf_sim_file)
}

# Summarize Shuffled Results
shuffled_summary <- summarize_random_results(shuffled_results)
shuffled_mean_selected <- shuffled_summary$selected_means
shuffled_mean_ties <- shuffled_summary$ties_means

# Calculate Error Rates for shuffled
shuffled_error_rates_selected <- calculate_error_rates(shuffled_mean_selected, n_barseq)

# Construct shuffled_comparison_df
shuffled_comparison_df <- data.frame(
  Mismatches = 0:max_mismatches,
  Real_Selected = real_selected_counts,
  Real_Dropped_Ties = real_ties_counts,
  Shuffled_Mean_Selected = round(shuffled_mean_selected),
  Shuffled_Mean_Dropped_Ties = round(shuffled_mean_ties),
  Shuffled_FPR_Selected = round(shuffled_error_rates_selected, 4)
)
shuffled_comparison_df$Avg_Shuffled_Selected_Per_BarSeq <- round(shuffled_mean_selected / n_barseq, 4)

# Print and save
print(shuffled_comparison_df)
write.csv(shuffled_comparison_df, file = "shuffled_comparison_df_FPR_error_resolved.csv", row.names = FALSE)

# Visualize Results (SHUFFLED)
p_shuf <- ggplot(shuffled_comparison_df, aes(x = Mismatches, y = Avg_Shuffled_Selected_Per_BarSeq)) +
  geom_line(color = "red", linewidth = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Shuffled *Selected* Match Rate per BarSeq Cell (post-resolver)",
    x = "Allowed Mismatches",
    y = "Avg Shuffled Selected per BarSeq"
  ) +
  theme_minimal()
print(p_shuf)
ggsave("BarSeq-MapSeq_Shuffled_FPR_resolved.pdf", plot = p_shuf, device = "pdf", width = 6, height = 6)
