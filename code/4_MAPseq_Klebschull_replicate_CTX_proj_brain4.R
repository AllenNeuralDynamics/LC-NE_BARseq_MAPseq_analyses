################# MAPseq Barcode checks for CTX slices replicating Klebschull paper ##################

#set working directory 
setwd('/scratch/BARseq_780346/')

#####################################################################################################################################################################################
#read in files which contain relevant information for projections such as projection matrix, cell IDs and cluster types
projection_matrix <- readr::read_csv("./MapSeq_matched_projections_1_mismatch.csv", col_names = TRUE)

head(projection_matrix)
LCNEneurons <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
colData(LCNEneurons)
cell_metadata <- as.data.frame(colData(LCNEneurons))

proj_index <- read_tsv("/data/780346_2025-06-11_00-00-00/MAPseq/M305_20251030_USEthis/M305sampleinfo.tsv")
head(proj_index)
sampleinfo <- read_excel("./M305_sample_information.xlsx", sheet = "Sample information", skip = 1,range = "A2:J122")
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

# Merge cell metadata (cell type, CCF info) into projection_matrix
# Select only the desired columns from cell_meta
cell_metadata_subset <- cell_metadata %>%
  select(uid, slice, barcode, CCF_AP, CCF_DV, CCF_ML, louvain_cluster)
# CellID in projection_matrix matches uid in cell_meta to join with projection_matrix and covert it to data frame
df <- projection_matrix %>%
  left_join(cell_metadata_subset, by = c("CellID" = "uid")) %>%
  as.data.frame()
# Add information about the hemispehre
df$inRH <- ifelse(df$CCF_ML > 228, 1, 0)
head(df)
print(colnames(df))

# Use cellID and cluster to annotate rows
df$CellID <- as.character(df$CellID)
df$louvain_cluster <- as.character(df$louvain_cluster)
df$row_id <- paste(df$CellID, df$louvain_cluster, sep=".")
rownames(df) <- make.unique(df$row_id)
head(df)

#check for non-unique cells based on CellID.louvain_cluster.unique pattern and count how many there are
rn <- rownames(df)
dot_counts <- lengths(regmatches(rn, gregexpr("\\.", rn)))
table(dot_counts)

id_col    <- "CellID"
drop_cols <- c("dist", "vbc_read", "louvain_cluster", "row_id", "barcode")

meta_key <- df %>%
  tibble::rownames_to_column("rowname") %>%
  dplyr::select(rowname, CellID, slice, barcode, vbc_read, louvain_cluster, CCF_AP, CCF_DV, CCF_ML, inRH)

df_dropped <- df %>% dplyr::select(all_of(c(id_col, drop_cols)))
df_dropped <- df_dropped %>% tibble::rownames_to_column("rowname")
df         <- df %>% dplyr::select(-all_of(c(id_col, drop_cols)))

print(colnames(df))

###################################################### barcode counts sanity check ###############################################################
# Plot histogram of raw value distributions for all columns to look at general barcode expression levels
# Specify metadata columns to exclude
exclude_cols <- c("slice", "CCF_AP", "CCF_DV", "CCF_ML", "inRH")
# Select only projection columns for plotting
proj_df <- df[, !(colnames(df) %in% exclude_cols)]

# Reshape to long format for plotting, account for possible NA and Inf values by setting those to zero
df_long <- proj_df %>%
  pivot_longer(everything(), names_to = "ROI", values_to = "value") %>%
  mutate(value = replace(value, is.na(value) | is.infinite(value), 0))
range(df_long$value)
p <- ggplot(df_long, aes(x = value))
p <- p + geom_histogram(binwidth = 1, fill = "blue", color = "blue", linewidth = 0.6) +
  labs(x = "Detected barcode counts", y = "Number of cells", title = "Raw barcode counts") +  theme_minimal()
p <- p + xlim(c(0, 15000))
p <- p + ylim(c(0, 100))
print(p)
ggsave("raw_barcode_counts.pdf", plot = p, device = "pdf", width = 12, height = 8)  

p <- p + xlim(c(0, 1500))
print(p)
ggsave("raw_barcode_counts_zoomin.pdf", plot = p, device = "pdf", width = 12, height =8)

#####################################################################################################################################################################################
# Organize rois, roi names and sample index
# Get all column names
col_names <- colnames(df)
# Helper to extract hemisphere and base ROI
extract_info <- function(name) {
  hemi <- ifelse(str_detect(name, "_RH"), "RH",
                 ifelse(str_detect(name, "_LH"), "LH", NA))
  # base ROI = stuff before _RH/_LH and the dot
  base_roi <- str_replace(name, "(_RH|_LH)\\..*$", "")
  
  # slide = number after the dot, e.g. LC_RH.5 -> 5
  slide <- as.integer(str_extract(name, "(?<=\\.)\\d+$"))
  return(data.frame(colname = name, base_roi = base_roi, hemisphere = hemi, slide = slide, stringsAsFactors = FALSE))
}
# Build metadata for all columns with _RH or _LH
meta <- do.call(rbind, lapply(col_names, extract_info)) %>%
  filter(!is.na(hemisphere))
# Original column index in df
meta$tube_idx <- match(meta$colname, col_names)

# Within each ROI/hemisphere, order by slide (then tube_idx for ties)
meta <- meta %>%
  group_by(hemisphere, base_roi) %>%
  arrange(slide, tube_idx, .by_group = TRUE) %>%
  ungroup()

# Anchor each ROI by its earliest slide, not earliest tube_idx
meta <- meta %>%
  group_by(hemisphere, base_roi) %>%
  mutate(min_slide = min(slide, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(hemisphere, min_slide, slide, tube_idx) %>%  # use slide ordering
  select(-min_slide)

# Assign plotting index within each hemisphere
meta <- meta %>%
  group_by(hemisphere) %>%
  mutate(idx = row_number()) %>%
  ungroup()

# Split RH and LH in original df order for convenience
RH_meta <- meta %>% filter(hemisphere == "RH") %>% arrange(tube_idx)
RH_tubes <- RH_meta$tube_idx
RH_rois  <- RH_meta$base_roi
RH_idx   <- RH_meta$idx

LH_meta <- meta %>% filter(hemisphere == "LH") %>% arrange(tube_idx)
LH_tubes <- LH_meta$tube_idx
LH_rois  <- LH_meta$base_roi
LH_idx   <- LH_meta$idx

# Sanity check ordering is correct and is tied to A-P slide order
as.data.frame(meta %>%
  group_by(hemisphere, base_roi) %>%
  summarise(first_slide = min(slide), .groups = "drop") %>%
  arrange(hemisphere, first_slide, base_roi))
nrow(meta)  # should be 106 since 3 spinal cord _SP samples are not considered here

############################################# identify high expressor cells computationally ##########################################
############################################# ipsilateral projections ################################################################
# Find cells with more than 75 barcode counts in any given cortical roi tube
# Get ctx column indices for RH and LH separately
ctx_RH_indices <- grep("ctx.*_RH", names(df))
#ctx_RH_indices <- ctx_RH_indices[-length(ctx_RH_indices)] #ctx is unique to RH, not present in LH
ctx_LH_indices <- grep("ctx.*_LH", names(df))

# --- RH ctx columns ---
df_RH <- df[df$inRH == 1, ]
df_subset_RH <- df_RH[, ctx_RH_indices, drop = FALSE]
df_subset_RH[] <- lapply(df_subset_RH, as.numeric) # ensure numeric
row_exceeds_RH <- apply(df_subset_RH, 1, function(x) any(x > 75))
df_subset_RH <- df_subset_RH[row_exceeds_RH, , drop = FALSE]
df_subset_RH <- df_subset_RH[rowSums(df_subset_RH) > 0, , drop = FALSE] # safety net

# Order columns to be consistent with your custom idx
RH_ctx_meta <- RH_meta %>%
  filter(str_detect(base_roi, "ctx")) %>%
  arrange(slide, base_roi, tube_idx)   # slide-major order
# Only keep columns that exist in df_subset_RH (in case of mismatch)
RH_ctx_cols <- RH_ctx_meta$colname[RH_ctx_meta$colname %in% colnames(df_subset_RH)]
df_subset_RH <- df_subset_RH[, RH_ctx_cols]
dim(df_subset_RH)

# --- LH ctx columns ---
df_LH <- df[df$inRH == 0, ]
df_subset_LH <- df_LH[, ctx_LH_indices, drop = FALSE]
df_subset_LH[] <- lapply(df_subset_LH, as.numeric) # ensure numeric
row_exceeds_LH <- apply(df_subset_LH, 1, function(x) any(x > 75))
df_subset_LH <- df_subset_LH[row_exceeds_LH, , drop = FALSE]
df_subset_LH <- df_subset_LH[rowSums(df_subset_LH) > 0, , drop = FALSE] # safety net

LH_ctx_meta <- LH_meta %>%
  filter(str_detect(base_roi, "ctx")) %>%
  arrange(slide, base_roi, tube_idx)
LH_ctx_cols <- LH_ctx_meta$colname[LH_ctx_meta$colname %in% colnames(df_subset_LH)]
df_subset_LH <- df_subset_LH[, LH_ctx_cols]
dim(df_subset_LH)

# Pivot to long format for plotting
df_subset_longRH <- df_subset_RH %>% 
  rownames_to_column(var = "Cell") %>% 
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         Source = "RH")
df_subset_longRH$ROI <- factor(df_subset_longRH$ROI, levels = RH_ctx_cols)

df_subset_longLH <- df_subset_LH %>% 
  rownames_to_column(var = "Cell") %>% 
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         Source = "LH")
df_subset_longLH$ROI <- factor(df_subset_longLH$ROI, levels = LH_ctx_cols)

plot_cells_in_batches <- function(df_long, batch_size = 16, save = FALSE, prefix = "plot_batch", title_prefix = "Cells") {
  unique_cells <- unique(df_long$Cell)
  n_cells <- length(unique_cells)
  n_batches <- ceiling(n_cells / batch_size)
  
  # Ensure ROI is a factor with the correct order
  if (!is.factor(df_long$ROI)) {
    df_long$ROI <- factor(df_long$ROI, levels = unique(df_long$ROI))
  }
  
  for (i in seq_len(n_batches)) {
    batch_cells <- unique_cells[((i - 1) * batch_size + 1):min(i * batch_size, n_cells)]
    batch_data <- df_long %>% filter(Cell %in% batch_cells)
    
    p <- ggplot(batch_data, aes(x = ROI, y = Value, group = Cell)) +
      geom_line() + 
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
      facet_wrap(~ Cell, scales = "free_y") + 
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            axis.title.x = element_text(size = 8),
            axis.title.y = element_text(size = 8),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()) +
      ylab('projection strength (barcode count)') +
      xlab('ctx projection target (ROIs)') +
      ggtitle(paste(title_prefix, "(Cells", min(batch_cells), "to", max(batch_cells), ")"))
    
    print(p)
    
    if (save) {
      ggsave(filename = sprintf("%s_%02d.pdf", prefix, i), plot = p, dpi = 300, width = 12, height = 8)
    }
  }
}

# Usage:
plot_cells_in_batches(df_subset_longRH, batch_size = 30, save = TRUE, prefix = "RH_ctx_cells", title_prefix = "RH somas projecting to RH ctx")
plot_cells_in_batches(df_subset_longLH, batch_size = 30, save = TRUE, prefix = "LH_ctx_cells", title_prefix = "LH somas projecting to LH ctx")

########################################################### sum across ctx sections to mimic Justus's slice data and Mathew's ExA-SPIM flow ##########################################################
head(df)
dim(df)
#Construct data table where all cortical ROIs summed by slide=slice for direct comparison to ExA-SPIM data
# Keep only cortical ROIs (anything with 'ctx' in the base ROI name)
ctx_meta <- meta %>%
  filter(str_detect(base_roi, "ctx"))
ctx_slides <- sort(unique(ctx_meta$slide))
# Build a ctx_df with rows = cells, cols = "slice<slide>"
ctx_df <- data.frame(row.names = rownames(df))

for (s in ctx_slides) {
  # all ctx columns on this slide (both hemispheres)
  cols <- ctx_meta$colname[ctx_meta$slide == s]
  cols_present <- cols[cols %in% colnames(df)]  # safety
  
  ctx_df[[paste0("slice", s)]] <- rowSums(df[, cols_present, drop = FALSE],
                                          na.rm = TRUE)
}
# # Define the ROI groups - manual version as for brain3 implementation
# roi_groups <- list(
#   slice1 = c("motor ctx_RH.3", "orb ctx_RH.3", "motor ctx_LH.3", "orb ctx_LH.3" ),  
#   slice2 = c("motor ctx_RH.4", "orb ctx_RH.4", "motor ctx_LH.4", "orb ctx_LH.4"),
#   slice3 = c("ctx 1_RH.5", "ctx 2_RH.5", "ctx 3_RH.5", "ctx 1_LH.5", "ctx 2_LH.5", "ctx 3_LH.5"),
#   slice4 = c("ctx 1_RH.6", "ctx 2_RH.6", "ctx 3_RH.6", "ctx 1_LH.6", "ctx 2_LH.6", "ctx 3_LH.6"),
#   slice5 = c("ctx 1_RH.7", "ctx 2_RH.7", "ctx 3_RH.7", "ctx 1_LH.7", "ctx 2_LH.7", "ctx 3_LH.7"),
#   slice6 = c("ctx 1_RH.8", "ctx 2_RH.8", "ctx 3_RH.8", "ctx 1_LH.8", "ctx 2_LH.8", "ctx 3_LH.8"),
#   slice7 = c("ctx 1_RH.9", "ctx 2_RH.9", "ctx 3_RH.9", "ctx 1_LH.9", "ctx 2_LH.9", "ctx 3_LH.9" ),
#   slice8 = c("ctx 1_RH.10", "ctx 2_RH.10", "ctx 3_RH.10", "ctx 1_LH.10", "ctx 2_LH.10", "ctx 3_LH.10" ),
#   slice9 = c("ctx_RH.11", "ctx_LH.11")
# )
# # Create a new data frame with summed values for each slice
# ctx_df <- data.frame(row.names = rownames(df))
# for (slice in names(roi_groups)) {
#   cols <- roi_groups[[slice]]
#   # Only keep columns that exist in df (in case some are missing)
#   cols_present <- cols[cols %in% colnames(df)]
#   ctx_df[[slice]] <- rowSums(df[, cols_present, drop=FALSE], na.rm=TRUE)
# }
head(ctx_df)
dim(ctx_df)

# Drop cells (rows) where all slices are zero
ctx_df_nonzero <- ctx_df[rowSums(ctx_df) > 0, ]
# Add cell IDs as a column for plotting
ctx_df_nonzero <- ctx_df_nonzero %>% 
  tibble::rownames_to_column(var = "Cell")
# Pivot to long format for plotting
ctx_long <- ctx_df_nonzero %>%
  tidyr::pivot_longer(
    cols = starts_with("slice"),
    names_to = "ROI",
    values_to = "Value"
  ) %>%
  mutate(
    lower = Value - sqrt(Value),
    upper = Value + sqrt(Value)
  )
# Ensure slice order is correct
ctx_long$ROI <- factor(ctx_long$ROI, levels = colnames(ctx_df))

# Plot function (similar to before)
plot_cells_in_batches(ctx_long, batch_size = 30, save = TRUE, prefix = "ctx_slices", title_prefix = "Cells projecting to ctx slices")

#####################################################################################################################################################################################
# Identify cells with 'duplicate' ctx projection profiles and check whether those are true duplicates in the full dataset or apparent duplicates based on subset of ctx regions
find_duplicate_cells <- function(df_long) {
  # Pivot to wide format: one row per cell, columns = ROIs
  df_wide <- df_long %>%
    select(Cell, ROI, Value) %>%
    pivot_wider(names_from = ROI, values_from = Value, values_fill = 0) %>%
    arrange(Cell)
  
  # Create a hash/signature for each cell's projection profile
  df_wide$profile_hash <- apply(df_wide %>% select(-Cell), 1, function(x) digest(x, algo = "md5"))
  
  # Group by hash and list cells with identical profiles
  dup_groups <- df_wide %>%
    group_by(profile_hash) %>%
    summarise(cells = list(Cell), n = n()) %>%
    filter(n > 1)
  
  print(dup_groups)
  return(dup_groups)
}

# Example usage:
dup_RH <- find_duplicate_cells(df_subset_longRH)
dup_LH <- find_duplicate_cells(df_subset_longLH)
dup_ctx <- find_duplicate_cells(ctx_long)

# Extract the cell IDs for each duplicate group
# For RH duplicates
dup_RH_cells <- dup_RH %>% unnest(cells)
print(dup_RH_cells)
# For LH duplicates
dup_LH_cells <- dup_LH %>% unnest(cells)
print(dup_LH_cells)
# For ctx duplicates
dup_ctx_cells <- dup_ctx %>% unnest(cells)
print(dup_ctx_cells)

# Function to compare full profiles for each duplicate group
compare_full_profiles <- function(dup_table, df) {
  for (i in seq_len(nrow(dup_table))) {
    cell_group <- dup_table$cells[[i]]
    if (length(cell_group) > 1) {
      profiles <- df[cell_group, , drop = FALSE]
      identical_profiles <- all(apply(profiles, 2, function(col) length(unique(col)) == 1))
      cat("Group", i, "Cells:", paste(cell_group, collapse = ", "), 
          "Identical full profiles:", identical_profiles, "\n")
    }
  }
}

# This comparison takes into account metadata differences
compare_full_profiles(dup_ctx, df)
compare_full_profiles(dup_RH, df)
compare_full_profiles(dup_LH, df)

# This comparison only checks numeric entries, not metadata
exclude_cols <- c("slice", "CCF_AP", "CCF_DV", "CCF_ML", "inRH")
full_proj_df <- df[, !(colnames(df) %in% exclude_cols)]
compare_full_profiles(dup_ctx, full_proj_df)
compare_full_profiles(dup_RH, full_proj_df)
compare_full_profiles(dup_LH, full_proj_df)

#####################################################################################################################################################################################
#####################################################################################################################################################################################
# Identify duplicate cells across the entire projection table - aka cells which have different metadata, but identical MAPseq profiles across all regions
# MAPseq-profile duplicates
head(full_proj_df)
proj_matrix <- as.data.frame(full_proj_df)
profile_hash <- apply(proj_matrix, 1, function(x) digest(x, algo = "md5"))
hash_tab <- table(profile_hash)
# number of groups with at least 2 identical rows
n_duplicate_groups <- sum(hash_tab > 1)
# total number of *extra* duplicated rows (beyond 1 per group)
n_extra_duplicate_rows <- sum(hash_tab[hash_tab > 1] - 1)
n_duplicate_groups
n_extra_duplicate_rows
dup_hashes <- names(hash_tab)[hash_tab > 1]
dup_df <- data.frame(
  CellID = rownames(proj_matrix),
  profile_hash = profile_hash,
  stringsAsFactors = FALSE
)
dup_df <- dup_df[dup_df$profile_hash %in% dup_hashes, ]
# inspect groups
dup_groups <- split(dup_df$CellID, dup_df$profile_hash)
dup_groups_mapseq <- dup_groups
length(dup_groups_mapseq)
head(dup_groups_mapseq)

# Identify duplicate BARseq barcodes (vbc_read) across rows
find_duplicate_vbc_groups <- function(meta_key, vbc_col = "vbc_read", min_group_size = 2) {
  stopifnot("rowname" %in% names(meta_key))
  stopifnot(vbc_col %in% names(meta_key))
  
  vbc_groups <- meta_key %>%
    dplyr::filter(!is.na(.data[[vbc_col]]), .data[[vbc_col]] != "") %>%
    dplyr::group_by(.data[[vbc_col]]) %>%
    dplyr::summarise(cells = list(rowname), n = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n >= min_group_size)
  
  # Convert to named list like your existing dup_groups
  out <- vbc_groups$cells
  names(out) <- as.character(vbc_groups[[vbc_col]])
  out
}
dup_vbc_groups <- find_duplicate_vbc_groups(meta_key, vbc_col = "vbc_read")
length(dup_vbc_groups)
head(dup_vbc_groups)

#####################################################################################################################################################################################
# Build a cell-centric registry: ONE ROW PER CELL 
# (Use all cells as base; later filter to only duplicated cells for plotting.)
# helper: cell -> semicolon-separated group ids
group_membership_table <- function(group_list, group_colname) {
  tibble::tibble(
    Cell = unlist(group_list, use.names = FALSE),
    group = rep(names(group_list), times = lengths(group_list))
  ) %>%
    dplyr::group_by(Cell) %>%
    dplyr::summarise(!!group_colname := paste(unique(group), collapse = ";"), .groups = "drop")
}

mapseq_membership <- group_membership_table(dup_groups_mapseq, rlang::sym("MAPseq_group"))
vbc_membership    <- group_membership_table(dup_vbc_groups,    rlang::sym("vbc_read_group"))

dup_registry <- tibble::tibble(Cell = rownames(df)) %>%
  dplyr::left_join(mapseq_membership, by = "Cell") %>%
  dplyr::left_join(vbc_membership,    by = "Cell") %>%
  dplyr::mutate(
    MAPseq_profile_dup = !is.na(MAPseq_group) & MAPseq_group != "",
    BARseq_vbc_dup     = !is.na(vbc_read_group) & vbc_read_group != "",
    dup_class = dplyr::case_when(
      MAPseq_profile_dup & BARseq_vbc_dup ~ "Both",
      MAPseq_profile_dup ~ "MAPseq_only",
      BARseq_vbc_dup     ~ "BARseq_only",
      TRUE ~ "None"
    )
  )

# Add coords
coords <- df[, c("CCF_ML", "CCF_DV", "CCF_AP", "inRH"), drop = FALSE] %>%
  tibble::rownames_to_column("Cell")
dup_registry <- dplyr::left_join(dup_registry, coords, by = "Cell")

# Checks
table(dup_registry$MAPseq_profile_dup, useNA = "ifany")
table(dup_registry$BARseq_vbc_dup, useNA = "ifany")
table(dup_registry$dup_class)

# Create plotting df
dup_plot <- dup_registry %>%
  dplyr::filter(MAPseq_profile_dup | BARseq_vbc_dup)

#####################################################################################################################################################################################
# DUPLICATE GROUP DIAGNOSTICS + PLOTTING (single run; faceted; pair table collapsed with flags)
# Build mapping between VBC groups and MAPseq hash groups (same partition, different order)
sig <- function(group_list) {
  vapply(group_list, function(x) paste(sort(x), collapse = "|"), character(1))
}

mapseq_sig_to_id <- setNames(names(dup_groups_mapseq), sig(dup_groups_mapseq))
vbc_sigs <- sig(dup_vbc_groups)

vbc_to_mapseq <- setNames(mapseq_sig_to_id[vbc_sigs], names(dup_vbc_groups))
stopifnot(all(!is.na(vbc_to_mapseq)))  # must map 1:1

# Create plotting df (faceted style; label includes both IDs)
plot_df <- lapply(names(dup_vbc_groups), function(vbc_id) {
  cells <- dup_vbc_groups[[vbc_id]]
  df_subset <- df[cells, c("CCF_ML", "CCF_DV", "CCF_AP"), drop = FALSE]
  df_subset$CellID <- rownames(df_subset)
  df_subset$group_hash <- paste0("VBC_", vbc_id, "\nMAPseqHash_", vbc_to_mapseq[[vbc_id]])
  df_subset
}) %>% dplyr::bind_rows()

plot_df$group_hash <- factor(plot_df$group_hash)

plot_df <- plot_df %>%
  dplyr::group_by(group_hash) %>%
  dplyr::mutate(cell_order = dplyr::row_number(),
                color = ifelse(cell_order == 1, "green", "red")) %>%
  dplyr::ungroup()

p <- ggplot(plot_df, aes(x = CCF_ML, y = CCF_DV, label = CellID, color = color)) +
  geom_point(size = 2) +
  geom_text_repel(nudge_y = 15, size = 2.4, show.legend = FALSE) +
  scale_color_manual(values = c("green", "red")) +
  facet_wrap(~ group_hash) +
  theme_minimal() +
  labs(
    title = "Duplicate-cell groups in coronal plane (ML × DV)",
    x = "CCF ML (medial–lateral)",
    y = "CCF DV (dorsal–ventral)",
    color = "Cell Order"
  ) +
  coord_fixed() +
  theme(legend.position = "bottom")

print(p)
ggsave("duplicate_cells_coronal_plane_plot.pdf", p, width = 12, height = 10)

# Pairwise distance table computed ONCE per group, then collapsed to one row per pair with flags
euclidean_distance <- function(coord1, coord2) {
  sqrt(sum((coord1 - coord2)^2))
}

compute_group_distances <- function(cell_ids, df) {
  coords <- df[cell_ids, c("CCF_AP", "CCF_DV", "CCF_ML")]
  pairs <- combn(cell_ids, 2, simplify = FALSE)
  distances <- sapply(pairs, function(pair) {
    coord1 <- as.numeric(coords[pair[1], ])
    coord2 <- as.numeric(coords[pair[2], ])
    euclidean_distance(coord1, coord2)
  })
  pair_df <- do.call(rbind, lapply(pairs, function(p) data.frame(cell1 = p[1], cell2 = p[2])))
  pair_df$distance_voxels  <- distances
  pair_df$distance_microns <- distances * 25
  return(pair_df)
}

# Build uncollapsed pair rows with group_hash
all_distances <- lapply(names(dup_vbc_groups), function(vbc_id) {
  cells <- dup_vbc_groups[[vbc_id]]
  if (length(cells) > 1) {
    dist_df <- compute_group_distances(cells, df)
    dist_df$group_hash <- paste0("VBC_", vbc_id, "|MAPseqHash_", vbc_to_mapseq[[vbc_id]])
    return(dist_df)
  } else {
    return(NULL)
  }
})
distance_table_raw <- do.call(rbind, all_distances)

# Collapse: one row per unordered pair
canon_pair <- function(a, b) {
  ifelse(a < b, paste(a, b, sep = "|"), paste(b, a, sep = "|"))
}

distance_table <- distance_table_raw %>%
  dplyr::mutate(pair_id = canon_pair(cell1, cell2)) %>%
  dplyr::group_by(pair_id) %>%
  dplyr::summarise(
    cell1 = dplyr::first(ifelse(cell1 < cell2, cell1, cell2)),
    cell2 = dplyr::first(ifelse(cell1 < cell2, cell2, cell1)),
    distance_voxels  = dplyr::first(distance_voxels),
    distance_microns = dplyr::first(distance_microns),
    
    # since groups are the same partition, any pair here is both kinds; keep explicit flags anyway
    MAPseq_pair = TRUE,
    BARseq_pair = TRUE,
    
    # keep traceability to the (shared) group id(s)
    group_hash = paste(unique(group_hash), collapse = ";"),
    .groups = "drop"
  )

print(distance_table)
write.csv(distance_table, "duplicate_pairs_distance_table_collapsed.csv", row.names = FALSE)

# Sanity table (minimal change; keep group_hash like your original)
sanity_check <- function(distance_table, df) {
  check_df <- data.frame()
  for (i in seq_len(nrow(distance_table))) {
    cell1 <- distance_table$cell1[i]
    cell2 <- distance_table$cell2[i]
    coord1 <- df[cell1, c("CCF_AP", "CCF_DV", "CCF_ML")]
    coord2 <- df[cell2, c("CCF_AP", "CCF_DV", "CCF_ML")]
    
    temp_df <- data.frame(
      group_hash = distance_table$group_hash[i],
      cell1 = cell1,
      AP1 = coord1$CCF_AP, DV1 = coord1$CCF_DV, ML1 = coord1$CCF_ML,
      cell2 = cell2,
      AP2 = coord2$CCF_AP, DV2 = coord2$CCF_DV, ML2 = coord2$CCF_ML,
      distance_voxels  = distance_table$distance_voxels[i],
      distance_microns = distance_table$distance_microns[i],
      MAPseq_pair = distance_table$MAPseq_pair[i],
      BARseq_pair = distance_table$BARseq_pair[i],
      stringsAsFactors = FALSE
    )
    check_df <- rbind(check_df, temp_df)
  }
  return(check_df)
}

sanity_table <- sanity_check(distance_table, df)
print(sanity_table)
write.csv(sanity_table, "duplicate_positions_sanity_table.csv", row.names = FALSE)

#####################################################################################################################################################################################
# QC + transcriptomic metrics for duplicate groups (VBC groups; includes mapped MAPseq hash ID)
####################################################################################################
# Load segmentation QC data
good_cells <- read.csv("LC_visualQC_barcoded_cells.csv", stringsAsFactors = FALSE)
# Standardize uid column as character (safer joins)
good_cells$uid <- as.character(good_cells$uid)
good_uids <- good_cells$uid[good_cells$good_barcoded == 1]
bad_uids  <- good_cells$uid[good_cells$good_barcoded == 0]
cat("Total cells in QC file:", nrow(good_cells), "\n")
cat("Good segmented cells (good_barcoded == 1):", length(good_uids), "\n")
cat("Bad segmented cells (good_barcoded == 0):", length(bad_uids), "\n")

# Helper: base uid from df rownames (e.g. "11_43_1941753.4" -> "11_43_1941753")
get_base_uid <- function(cell_id) sub("\\.[0-9]+$", "", cell_id)

# Segmentation QC lookup for a vector of cells (rownames(df))
check_segmentation <- function(cells, qc_df = good_cells) {
  base_uid <- get_base_uid(cells)
  
  # match base_uid -> qc_df
  m <- match(base_uid, qc_df$uid)
  
  status <- ifelse(
    is.na(m), "Not found in QC list",
    ifelse(qc_df$good_barcoded[m] == 1, "Good", "Bad (badly segmented)")
  )
  
  data.frame(
    CellID = cells,
    Base_UID = base_uid,
    Segmentation_Status = status,
    stringsAsFactors = FALSE
  )
}

# Transcriptomic metrics
#    - total_umis: sum of COUNTS
#    - sparsity: fraction of genes with 0 counts
get_transcriptomic_metrics <- function(base_uid, sce = LCNEneurons, assay_name = "counts") {
  if (!base_uid %in% colnames(sce)) {
    return(list(total_umis = NA_real_, sparsity = NA_real_))
  }
  
  x <- assay(sce, assay_name)[, base_uid]
  
  # If assay returns sparse matrix column, keep numeric ops safe:
  x <- as.numeric(x)
  
  total_umis <- sum(x, na.rm = TRUE)
  sparsity   <- mean(x == 0, na.rm = TRUE)
  
  list(
    total_umis = round(total_umis, 2),
    sparsity   = round(sparsity, 3)
  )
}

check_transcriptomics <- function(cells, sce = LCNEneurons, assay_name = "counts") {
  base_uid <- get_base_uid(cells)
  metrics <- lapply(base_uid, get_transcriptomic_metrics, sce = sce, assay_name = assay_name)
  
  data.frame(
    CellID = cells,
    Base_UID = base_uid,
    Total_UMIs = vapply(metrics, `[[`, numeric(1), "total_umis"),
    Sparsity   = vapply(metrics, `[[`, numeric(1), "sparsity"),
    stringsAsFactors = FALSE
  )
}

#  Build one tidy table for ALL duplicate groups
#    - group ids: both VBC id and mapped MAPseq hash id
#    - ONE ROW PER CELL within each duplicate group 
dup_qc_table <- dplyr::bind_rows(lapply(names(dup_vbc_groups), function(vbc_id) {
  cells <- dup_vbc_groups[[vbc_id]]
  
  seg  <- check_segmentation(cells, qc_df = good_cells)
  tran <- check_transcriptomics(cells, sce = LCNEneurons, assay_name = "counts")
  
  out <- dplyr::left_join(seg, tran, by = c("CellID", "Base_UID"))
  
  out$VBC_group_id    <- vbc_id
  out$MAPseq_group_id <- unname(vbc_to_mapseq[[vbc_id]])
  
  # match your plotting facet label style (useful for debugging/merging)
  out$Group_Label <- paste0("VBC_", vbc_id, " | MAPseqHash_", unname(vbc_to_mapseq[[vbc_id]]))
  
  # optional: add coords from df (since df rownames are CellID)
  coords <- df[cells, c("CCF_ML", "CCF_DV", "CCF_AP", "inRH"), drop = FALSE]
  coords <- tibble::rownames_to_column(as.data.frame(coords), "CellID")
  
  out <- dplyr::left_join(out, coords, by = "CellID")
  
  out
}))

print(head(dup_qc_table, 10))
write.csv(dup_qc_table, "duplicate_groups_detailed.csv", row.names = FALSE)

# Hardcoded manual annotations: Named vectors for CellID -> Doublet_status and Final_selection
# Based on the updated duplicate_groups_detailed.csv. Keys are CellID.
manual_doublet_status <- c(
  "10_40_2680973.4" = "same",
  "10_40_2790931.4" = "same",
  "9_34_1280890.4" = "same",
  "9_34_1390821.4" = "same",
  "5_18_831314.4" = "different",
  "5_19_1590526.2" = "different",
  "5_20_2412260.1" = "same",
  "5_20_2420371.1" = "same",
  "9_35_2160997.3" = "same",
  "9_35_2271333.3" = "same",
  "9_34_1280954.2" = "same",
  "9_34_1390899.2" = "same",
  "11_41_401501.4" = "same",
  "11_41_401512.3" = "same",
  "11_43_1941572.3" = "different",
  "11_43_1941700.4" = "different",
  "11_43_1941708.4" = "different",
  "11_43_1941664.2" = "different",
  "11_43_1941681.4" = "different",
  "11_43_1941644.2" = "different",
  "11_43_1941650.2" = "different",
  "11_44_2711431.2" = "different",
  "9_35_2160966.3" = "same",
  "9_35_2271284.3" = "same",
  "5_20_2420651.3" = "same",
  "5_20_2420670.3" = "same",
  "9_35_2161030.2" = "same",
  "9_35_2271374.2" = "same"
)

manual_final_selection <- c(
  "10_40_2680973.4" = "discard",
  "10_40_2790931.4" = "keep",
  "9_34_1280890.4" = "keep",
  "9_34_1390821.4" = "discard",
  "5_18_831314.4" = "keep",
  "5_19_1590526.2" = "discard",
  "5_20_2412260.1" = "keep",
  "5_20_2420371.1" = "discard",
  "9_35_2160997.3" = "discard",
  "9_35_2271333.3" = "keep",
  "9_34_1280954.2" = "keep",
  "9_34_1390899.2" = "discard",
  "11_41_401501.4" = "discard",
  "11_41_401512.3" = "discard",
  "11_43_1941572.3" = "discard",
  "11_43_1941700.4" = "discard",
  "11_43_1941708.4" = "keep",
  "11_43_1941664.2" = "discard",
  "11_43_1941681.4" = "discard",
  "11_43_1941644.2" = "discard",
  "11_43_1941650.2" = "discard",
  "11_44_2711431.2" = "discard",
  "9_35_2160966.3" = "discard",
  "9_35_2271284.3" = "keep",
  "5_20_2420651.3" = "keep",
  "5_20_2420670.3" = "discard",
  "9_35_2161030.2" = "discard",
  "9_35_2271374.2" = "keep"
)

# Preconditions
stopifnot(exists("dup_qc_table"))
stopifnot("CellID" %in% names(dup_qc_table))
stopifnot(exists("manual_doublet_status"))
stopifnot(exists("manual_final_selection"))

# Normalize keys (force character; trim whitespace)
dup_qc_table <- dup_qc_table %>%
  mutate(CellID = trimws(as.character(CellID)))
mds <- as.character(manual_doublet_status)
names(mds) <- trimws(as.character(names(manual_doublet_status)))
mfs <- as.character(manual_final_selection)
names(mfs) <- trimws(as.character(names(manual_final_selection)))

stopifnot(!is.null(names(mds)), !is.null(names(mfs)))
stopifnot(!any(is.na(names(mds))), !any(is.na(names(mfs))))
stopifnot(!any(names(mds) == ""), !any(names(mfs) == ""))

# Literal overlap sanity check (hard-stop if 0)
n_overlap_doublet <- sum(dup_qc_table$CellID %in% names(mds))
n_overlap_final   <- sum(dup_qc_table$CellID %in% names(mfs))

cat("Overlap (Doublet_status):", n_overlap_doublet, "of", nrow(dup_qc_table), "\n")
cat("Overlap (Final_selection):", n_overlap_final,   "of", nrow(dup_qc_table), "\n")

if (n_overlap_doublet == 0 || n_overlap_final == 0) {
  stop(
    "0 overlap detected. Example dup_qc_table CellIDs: ",
    paste(head(dup_qc_table$CellID, 10), collapse = ", "),
    " | Example manual keys: ",
    paste(head(names(mds), 10), collapse = ", ")
  )
}

# Attach annotations (match-safe; avoids factor indexing pitfalls)
dup_qc_table$Doublet_status  <- unname(mds[dup_qc_table$CellID])
dup_qc_table$Final_selection <- unname(mfs[dup_qc_table$CellID])

dup_qc_table$Doublet_status[is.na(dup_qc_table$Doublet_status)  | dup_qc_table$Doublet_status  == ""] <- "not_reviewed"
dup_qc_table$Final_selection[is.na(dup_qc_table$Final_selection) | dup_qc_table$Final_selection == ""] <- "not_reviewed"

# Validate / summarize
cat("\nDoublet_status counts:\n")
print(table(dup_qc_table$Doublet_status, useNA = "ifany"))
cat("\nFinal_selection counts:\n")
print(table(dup_qc_table$Final_selection, useNA = "ifany"))

missing_doublet <- dup_qc_table %>% filter(Doublet_status == "not_reviewed") %>% pull(CellID)
missing_final   <- dup_qc_table %>% filter(Final_selection == "not_reviewed") %>% pull(CellID)
cat("\nCells missing Doublet_status:", length(missing_doublet), "\n")
if (length(missing_doublet)) print(missing_doublet)
cat("\nCells missing Final_selection:", length(missing_final), "\n")
if (length(missing_final)) print(missing_final)

extra_doublet <- setdiff(names(mds), unique(dup_qc_table$CellID))
extra_final   <- setdiff(names(mfs), unique(dup_qc_table$CellID))
cat("\nManual Doublet_status entries not in dup_qc_table:", length(extra_doublet), "\n")
if (length(extra_doublet)) print(extra_doublet)
cat("\nManual Final_selection entries not in dup_qc_table:", length(extra_final), "\n")
if (length(extra_final)) print(extra_final)

write.csv(dup_qc_table, "dup_qc_table_with_manual_annotations.csv", row.names = FALSE)

# Generate blacklist from manual selections
manual_blacklist <- dup_qc_table$CellID[dup_qc_table$Final_selection == "discard"]
blacklist_df <- data.frame(row_name = manual_blacklist)
write.csv(blacklist_df, "duplicate_blacklist_informed_final.csv", row.names = FALSE)

#####################################################################################################################################################################################
#combine RH and LH matrices
df_subset_long_combined <- rbind(df_subset_longRH, df_subset_longLH)
df_subset_long_combined <- df_subset_long_combined %>%
  mutate(Cell_Source = paste(Cell, Source, sep = "_")) 
#plot all examples on the same figure
p <- ggplot(df_subset_long_combined, aes(x = ROI, y = Value, group = Cell_Source, color = Source)) +
  geom_line() + geom_ribbon(aes(ymin = lower, ymax = upper, fill = Source), colour = NA, alpha = 0.5) +
  facet_wrap(~ Cell_Source, scales = "free_y") +  # Use new column for faceting
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 1),# Set x-axis text size
        axis.text.y = element_text(size = 2),# Set y-axis text size
        axis.title.x = element_text(size = 8),  # Set x-axis label text size
        axis.title.y = element_text(size = 8),  # Set y-axis label text size
        strip.text = element_text(size = 3),  # Reduce facet label font size
        panel.grid.major.x = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +  # Remove minor grid lines
  ylab('projection strength (barcode count)') +
  xlab('ctx projection target (ROIs)') +
  ggtitle('RH and LH somas projecting to ctx') +
  scale_color_manual(values = c("RH" = "#f58231", "LH" = "#4363d8")) +  # Set colors for RH and LH
  scale_fill_manual(values = c("RH" = "#f58231", "LH" = "#4363d8"))  # Set fill colors for RH and LH
print(p)

ggsave("RH_LH_rawcount_example_cortex_projections.pdf", plot = p, device = "pdf",width = 14, height = 10)

range_df <- df_subset_long_combined %>%
  group_by(Cell_Source, Source) %>%
  summarise(range_value = max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE), .groups = 'drop')

# Add log10 of range (add small value to avoid log(0))
range_df <- range_df %>%
  mutate(log_range = log10(range_value + 1e-6))

# Plot histogram of log-transformed ranges
hist_plot <- ggplot(range_df, aes(x = log_range, fill = Source)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  scale_fill_manual(values = c("RH" = "#f58231", "LH" = "#4363d8")) +
  theme_minimal() +
  xlab("log10(Range of projection strength across ctx ROIs)") +
  ylab("Number of cells") +
  ggtitle("Distribution of log10(range) of projection strength by Source")
print(hist_plot)

ggsave("RH_LH_rawcount_hist_ctx.pdf", plot = hist_plot, device = "pdf",width = 8, height = 8)


################################### check for cross-hemispheric projections ###########################################
################################### contralateral projections #########################################################
# --- Cross-hemispheric projections: LH somas projecting to RH ctx ---
# Select LH somas (inRH == 0), look for ctx ROIs in RH
ctx_RH_cols <- RH_meta %>% filter(str_detect(base_roi, "ctx")) %>% pull(colname)
ctx_RH_cols <- ctx_RH_cols[-length(ctx_RH_cols)] #ctx is unique to RH, not present in LH
df_LH <- df[df$inRH == 0, ctx_RH_cols, drop = FALSE]
df_LH[] <- lapply(df_LH, as.numeric)
row_exceeds_LH_cross_RH <- apply(df_LH, 1, function(x) any(x > 50))
df_LH_cross <- df_LH[row_exceeds_LH_cross_RH, , drop = FALSE]
df_LH_cross <- df_LH_cross[rowSums(df_LH_cross) > 0, , drop = FALSE]
dim(df_LH_cross)

# Pivot for plotting
df_LH_cross_long <- df_LH_cross %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value))
df_LH_cross_long$ROI <- factor(df_LH_cross_long$ROI, levels = ctx_RH_cols)

p_LH_cross <- ggplot(df_LH_cross_long, aes(x = ROI, y = Value, group = Cell)) +
  geom_line(color = "#f58231") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "#f58231") +
  facet_wrap(~ Cell, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2),
        axis.text.y = element_text(size = 4),# Set y-axis text size
        axis.title.x = element_text(size = 8),  # Set x-axis label text size
        axis.title.y = element_text(size = 8),  # Set y-axis label text size
        strip.text = element_text(size = 5),  # Reduce facet label font size
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab('projection strength (barcode count)') +
  xlab('ctx projection target (ROIs)') +
  ggtitle('LH somas projecting to RH ctx')
print(p_LH_cross)

ggsave("LHsoma_rawcount_RHcortex_projections.pdf", plot = p_LH_cross, device = "pdf", width = 12, height = 10)

# --- Cross-hemispheric projections: RH somas projecting to LH ctx ---
ctx_LH_cols <- LH_meta %>% filter(str_detect(base_roi, "ctx")) %>% pull(colname)
df_RH <- df[df$inRH == 1, ctx_LH_cols, drop = FALSE]
df_RH[] <- lapply(df_RH, as.numeric)
row_exceeds_RH_cross_LH <- apply(df_RH, 1, function(x) any(x > 50))
df_RH_cross <- df_RH[row_exceeds_RH_cross_LH, , drop = FALSE]
df_RH_cross <- df_RH_cross[rowSums(df_RH_cross) > 0, , drop = FALSE]
dim(df_RH_cross)

df_RH_cross_long <- df_RH_cross %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value))
df_RH_cross_long$ROI <- factor(df_RH_cross_long$ROI, levels = ctx_LH_cols)

p_RH_cross <- ggplot(df_RH_cross_long, aes(x = ROI, y = Value, group = Cell)) +
  geom_line(color = "#4363d8") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "#4363d8") +
  facet_wrap(~ Cell, scales = "free_y") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2),
        axis.text.y = element_text(size = 4),# Set y-axis text size
        axis.title.x = element_text(size = 8),  # Set x-axis label text size
        axis.title.y = element_text(size = 8),  # Set y-axis label text size
        strip.text = element_text(size = 5),  # Reduce facet label font size
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  ylab('projection strength (barcode count)') +
  xlab('ctx projection target (ROIs)') +
  ggtitle('RH somas projecting to LH ctx')
print(p_RH_cross)

ggsave("RHsoma_rawcount_LHcortex_projections.pdf", plot = p_RH_cross, device = "pdf", width = 12, height = 10)


############################################### check ipsi- and contralateral side by side ###################################
# --- LH somata: contralateral (RH ctx) and ipsilateral (LH ctx) projections ---
row_names_LH_crossing_RH <- rownames(df_LH_cross)
ctx_RH_cols <- RH_meta %>% filter(str_detect(base_roi, "ctx")) %>% pull(colname)
ctx_RH_cols <- ctx_RH_cols[-length(ctx_RH_cols)] # drop unique RH ctx column if needed
ctx_LH_cols <- LH_meta %>% filter(str_detect(base_roi, "ctx")) %>% pull(colname)

df_subset_contra_LH <- df[row_names_LH_crossing_RH, ctx_RH_cols, drop = FALSE]
df_subset_ipsi_LH   <- df[row_names_LH_crossing_RH, ctx_LH_cols, drop = FALSE]
df_subset_contra_LH[] <- lapply(df_subset_contra_LH, as.numeric)
df_subset_ipsi_LH[]   <- lapply(df_subset_ipsi_LH, as.numeric)

RH_ctx_meta <- RH_meta %>% filter(colname %in% ctx_RH_cols) %>% arrange(idx)
LH_ctx_meta <- LH_meta %>% filter(colname %in% ctx_LH_cols) %>% arrange(idx)
df_subset_contra_LH <- df_subset_contra_LH[, RH_ctx_meta$colname]
df_subset_ipsi_LH   <- df_subset_ipsi_LH[, LH_ctx_meta$colname]

df_subset_long_contra_LH <- df_subset_contra_LH %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         ProjSide = "RH") # Contralateral

df_subset_long_ipsi_LH <- df_subset_ipsi_LH %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         ProjSide = "LH") # Ipsilateral

df_combined_LH <- bind_rows(df_subset_long_contra_LH, df_subset_long_ipsi_LH)
df_combined_LH$ROI <- factor(df_combined_LH$ROI, levels = unique(c(RH_ctx_meta$colname, LH_ctx_meta$colname)))

# --- RH somata: contralateral (LH ctx) and ipsilateral (RH ctx) projections ---
row_names_RH_crossing_LH <- rownames(df_RH_cross)
df_subset_contra_RH <- df[row_names_RH_crossing_LH, ctx_LH_cols, drop = FALSE]
df_subset_ipsi_RH   <- df[row_names_RH_crossing_LH, ctx_RH_cols, drop = FALSE]
df_subset_contra_RH[] <- lapply(df_subset_contra_RH, as.numeric)
df_subset_ipsi_RH[]   <- lapply(df_subset_ipsi_RH, as.numeric)

df_subset_contra_RH <- df_subset_contra_RH[, LH_ctx_meta$colname]
df_subset_ipsi_RH   <- df_subset_ipsi_RH[, RH_ctx_meta$colname]

df_subset_long_contra_RH <- df_subset_contra_RH %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         ProjSide = "LH") # Contralateral

df_subset_long_ipsi_RH <- df_subset_ipsi_RH %>%
  rownames_to_column(var = "Cell") %>%
  pivot_longer(cols = -Cell, names_to = "ROI", values_to = "Value") %>%
  mutate(lower = Value - sqrt(Value),
         upper = Value + sqrt(Value),
         ProjSide = "RH") # Ipsilateral

df_combined_RH <- bind_rows(df_subset_long_contra_RH, df_subset_long_ipsi_RH)
df_combined_RH$ROI <- factor(df_combined_RH$ROI, levels = unique(c(LH_ctx_meta$colname, RH_ctx_meta$colname)))

# --- Plotting for LH somata ---
colors <- c("RH" = "#f58231", "LH" = "#4363d8")
p_LH <- ggplot(df_combined_LH, aes(x = ROI, y = Value, group = Cell)) +
  geom_line(aes(color = ProjSide)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ProjSide), alpha = 0.5) +
  facet_grid(ProjSide ~ Cell, scales = "fixed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 5, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab('projection strength (barcode count)') +
  xlab('ctx projection target (ROIs)') +
  ggtitle('LH soma: ipsi- and contralateral cortex projections')
print(p_LH)

ggsave("LHsoma_rawcount_RH-LHprojections_side_by_side.pdf", plot = p_LH, device = "pdf", width = 16, height = 9)

# --- Plotting for RH somata ---
p_RH <- ggplot(df_combined_RH, aes(x = ROI, y = Value, group = Cell)) +
  geom_line(aes(color = ProjSide)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ProjSide), alpha = 0.5) +
  facet_grid(ProjSide ~ Cell, scales = "fixed") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 2),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        strip.text = element_text(size = 4, angle = 45),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  ylab('projection strength (barcode count)') +
  xlab('ctx projection target (ROIs)') +
  ggtitle('RH soma: ipsi- and contralateral cortex projections')
print(p_RH)

ggsave("RHsoma_rawcount_LH-RHprojections_side_by_side.pdf", plot = p_RH, device = "pdf", width = 16, height = 9)


########################################################### plot histograms with highest projection density regions for raw and log norm ##########################################################
# Drop unwanted columns
drop_cols <- c("slice", "CCF_AP", "CCF_DV", "CCF_ML", "inRH")
df_clean <- df[ , !(names(df) %in% drop_cols)]
# Convert all columns to numeric
df_numeric <- as.data.frame(lapply(df_clean, as.numeric))
# Log-normalize the full dataframe
df_log <- log10(1 + df_numeric * 100)
# New normalization: per-neuron (row) sum, then log1p, then max normalize
row_sums <- rowSums(df_numeric)
df_norm <- sweep(df_numeric, 1, row_sums, "/") # divide each row by its sum
df_norm[is.na(df_norm)] <- 0 # handle division by zero
df_norm_log <- log1p(df_norm) # safer log transform
max_val <- max(df_norm_log, na.rm = TRUE)
df_norm_log_max <- df_norm_log / max_val
# Compute column sums for all three metrics
raw_sums <- colSums(df_numeric)
log_sums <- colSums(df_log)
norm_sums <- colSums(df_norm_log_max)
# Get top 50 and bottom 10 regions for each metric
top50_raw <- names(sort(raw_sums, decreasing = TRUE))[1:50]
bottom10_raw <- rev(names(sort(raw_sums, decreasing = FALSE))[1:10])
regions_raw <- c(top50_raw, bottom10_raw)
top50_log <- names(sort(log_sums, decreasing = TRUE))[1:50]
bottom10_log <- rev(names(sort(log_sums, decreasing = FALSE))[1:10])
regions_log <- c(top50_log, bottom10_log)
top50_norm <- names(sort(norm_sums, decreasing = TRUE))[1:50]
bottom10_norm <- rev(names(sort(norm_sums, decreasing = FALSE))[1:10])
regions_norm <- c(top50_norm, bottom10_norm)

# Calculate means for all regions
all_means_raw <- colMeans(df_numeric)
all_means_log <- colMeans(df_log)
all_means_norm <- colMeans(df_norm_log_max)
# Reorder means to match selected regions
means_raw <- all_means_raw[regions_raw]
means_log <- all_means_log[regions_log]
means_norm <- all_means_norm[regions_norm]

# Plot barplots with horizontal bars, light blue fill, and red dashed separator
par(mfrow = c(1, 3), mar = c(6, 12, 4, 2)) # 3 plots side by side
offset <- 0.65 # increase this value to move labels further left/up
# Raw data barplot
bp1 <- barplot(means_raw, horiz = TRUE, las = 1, main = "Raw Data", 
               xlab = "Mean Value", names.arg = regions_raw, cex.names = 1.2, col = "lightblue")
abline(h = 50.5, col = "red", lty = 2, lwd = 2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 25, labels = "top50", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 55, labels = "bottom10", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
mtext("Raw BC counts", side = 1, line = 4, cex = 1)
# Log-normalized data barplot
bp2 <- barplot(means_log, horiz = TRUE, las = 1, main = "Log10(1 + 100x)", 
               xlab = "Mean Value", names.arg = regions_log, cex.names = 1.2, col = "lightblue")
abline(h = 50.5, col = "red", lty = 2, lwd = 2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 25, labels = "top50", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 55, labels = "bottom10", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
mtext("Log-normalized BC counts", side = 1, line = 4, cex = 1)
# Per-neuron normalized, log, max-normalized barplot
bp3 <- barplot(means_norm, horiz = TRUE, las = 1, main = "RowNorm+log1p+MaxNorm", 
               xlab = "Mean Value", names.arg = regions_norm, cex.names = 1.2, col = "lightblue")
abline(h = 50.5, col = "red", lty = 2, lwd = 2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 25, labels = "top50", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
text(x = par("usr")[1] - offset * diff(par("usr")[1:2]), y = 55, labels = "bottom10", srt = 90, adj = 0.5, xpd = TRUE, col = "blue", cex = 1.2)
mtext("Row-normalized, log1p, max-normalized", side = 1, line = 4, cex = 1)
dev.copy(pdf, "top_bottom_innervated_areas_by_norm_method.pdf", width = 12, height = 8)
dev.off()

