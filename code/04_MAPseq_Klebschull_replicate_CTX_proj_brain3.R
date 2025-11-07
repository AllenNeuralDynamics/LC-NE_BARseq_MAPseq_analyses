################# MAPseq Barcode checks for CTX slices replicating Klebschull paper ##################

#set working directory 
setwd('/scratch/BARseq_780345/')

#####################################################################################################################################################################################
#read in files which contain relevant information for projections such as projection matrix, cell IDs and cluster types
projection_matrix <- readr::read_csv("./MapSeq_matched_projections_1_mismatch.csv", col_names = TRUE)

head(projection_matrix)
LCNEneurons <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
colData(LCNEneurons)
cell_metadata <- as.data.frame(colData(LCNEneurons))

proj_index <- read_tsv("sampleinfo.tsv")
head(proj_index)
sampleinfo <- read_excel("/data/780345_2025-02-20_00-00-00/MAPseq/M295_20250729_USEthis/M295_20250721.sampleinfo.xlsx", sheet = "Sample information")
head(sampleinfo)
colnames(sampleinfo) <- c("usertube", "ourtube", "samplename", "siteinfo", "QC_qPCR", "rtprimer", "brain")
proj_index$MapSeqV1_tube <- sampleinfo$rtprimer[match(proj_index$usertube, sampleinfo$usertube)]

# Format the data by putting it together
# Get the BC columns and match them to brain region names and samplename
bc_cols <- grep("^BC\\d+$", colnames(projection_matrix), value = TRUE)
rtprimer_numbers <- as.integer(sub("BC", "", bc_cols))
# For each BC column, find the corresponding region and samplename
region_names <- proj_index$region[match(rtprimer_numbers, proj_index$MapSeqV1_tube)] 
samplenames  <- proj_index$samplename[match(rtprimer_numbers, proj_index$MapSeqV1_tube)] 
# Combine region and samplename for new column names
region_samplename <- paste(region_names, samplenames, sep = "_")
# Make unique in case of duplicates
region_samplename_unique <- make.unique(region_samplename)
# Rename BC columns to region_samplename
colnames(projection_matrix)[match(bc_cols, colnames(projection_matrix))] <- region_samplename_unique
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
df <- df[, !(colnames(df) %in% c("CellID", "dist", "vbc_read", "louvain_cluster", "row_id"))]
print(colnames(df))

###################################################### barcode counts sanity check ###############################################################
# Plot histogram of raw value distributions for all columns to look at general barcode expression levels
# Specify metadata columns to exclude
exclude_cols <- c("slice", "barcode","CCF-AP", "CCF_DV", "CCF_ML", "inRH")
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
  base_roi <- str_replace(name, "(_RH|_LH).*", "")
  return(data.frame(colname = name, base_roi = base_roi, hemisphere = hemi, stringsAsFactors = FALSE))
}
# Build metadata for all columns with _RH or _LH
meta <- do.call(rbind, lapply(col_names, extract_info)) %>%
  filter(!is.na(hemisphere))
# Add tube index (original column index)
meta$tube_idx <- match(meta$colname, col_names)
# Group by hemisphere and base_roi, then arrange by tube_idx (original order)
meta <- meta %>%
  group_by(hemisphere, base_roi) %>%
  arrange(tube_idx, .by_group = TRUE) %>%
  ungroup()
# Now, arrange meta by the minimum tube_idx for each base_roi group (to anchor to first appearance)
meta <- meta %>%
  group_by(hemisphere, base_roi) %>%
  mutate(min_idx = min(tube_idx)) %>%
  ungroup() %>%
  arrange(hemisphere, min_idx, tube_idx) %>%
  select(-min_idx)
# Assign new indices within each hemisphere
meta <- meta %>%
  group_by(hemisphere) %>%
  mutate(idx = row_number()) %>%
  ungroup()
# For RH columns, in original df order
RH_meta <- meta %>% filter(hemisphere == "RH") %>% arrange(tube_idx)
RH_tubes <- RH_meta$tube_idx
RH_rois  <- RH_meta$base_roi
RH_idx   <- RH_meta$idx  # This is the plotting/grouping index
# For LH columns, in original df order
LH_meta <- meta %>% filter(hemisphere == "LH") %>% arrange(tube_idx)
LH_tubes <- LH_meta$tube_idx
LH_rois  <- LH_meta$base_roi
LH_idx   <- LH_meta$idx

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
RH_ctx_meta <- RH_meta %>% filter(str_detect(base_roi, "ctx"))
RH_ctx_meta <- RH_ctx_meta %>% arrange(idx)
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

LH_ctx_meta <- LH_meta %>% filter(str_detect(base_roi, "ctx"))
LH_ctx_meta <- LH_ctx_meta %>% arrange(idx)
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
plot_cells_in_batches(df_subset_longRH, batch_size = 20, save = TRUE, prefix = "RH_ctx_cells", title_prefix = "RH somas projecting to RH ctx")
plot_cells_in_batches(df_subset_longLH, batch_size = 20, save = TRUE, prefix = "LH_ctx_cells", title_prefix = "LH somas projecting to LH ctx")

########################################################### sum across ctx sections to mimic Justus's slice data and Mathew's ExA-SPIM flow ##########################################################
head(df)
dim(df)
# Define the ROI groups
roi_groups <- list(
  slice1 = c("motor ctx_RH", "orb ctx_RH", "motor ctx_LH", "orb ctx_LH" ),  
  slice2 = c("motor ctx_RH.1", "orb ctx_RH.1", "motor ctx_LH.1", "orb ctx_LH.1"),
  slice3 = c("ctx 1_RH", "ctx 2_RH", "ctx 3_RH", "ctx 1_LH", "ctx 2_LH", "ctx 3_LH"),
  slice4 = c("ctx 1_RH.1", "ctx 2_RH.1", "ctx 3_RH.1", "ctx 1_LH.1", "ctx 2_LH.1", "ctx 3_LH.1"),
  slice5 = c("ctx 1_RH.2", "ctx 2_RH.2", "ctx 3_RH.2", "ctx 1_LH.2", "ctx 2_LH.2", "ctx 3_LH.2"),
  slice6 = c("ctx 1_RH.3", "ctx 2_RH.3", "ctx 3_RH.3", "ctx 1_LH.3", "ctx 2_LH.3", "ctx 3_LH.3"),
  slice7 = c("ctx 1_RH.4", "ctx 2_RH.4", "ctx 3_RH.4", "ctx 1_LH.4", "ctx 2_LH.4", "ctx 3_LH.4" ),
  slice8 = c("ctx_RH")
)
# Create a new data frame with summed values for each slice
ctx_df <- data.frame(row.names = rownames(df))
for (slice in names(roi_groups)) {
  cols <- roi_groups[[slice]]
  # Only keep columns that exist in df (in case some are missing)
  cols_present <- cols[cols %in% colnames(df)]
  ctx_df[[slice]] <- rowSums(df[, cols_present, drop=FALSE], na.rm=TRUE)
}
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
ctx_long$ROI <- factor(ctx_long$ROI, levels = paste0("slice", 1:8))

# Plot function (similar to before)
plot_cells_in_batches(ctx_long, batch_size = 25, save = TRUE, prefix = "ctx_slices", title_prefix = "Cells projecting to ctx slices")

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

compare_full_profiles(dup_ctx, df)
compare_full_profiles(dup_RH, df)
compare_full_profiles(dup_LH, df)

# Estimate how many duplicates can be expected by chance
simulate_duplicate_probability_sparse <- function(N = 517, K = 20, 
                                                  p_zero = 0.95, 
                                                  value_range = 1:10, 
                                                  n_sim = 100) {
  dup_counts <- numeric(n_sim)
  for (i in seq_len(n_sim)) {
    mat <- matrix(
      ifelse(runif(N * K) < p_zero, 0, sample(value_range, N * K, replace = TRUE)),
      nrow = N, ncol = K
    )
    profiles <- apply(mat, 1, paste, collapse = "_")
    dup_counts[i] <- sum(duplicated(profiles))
  }
  mean(dup_counts)
}

# Example usage for ctx subset:
set.seed(42)
simulate_duplicate_probability_sparse(N = 517, K = 20, p_zero = 0.95, value_range = 1:10, n_sim = 100)

simulate_duplicate_probability_sparse(N = 517, K = 120, p_zero = 0.95, value_range = 1:10, n_sim = 100)

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

ggsave("RH_LH_rawcount_example_cortex_projections.png", plot = p, device = "png", dpi = 1200, width = 14, height = 10, units = "in")
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

ggsave("RH_LH_rawcount_hist_ctx.png", plot = hist_plot, device = "png", dpi = 300, width = 8, height = 8, units = "in")
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
ggsave("LHsoma_rawcount_RHcortex_projections.png", plot = p_LH_cross, device = "png", dpi = 1200, width = 12, height = 10, units = "in")
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
ggsave("RHsoma_rawcount_LHcortex_projections.png", plot = p_RH_cross, device = "png", dpi = 1200, width = 12, height = 10, units = "in")
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
  facet_grid(ProjSide ~ Cell, scales = "free_y") +
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
ggsave("LHsoma_rawcount_RH-LHprojections_side_by_side.png", plot = p_LH, device = "png", dpi = 1200, width = 16, height = 9, units = "in")
ggsave("LHsoma_rawcount_RH-LHprojections_side_by_side.pdf", plot = p_LH, device = "pdf", width = 16, height = 9)

# --- Plotting for RH somata ---
p_RH <- ggplot(df_combined_RH, aes(x = ROI, y = Value, group = Cell)) +
  geom_line(aes(color = ProjSide)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = ProjSide), alpha = 0.5) +
  facet_grid(ProjSide ~ Cell, scales = "free_y") +
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
ggsave("RHsoma_rawcount_LH-RHprojections_side_by_side.png", plot = p_RH, device = "png", dpi = 1200, width = 16, height = 9, units = "in")
ggsave("RHsoma_rawcount_LH-RHprojections_side_by_side.pdf", plot = p_RH, device = "pdf", width = 16, height = 9)


########################################################### plot histograms with highest projection density regions for raw and log norm ##########################################################
# Drop unwanted columns
drop_cols <- c("slice", "barcode", "CCF_DV", "CCF_ML", "inRH")
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
offset <- 0.45 # increase this value to move labels further left/up
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
quartz.save("top_bottom_innervated_areas_by_norm_method.png", type = "png", dpi = 600, bg = "white")
