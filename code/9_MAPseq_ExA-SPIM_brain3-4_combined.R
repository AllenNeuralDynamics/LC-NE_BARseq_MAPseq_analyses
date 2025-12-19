# Load functions which handle pre-processing or organizing of the data
source("~/capsule/code/01_loaders_brain3-4_combined.R")

brain3 <- load_data_brain3()
names(brain3)
# Extract components
brain3_raw_proj_matrix <- as.data.frame(brain3$proj_matrix_raw)
# Preserve original rownames
rownames(brain3_raw_proj_matrix) <- rownames(brain3$proj_matrix_raw)
brain3_inRH_lookup <- brain3$inRH_lookup
brain3_metadata <- brain3$metadata

brain4 <- load_data_brain4()
names(brain4)
# Extract components
brain4_raw_proj_matrix <- as.data.frame(brain4$proj_matrix_raw)
# Preserve original rownames
rownames(brain4_raw_proj_matrix) <- rownames(brain4$proj_matrix_raw)
brain4_inRH_lookup <- brain4$inRH_lookup
brain4_metadata <- brain4$metadata

#set working directory to save things to
setwd('/scratch/BARseq_780345-780346_combined/')

# Ensure ROI compatibility between brain3 and brain4 ROIs 
colnames(brain3_raw_proj_matrix)
colnames(brain4_raw_proj_matrix)
# Correct Amyg-GPe spelling
colnames(brain4_raw_proj_matrix) <- stringr::str_replace_all(colnames(brain4_raw_proj_matrix), "amyg-Gpe", "amyg-GPe")
# Sum cerebellum sections into one sample
groups <- c("RH.12", "LH.12")
for (group in groups) {
  i_col <- paste0("cerebellum I_", group)
  ii_col <- paste0("cerebellum II_", group)
  new_col <- paste0("cerebellum_", group)
  if (i_col %in% colnames(brain4_raw_proj_matrix) && ii_col %in% colnames(brain4_raw_proj_matrix)) {
    # Sum the projections
    brain4_raw_proj_matrix[, new_col] <- brain4_raw_proj_matrix[, i_col] + brain4_raw_proj_matrix[, ii_col]
    # Remove the original I and II columns
    brain4_raw_proj_matrix <- brain4_raw_proj_matrix[, !(colnames(brain4_raw_proj_matrix) %in% c(i_col, ii_col))]
    cat("Summed", i_col, "+", ii_col, "into", new_col, "\n")
  } else {
    cat("Warning: Columns for", group, "not found\n")
  }
}

# Sum columns with corresponding base ROI names to enhance dataset compatibility - brain3
dim(brain3_raw_proj_matrix)
head(brain3_raw_proj_matrix)
brain3_result <- sum_by_base_roi(brain3_raw_proj_matrix)
brain3_summed <- brain3_result$summed_matrix
dim(brain3_summed)
head(brain3_summed)
brain3_mapping <- brain3_result$mapping
cat("Brain3 mapping (what was combined into each base ROI):\n")
for (base in names(brain3_mapping)) {
  cat(base, ":", paste(brain3_mapping[[base]], collapse = ", "), "\n")
}

# Sum columns with corresponding base ROI names to enhance dataset compatibility - brain4
dim(brain4_raw_proj_matrix)
head(brain4_raw_proj_matrix)
brain4_result <- sum_by_base_roi(brain4_raw_proj_matrix)
brain4_summed <- brain4_result$summed_matrix
dim(brain4_summed)
head(brain4_summed)
brain4_mapping <- brain4_result$mapping
cat("\nBrain4 mapping:\n")
for (base in names(brain4_mapping)) {
  cat(base, ":", paste(brain4_mapping[[base]], collapse = ", "), "\n")
}

# Apply ipsi-contra to brain3 summed matrix
ipsi_contra_brain3 <- create_ipsi_contra_from_matrix(brain3_summed, brain3_inRH_lookup, "brain3 summed")
head(ipsi_contra_brain3)

# Apply ipsi-contra to brain4 summed matrix
ipsi_contra_brain4 <- create_ipsi_contra_from_matrix(brain4_summed, brain4_inRH_lookup, "brain4 summed")
head(ipsi_contra_brain4)

# Get column names from the ipsi-contra matrices
brain3_cols <- colnames(ipsi_contra_brain3)
brain4_cols <- colnames(ipsi_contra_brain4)
# Find shared regions (exact matches)
shared_regions <- intersect(brain3_cols, brain4_cols)
# Find unique regions
unique_brain3 <- setdiff(brain3_cols, brain4_cols)
unique_brain4 <- setdiff(brain4_cols, brain3_cols)
# Print results
cat("Shared regions (", length(shared_regions), "):\n")
if (length(shared_regions) > 0) {
  cat(paste(sort(shared_regions), collapse = ", "), "\n\n")
}
cat("Unique to brain3 (", length(unique_brain3), "):\n")
if (length(unique_brain3) > 0) {
  cat(paste(sort(unique_brain3), collapse = ", "), "\n\n")
}
cat("Unique to brain4 (", length(unique_brain4), "):\n")
if (length(unique_brain4) > 0) {
  cat(paste(sort(unique_brain4), collapse = ", "), "\n\n")
}

# Subset to shared regions
brain3_shared <- ipsi_contra_brain3[, shared_regions, drop = FALSE]
head(brain3_shared)
brain4_shared <- ipsi_contra_brain4[, shared_regions, drop = FALSE]
head(brain4_shared)

# Normalize each subset
brain3_shared_norm <- normalize_projection_matrix(brain3_shared, "brain3 shared ipsi-contra")
head(brain3_shared_norm)
brain4_shared_norm <- normalize_projection_matrix(brain4_shared, "brain4 shared ipsi-contra")
head(brain4_shared_norm)

# Combine the normalized subsets
combined_norm <- as.data.frame(rbind(brain3_shared_norm, brain4_shared_norm))
head(combined_norm)
cat("Combined normalized matrix dimensions:", dim(combined_norm), "\n")

# Preserve rownames before conversion
original_rownames <- rownames(combined_norm)
# Convert to numeric matrix
combined_norm <- as.matrix(combined_norm)
combined_norm <- apply(combined_norm, 2, as.numeric)  # Ensure numeric
combined_norm <- as.matrix(combined_norm)
# Reassign rownames
rownames(combined_norm) <- original_rownames
head(combined_norm) 

# Combine metadata from brain3 and brain4
combined_metadata <- rbind(brain3_metadata, brain4_metadata)

# Check if combined_metadata is in the same order as combined_norm based on rownames
order_match <- all(rownames(combined_norm) == combined_metadata$row_id)
cat("Are the rownames of combined_norm and row_id in combined_metadata in the same order?", order_match, "\n")
if (!order_match) {
  combined_metadata <- combined_metadata[match(rownames(combined_norm), combined_metadata$row_id), ]
  cat("Reordered combined_metadata to match combined_norm rownames.\n")
}
head(combined_metadata)
# Save combined_metadata to a CSV file in the working directory
write.csv(combined_metadata, file = "MAPseq_combined_metadata.csv", row.names = FALSE)

##########################################################################################################################################
# Working with combined ipsilateral and contralateral matrices
full_df <- combined_norm

# Function to extract base region name from column names (for matching only)
get_base_region <- function(col_name) {
  # Remove -contra or -ipsi, optionally followed by .digits
  name <- sub("-(contra|ipsi)(\\.\\d+)?$", "", col_name)
  # Then remove trailing .digits if any (for SP form, but adjust as needed)
  name <- sub("\\.\\d+$", "", name)
  name
}

region_map <- list(
  OLF = c("olf.bulb", "AON"),
  Isocortex = c("motor.ctx", "orb.ctx", "ctx.1", "ctx.2", "ctx.3", "ctx", "cc"),
  HPF = c("hippocampus"),
  CTXsp = c("amygdala", "amyg.GPe"),
  CNU = c("CPu", "NAc", "septum", "BNST"),
  TH = c("thalamus"),
  HY = c("hypothalamus"),
  MB = c("midbrain"),
  CB = c("cerebellum"),
  P = c("hindbrain"),
  MY = c("medulla"),
  SP = c("sp.cord.1_SP", "sp.cord.2_SP", "sp.cord.3_SP")
)

# Initialize an empty data frame for the combined regions
combined_df <- data.frame(matrix(ncol = length(region_map), nrow = nrow(full_df)))
colnames(combined_df) <- names(region_map)
rownames(combined_df) <- rownames(full_df)

# Sum the values for each region group based on region_map
for (region in names(region_map)) {
  # Find columns where the base region name (after removing suffixes) matches
  matching_cols <- colnames(full_df)[
    get_base_region(colnames(full_df)) %in% region_map[[region]]
  ]
  
  if (length(matching_cols) > 0) {
    combined_df[[region]] <- rowSums(full_df[, matching_cols, drop = FALSE], na.rm = TRUE)
  } else {
    combined_df[[region]] <- 0  # If no matching columns, assign zero
  }
}

# View the combined data frame
head(combined_df)

#Keep track of which columns contributed to each region
for (region in names(region_map)) {
  matching_cols <- colnames(full_df)[
    get_base_region(colnames(full_df)) %in% region_map[[region]]
  ]
  cat("Region", region, "includes:", paste(matching_cols, collapse = ", "), "\n")
}

# Normalize each row by its sum (row-wise normalization)
normalized_df <- combined_df / rowSums(combined_df)

# Convert to numeric matrix (if not already)
normalized_mat <- as.matrix(normalized_df)
rownames(normalized_mat) <- rownames(normalized_df)

# Find the region with the maximum projection for each neuron
top_proj_indices <- apply(normalized_mat, 1, function(x) {
  idx <- which.max(x)
  if(length(idx) > 1) idx <- idx[1] # take the first if there are ties
  return(idx)
})
top_proj_indices <- as.integer(top_proj_indices)
top_proj <- colnames(normalized_mat)[top_proj_indices]

# Define region order (adjust to match your preferred order)
region_order <- c("OLF","Isocortex","HPF","CTXsp","CNU","TH","HY","MB","CB","P","MY","SP")
top_proj_factor <- factor(top_proj, levels = region_order)

# Sort neurons by their top projection region
sorted_indices <- order(top_proj_factor)
normalized_mat_sorted <- normalized_mat[sorted_indices, ]

# Plot heatmap
p <- pheatmap(
  t(normalized_mat_sorted),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(256, option = "magma"),
  fontsize_row = 14,
  fontsize_col = 2,
  main = "Normalized Projection Strengths (Sorted by Top Projection)",
  labels_col = rownames(normalized_mat_sorted)
)

# PDF output (vector format, great for publications)
pdf("sorted_proj_heatmap_ipsi-contra.pdf", width = 10, height = 6)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

# PNG output (raster, with DPI)
png("sorted_proj_heatmap_ipsi-contra.png", width = 3000, height = 1800, res = 300)
grid::grid.newpage()
grid::grid.draw(p$gtable)
dev.off()

corr_matrix <- cor(normalized_mat_sorted, method = "spearman", use = "pairwise.complete.obs")
head(corr_matrix)

icefire <- c(
  "#bce6da", "#b9e4d9", "#b7e2d9", "#b4e0d8", "#b1ded7", "#afdcd7", "#acdad6", "#a9d8d5", "#a6d6d5", "#a3d4d4", 
  "#a0d2d3", "#9ed1d3", "#9bcfd2", "#98cdd2", "#95cbd1", "#92c9d1", "#8fc7d0", "#8cc5d0", "#89c4cf", "#86c2cf", 
  "#83c0cf", "#80bece", "#7dbcce", "#7abbce", "#77b9cd", "#74b7cd", "#71b5cd", "#6eb4cd", "#6bb2cd", "#68b0cd", 
  "#65aecd", "#62accd", "#5faacd", "#5ca9cd", "#5aa7cd", "#57a5cd", "#55a3cd", "#52a1cd", "#509fcd", "#4d9dcd", 
  "#4b9bcd", "#4999cd", "#4797cd", "#4595cd", "#4393ce", "#4191ce", "#3f8fce", "#3d8dce", "#3b8bce", "#3a89cf", 
  "#3987cf", "#3885cf", "#3783cf", "#3680d0", "#367ed0", "#367ccf", "#3779cf", "#3777cf", "#3875ce", "#3972ce", 
  "#3b70cd", "#3c6ecc", "#3d6bca", "#3f69c9", "#4067c7", "#4264c5", "#4362c3", "#4460c0", "#455ebe", "#465bbb", 
  "#4759b9", "#4857b6", "#4855b3", "#4953b0", "#4952ac", "#4950a9", "#494ea5", "#494da1", "#494b9d", "#484a99", 
  "#474895", "#474791", "#46468d", "#45448a", "#444386", "#434282", "#42407e", "#413f7a", "#403e77", "#3f3d73", 
  "#3e3c70", "#3d3a6c", "#3b3969", "#3a3865", "#393762", "#38365f", "#36345c", "#353359", "#343255", "#333152", 
  "#323050", "#302f4d", "#2f2e4a", "#2e2d47", "#2d2c44", "#2c2b42", "#2b2a3f", "#29293d", "#28283b", "#272738", 
  "#262636", "#262534", "#252432", "#242330", "#23222e", "#22222c", "#22212a", "#212129", "#212027", "#201f26", 
  "#201f24", "#1f1f23", "#1f1e22", "#1f1e21", "#1f1e20", "#1e1e1f", "#1e1e1e", "#1e1e1e", "#1f1e1e", "#201d1e", 
  "#211d1e", "#221d1e", "#231d1e", "#241d1e", "#251d1e", "#271d1f", "#281d1f", "#2a1e20", "#2b1e20", "#2d1e21", 
  "#2f1e21", "#301f22", "#321f23", "#342024", "#362024", "#382025", "#3a2126", "#3c2127", "#3e2228", "#412229", 
  "#43232a", "#45232b", "#48242c", "#4a242d", "#4c252e", "#4f252f", "#512630", "#542731", "#572732", "#592833", 
  "#5c2835", "#5e2836", "#612937", "#642938", "#672a39", "#692a3a", "#6c2b3b", "#6f2b3c", "#722b3d", "#752c3d", 
  "#782c3e", "#7b2c3f", "#7e2c40", "#812d41", "#842d41", "#872d42", "#8a2d42", "#8d2d43", "#902d43", "#932e44", 
  "#962e44", "#992e44", "#9c2e44", "#9f2e44", "#a22e44", "#a52f44", "#a82f43", "#ab3043", "#ad3142", "#b03142", 
  "#b33241", "#b63340", "#b8343f", "#bb353f", "#bd373e", "#c0383d", "#c2393c", "#c43b3b", "#c73d3a", "#c93e39", 
  "#cb4038", "#cd4238", "#cf4437", "#d14636", "#d34936", "#d54b35", "#d74d34", "#d85034", "#da5234", "#dc5433", 
  "#dd5733", "#df5933", "#e05c32", "#e25e32", "#e36133", "#e46433", "#e56633", "#e76934", "#e86c35", "#e96f36", 
  "#ea7237", "#eb7539", "#eb783b", "#ec7a3d", "#ed7d3f", "#ee8042", "#ee8345", "#ef8648", "#ef894b", "#f08c4e", 
  "#f18f51", "#f19154", "#f29457", "#f2975b", "#f39a5e", "#f49d62", "#f49f66", "#f5a269", "#f5a56d", "#f6a871", 
  "#f6aa74", "#f7ad78", "#f7b07c", "#f8b27f", "#f9b583", "#f9b887", "#fabb8b", "#fabd8e", "#fbc092", "#fbc396", 
  "#fcc699", "#fcc89d", "#fdcba1", "#fdcea4", "#fed1a8", "#fed3ac")

# Number of colors in your palette
n_colors <- length(icefire)
# Create breaks from -1 to 1, with zero exactly in the middle
breaks <- seq(-1, 1, length.out = n_colors + 1)

p <- pheatmap(
  corr_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = icefire,
  breaks = breaks,
  main = "Co-innervation (Spearman Correlation)",
  fontsize_row = 14,
  fontsize_col = 14,
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("-1", "0", "1")
)
# Add rotated text next to the color bar
grid.text("Rank correlation", x = 0.992, y = 0.83, rot = 90, gp = gpar(fontsize = 11))

# PDF output (vector format, great for publications)
pdf("rank_corr_ipsi-contra.pdf", width = 8, height = 8)
grid::grid.newpage()
grid::grid.draw(p$gtable)
grid.text("Rank correlation", x = 0.992, y = 0.83, rot = 90, gp = gpar(fontsize = 11))
dev.off()

# PNG output (raster, with DPI)
png("rank_corr_ipsi-contra.png", width = 2400, height = 2400, res = 300)
grid::grid.newpage()
grid::grid.draw(p$gtable)
grid.text("Rank correlation", x = 0.992, y = 0.83, rot = 90, gp = gpar(fontsize = 11))
dev.off()

# Create a data frame with cell information and their top projections
cell_top_projections <- data.frame(
  cell_id = rownames(normalized_df),
  top_projection = top_proj,
  top_projection_strength = apply(normalized_mat, 1, max),
  stringsAsFactors = FALSE
)
# View the results
head(cell_top_projections)
dim(cell_top_projections)
write.csv(cell_top_projections, "cell_top_projections.csv", row.names = FALSE)
# Summary of cells by top projection region
table(cell_top_projections$top_projection)

# Use row_id from metadata which matches cell_id exactly
cell_top_projections_with_coords <- merge(
  cell_top_projections, 
  combined_metadata[, c("row_id", "CCF_DV", "CCF_ML", "CCF_AP")], 
  by.x = "cell_id", 
  by.y = "row_id", 
  all.x = TRUE
)
# View the result
head(cell_top_projections_with_coords)
dim(cell_top_projections_with_coords)
write.csv(cell_top_projections_with_coords, "cell_top_projections_with_coords.csv", row.names = FALSE)

