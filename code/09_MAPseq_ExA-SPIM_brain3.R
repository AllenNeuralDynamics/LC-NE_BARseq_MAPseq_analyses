library(tidyverse)
library(viridis)
library(pheatmap)
library(grid)
library(RColorBrewer)
library(ggplot2)
library(gplots)

setwd('/Users/polina.kosillo/Seurat_projects/BarSeq/780345_full_brain_dataset')
source("loaders_brain3_v3.R")

# setwd('/Users/polina.kosillo/Seurat_projects/BarSeq/780346_full_brain_dataset')
# source("loaders_brain4.R")

dat <- load_data()
names(dat)
raw_matrix <- dat$proj_matrix_raw           # Raw counts 
#log_matrix <- dat$mat_log_ordered          # Log-transformed 
#rownorm_matrix <- dat$mat_rownorm_ordered  # Row-normalized 
inRH_info <- dat$inRH_lookup               # Hemisphere lookup
metadata <- dat$metadata                   # Cell metadata

cat("Raw projection matrix dimensions:", dim(raw_matrix), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

# Convert projection matrices to ipsi-contra format
ipsi_contra_raw <- create_ipsi_contra_from_raw(dat$proj_matrix_raw, dat$inRH_lookup)
#ipsi_contra_log <- create_ipsi_contra_from_log(dat$mat_log_ordered, dat$inRH_lookup)
#ipsi_contra_rownorm <- create_ipsi_contra_from_rownorm(dat$mat_rownorm_ordered, dat$inRH_lookup)

###################### Only use ipsilateral data ##########################
# # Only keep 'ipsi' or 'sp.cord' columns, drop contralateral data
# ipsi_df <- ipsi_contra_raw[, grepl("ipsi|sp.cord", colnames(ipsi_contra_raw))]
# # Remove '-ipsi' from column names
# colnames(ipsi_df) <- gsub("-ipsi", "", colnames(ipsi_df))
# colnames(ipsi_df)
# 
# region_map <- list(
#   OLF = c("olf.bulb", "olf.bulb.1", "AON"),
#   Isocortex = c("motor.ctx", "orb.ctx", "ctx.1", "ctx.2", "ctx.3", "ctx.1.1", "ctx.2.1", "ctx.3.1", "ctx.1.2", "ctx.2.2", "ctx.3.2", "ctx.2.3", "ctx.3.3", "cc", "cc.1", "cc.2"),
#   HPF = c("hippocampus", "hippocampus.1"),
#   CTXsp = c("ctx.1.3", "amyg.GPe", "ctx.1.4", "amygdala"),
#   CNU = c("CPu", "septum", "NAc", "CPu.1", "septum.1", "NAc.1", "CPu.2", "septum.2", "BNST"),
#   TH = c("thalamus","thalamus.1"),
#   HY = c("hypothalamus", "hypothalamus.1"),
#   MB = c("midbrain", "midbrain.1", "midbrain.2"),
#   CB = c("cerebellum"),
#   P = c( "hindbrain", "hindbrain.1", "hindbrain.2"),
#   MY = c("medulla"),
#   SP = c( "sp.cord.1_SP", "sp.cord.2_SP", "sp.cord.3_SP")
# )

# region_map <- list(
#   OLF = c("olf.bulb", "olf.bulb.1", "AON", "AON.1"),
#   Isocortex = c("motor.ctx", "orb.ctx", "ctx.1", "ctx.2", "ctx.3", "ctx.1.1", "ctx.2.1", "ctx.3.1", "ctx.1.2", "ctx.2.2", "ctx.3.2", "ctx.2.3", "ctx.3.3", "cc", "cc.1", "ctx.3.4", "ctx.2.4", "ctx.1.4", "ctx.3.5", "ctx" ),
#   HPF = c("hippocampus", "hippocampus.1", "hippocampus.2"),
#   CTXsp = c("ctx.1.3", "amyg.GPe", "ctx.1.4", "amygdala"),
#   CNU = c("CPu", "septum", "NAc", "CPu.1", "septum.1", "NAc.1", "CPu.2", "septum.2", "BNST"),
#   TH = c("thalamus","thalamus.1"),
#   HY = c("hypothalamus", "hypothalamus.1"), 
#   MB = c("midbrain", "midbrain.1", "midbrain.2"), 
#   CB = c("cerebellum.I", "cerebellum.II"),
#   P = c( "hindbrain", "hindbrain.1", "hindbrain.2"),
#   MY = c("medulla", "medulla.1"),
#   SP = c( "sp.cord.1_SP", "sp.cord.2_SP", "sp.cord.3_SP")
# )

# # Initialize an empty data frame for the combined regions
# combined_df <- data.frame(matrix(ncol = length(region_map), nrow = nrow(ipsi_df)))
# colnames(combined_df) <- names(region_map)
# rownames(combined_df) <- rownames(ipsi_df)
# # Sum the values for each region group based on region_map
# for (region in names(region_map)) {
#   matching_cols <- intersect(region_map[[region]], colnames(ipsi_df))  # Ensure columns exist in the data
#   if (length(matching_cols) > 0) {
#     combined_df[[region]] <- rowSums(ipsi_df[, matching_cols, drop = FALSE], na.rm = TRUE)
#   } else {
#     combined_df[[region]] <- 0  # If no matching columns, assign zero
#   }
# }
# # View the combined data frame
# head(combined_df)

# Working with combined ipsilateral and contralateral matrices
# Work with the full dataframe (keeping original column names)
full_df <- ipsi_contra_raw

# Function to extract base region name from column names (for matching only)
get_base_region <- function(col_name) {
  # Remove -contra/-ipsi and any trailing numbers for matching purposes
  gsub("(-contra|-ipsi)(\\.\\d+)?$", "", col_name)
}

region_map <- list(
  OLF = c("olf.bulb", "olf.bulb.1", "AON"),
  Isocortex = c("motor.ctx", "orb.ctx", "ctx.1", "ctx.2", "ctx.3", "ctx.1.1", "ctx.2.1", "ctx.3.1", "ctx.1.2", "ctx.2.2", "ctx.3.2", "ctx.2.3", "ctx.3.3", "cc", "cc.1", "cc.2"),
  HPF = c("hippocampus", "hippocampus.1"),
  CTXsp = c("ctx.1.3", "amyg.GPe", "ctx.1.4", "amygdala"),
  CNU = c("CPu", "septum", "NAc", "CPu.1", "septum.1", "NAc.1", "CPu.2", "septum.2", "BNST"),
  TH = c("thalamus","thalamus.1"),
  HY = c("hypothalamus", "hypothalamus.1"),
  MB = c("midbrain", "midbrain.1", "midbrain.2"),
  CB = c("cerebellum"),
  P = c( "hindbrain", "hindbrain.1", "hindbrain.2"),
  MY = c("medulla"),
  SP = c( "sp.cord.1_SP", "sp.cord.2_SP", "sp.cord.3_SP")
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
write.csv(cell_top_projections, "cell_top_projections.csv", row.names = FALSE)
# Summary of cells by top projection region
table(cell_top_projections$top_projection)

# Use row_id from metadata which matches cell_id exactly
cell_top_projections_with_coords <- merge(
  cell_top_projections, 
  metadata[, c("row_id", "CCF_DV", "CCF_ML", "CCF_AP")], 
  by.x = "cell_id", 
  by.y = "row_id", 
  all.x = TRUE
)
# View the result
head(cell_top_projections_with_coords)
write.csv(cell_top_projections_with_coords, "cell_top_projections_with_coords.csv", row.names = FALSE)
