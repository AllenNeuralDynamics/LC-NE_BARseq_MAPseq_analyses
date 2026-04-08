# Load functions which handle pre-processing or organizing of the data
source("~/capsule/code/01_loaders_brain3.R")

#set working directory 
setwd('/results/BARseq_780345/')

dat <- load_data()
names(dat)
raw_matrix <- dat$proj_matrix_raw           # Raw counts 
log_matrix <- dat$mat_log_ordered          # Log-transformed 
rownorm_matrix <- dat$mat_rownorm_ordered  # Row-normalized 
inRH_info <- dat$inRH_lookup               # Hemisphere lookup
metadata <- dat$metadata                   # Cell metadata

cat("Raw projection matrix dimensions:", dim(raw_matrix), "\n")
cat("Metadata dimensions:", dim(metadata), "\n")

# Create a custom color palette for the inRH column
inRH_colors <- ifelse(inRH_info$inRH == 1, "red", "blue")
numeric_matrix <- data.matrix(log_matrix)
stopifnot(all(rownames(log_matrix) == inRH_info$cell_id))
# Sanity check RH and LH projections map
# Plot the heatmap with the custom color palette for the inRH column
heatmap.2(numeric_matrix, 
          distfun = function(x) proxy::dist(x, method = "cosine"), 
          hclustfun = function(x) hclust(x, method = "ward.D"), 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "none", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.6, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,  #rev(brewer.pal(11, "RdYlBu")),
          RowSideColors = inRH_colors,
          main = "Log transformed projections map with inRH column")


# Convert projection matrices to ipsi-contra format
ipsi_contra_raw <- create_ipsi_contra_from_raw(dat$proj_matrix_raw, dat$inRH_lookup)
ipsi_contra_log <- create_ipsi_contra_from_log(dat$mat_log_ordered, dat$inRH_lookup)
ipsi_contra_rownorm <- create_ipsi_contra_from_rownorm(dat$mat_rownorm_ordered, dat$inRH_lookup)

#Sanity check ipsi-contra projections map
#plot the heatmap to visualize ipsi and contralateral patterns
numeric_matrix <- data.matrix(ipsi_contra_log)
heatmap.2(numeric_matrix, 
          distfun = function(x) proxy::dist(x, method = "cosine"), 
          hclustfun = function(x) hclust(x, method = "ward.D"), 
          scale = "none", 
          Colv = NA, 
          Rowv = TRUE, 
          dendrogram = "none", 
          trace = "none", 
          density.info = "histogram", 
          cexRow = 0.3, 
          cexCol = 0.6, 
          keysize = 1,
          margins = c(10, 8), 
          col = viridis,  #rev(brewer.pal(11, "RdYlBu")),
          main = "Log transformed projections map")


#Sort by top projection
# Assign a numeric order to each region (column)
region_order <- setNames(seq_along(colnames(ipsi_contra_log)), colnames(ipsi_contra_log))
# Find the top projection region for each cell
top_proj_region <- colnames(ipsi_contra_log)[apply(ipsi_contra_log, 1, which.max)]
# Map each cell's top region to its numeric order
top_proj_order <- region_order[top_proj_region]
# Add as new columns
ipsi_contra_log$top_proj_order <- top_proj_order
ipsi_contra_log$top_proj_region <- top_proj_region
# Sort the data frame by top_proj_order
ipsi_contra_log_sorted <- ipsi_contra_log[order(ipsi_contra_log$top_proj_order), ]
projection_cols <- setdiff(colnames(ipsi_contra_log_sorted), c("top_proj_order", "top_proj_region"))
heatmap_matrix <- as.matrix(ipsi_contra_log_sorted[, projection_cols])
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) <- rownames(ipsi_contra_log_sorted)
n_cells <- nrow(heatmap_matrix)
main_title <- paste0(
  "Cell-by-Region Projections\n(Cells Sorted by Top Projection)\n",
  "n = ", n_cells, " cells")

heatmap.2(
  heatmap_matrix,
  col = viridis(256),
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = main_title,
  margins = c(10, 10),
  cexRow = 0.2,
  cexCol = 0.5,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  key.title = "Log10(1 + 100*Count)"
)
dev.copy(pdf, "Ipsi-contra projections heatmap top region sorted.pdf", width = 12, height = 10)
dev.off()

# Select ipsi columns and spinal cord columns
ipsi_cols <- grep("-ipsi$", colnames(ipsi_contra_log), value = TRUE)
spinal_cols <- grep("_SP$", colnames(ipsi_contra_log), value = TRUE)
cols_to_plot <- c(ipsi_cols, spinal_cols[spinal_cols %in% colnames(ipsi_contra_log)])
ipsi_matrix <- ipsi_contra_log[, cols_to_plot, drop = FALSE]
# Assign a numeric order to each ipsi region (column)
region_order_ipsi <- setNames(seq_along(colnames(ipsi_matrix)), colnames(ipsi_matrix))
# Find the top projection region for each cell (among ipsi only)
top_proj_region_ipsi <- colnames(ipsi_matrix)[apply(ipsi_matrix, 1, which.max)]
# Map each cell's top region to its numeric order
top_proj_order_ipsi <- region_order_ipsi[top_proj_region_ipsi]
# Add as new columns
ipsi_matrix$top_proj_order <- top_proj_order_ipsi
ipsi_matrix$top_proj_region <- top_proj_region_ipsi
# Sort the data frame by top_proj_order
ipsi_matrix_sorted <- ipsi_matrix[order(ipsi_matrix$top_proj_order), ]
projection_cols_ipsi <- setdiff(colnames(ipsi_matrix_sorted), c("top_proj_order", "top_proj_region"))
heatmap_matrix_ipsi <- as.matrix(ipsi_matrix_sorted[, projection_cols_ipsi])
heatmap_matrix_ipsi <- apply(heatmap_matrix_ipsi, 2, as.numeric)
rownames(heatmap_matrix_ipsi) <- rownames(ipsi_matrix_sorted)
n_cells <- nrow(heatmap_matrix_ipsi)
main_title <- paste0(
  "Cell-by-Region Projections\n(Cells Sorted by Top Ipsi Projection)\n",
  "n = ", n_cells, " cells"
)
heatmap.2(
  heatmap_matrix_ipsi,
  col = viridis(256),
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = main_title,
  margins = c(10, 10),
  cexRow = 0.2,
  cexCol = 0.7,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  key.title = "Log10(1 + 100*Count)"
)
dev.copy(pdf, "Ipsi-only projections heatmap top region sorted.pdf", width = 10, height = 8)
dev.off()

################################## Regional grouping and probability analysis #########################################
# Check if any rows sum to zero (should be none due to filtering in loader)
zero_rows <- rowSums(ipsi_contra_raw) == 0
cat("Number of zero-sum rows found:", sum(zero_rows), "\n")

# Drop rows which sum up to zero
if(sum(zero_rows) > 0) {
  cat("Removing", sum(zero_rows), "zero-sum rows...\n")
  ipsi_contra_clean <- ipsi_contra_raw[!zero_rows, ]
} else {
  cat("No zero-sum rows found - using data as-is\n")
  ipsi_contra_clean <- ipsi_contra_raw
}
# after cleaning
source_matrix <- ipsi_contra_clean

region_map <- list(
  olf_bulb_contra = c("olf.bulb-contra", "olf.bulb-contra.1"),
  AON_contra = c("AON-contra"),
  frontal_ctx_contra = c("orb.ctx-contra", "motor.ctx-contra"),
  dorsal_ctx_contra = c("ctx.3-contra", "ctx.3-contra.1", "ctx.3-contra.2", "ctx.3-contra.3"),
  middle_ctx_contra = c("ctx.2-contra", "ctx.2-contra.1", "ctx.2-contra.2", "ctx.2-contra.3"),
  ventral_ctx_contra = c("ctx.1-contra", "ctx.1-contra.1", "ctx.1-contra.2", "ctx.1-contra.3", "ctx.1-contra.4"),
  corpus_callosum_contra = c("cc-contra", "cc-contra.1", "cc-contra.2"),
  striatum_contra = c("CPu-contra", "CPu-contra.1", "CPu-contra.2"),
  NAc_contra = c("NAc-contra", "NAc-contra.1"),
  septum_BNST_contra = c("septum-contra", "septum-contra.1", "septum-contra.2", "BNST-contra"),
  hippocampus_contra = c("hippocampus-contra", "hippocampus-contra.1"),
  amygdala_contra = c("amygdala-contra", "amyg.GPe-contra"),
  thalamus_contra = c("thalamus-contra", "thalamus-contra.1"),
  hypothalamus_contra = c("hypothalamus-contra", "hypothalamus-contra.1"),
  midbrain_contra = c("midbrain-contra", "midbrain-contra.1", "midbrain-contra.2"),
  hindbrain_contra = c("hindbrain-contra", "hindbrain-contra.1", "hindbrain-contra.2"),
  cerebellum_contra = c("cerebellum-contra"),
  medulla_contra = c("medulla-contra"),
  
  olf_bulb_ipsi = c("olf.bulb-ipsi", "olf.bulb-ipsi.1"),
  AON_ipsi = c("AON-ipsi"),
  frontal_ctx_ipsi = c("orb.ctx-ipsi", "motor.ctx-ipsi"),
  dorsal_ctx_ipsi = c("ctx.3-ipsi", "ctx.3-ipsi.1", "ctx.3-ipsi.2", "ctx.3-ipsi.3"),
  middle_ctx_ipsi = c("ctx.2-ipsi", "ctx.2-ipsi.1", "ctx.2-ipsi.2", "ctx.2-ipsi.3"),
  ventral_ctx_ipsi = c("ctx.1-ipsi", "ctx.1-ipsi.1", "ctx.1-ipsi.2", "ctx.1-ipsi.3", "ctx.1-ipsi.4"),
  corpus_callosum_ipsi = c("cc-ipsi", "cc-ipsi.1", "cc-ipsi.2"),
  striatum_ipsi = c("CPu-ipsi", "CPu-ipsi.1", "CPu-ipsi.2"),
  NAc_ipsi = c("NAc-ipsi", "NAc-ipsi.1"),
  septum_BNST_ipsi = c("septum-ipsi", "septum-ipsi.1", "septum-ipsi.2", "BNST-ipsi"),
  hippocampus_ipsi = c("hippocampus-ipsi", "hippocampus-ipsi.1"),
  amygdala_ipsi = c("amygdala-ipsi", "amyg.GPe-ipsi"),
  thalamus_ipsi = c("thalamus-ipsi", "thalamus-ipsi.1"),
  hypothalamus_ipsi = c("hypothalamus-ipsi", "hypothalamus-ipsi.1"),
  midbrain_ipsi = c("midbrain-ipsi", "midbrain-ipsi.1", "midbrain-ipsi.2"),
  hindbrain_ipsi = c("hindbrain-ipsi", "hindbrain-ipsi.1", "hindbrain-ipsi.2"),
  cerebellum_ipsi = c("cerebellum-ipsi"),
  medulla_ipsi = c("medulla-ipsi"),
  
  spinal_cord_cervical = c("sp.cord.1_SP"),
  spinal_cord_thoracic = c( "sp.cord.2_SP", "sp.cord.3_SP")
)

################################## Regional grouping and probability analysis #########################################
# Starting from ipsi_contra_raw: cell × fine-grained region MAPseq counts (projection strength).
# Here we pool columns into anatomically-defined region groups (region_map)
# and SUM their raw counts per cell.
# Result: combined_df is cell × grouped-region with RAW counts (projection strength),
#         not yet normalized, but aggregated across subcolumns for each group.
combined_df <- data.frame(matrix(ncol = length(region_map), nrow = nrow(source_matrix)))
colnames(combined_df) <- names(region_map)
rownames(combined_df) <- rownames(source_matrix)
# Sum the values for each region group based on region_map
for (region in names(region_map)) {
  matching_cols <- intersect(region_map[[region]], colnames(source_matrix))  # Ensure columns exist in the data
  if (length(matching_cols) > 0) {
    combined_df[[region]] <- rowSums(source_matrix[, matching_cols, drop = FALSE], na.rm = TRUE)
  } else {
    combined_df[[region]] <- 0  # If no matching columns, assign zero
  }
}
# View the combined data frame which shows summed raw projection strength per grouped region (hemisphere-specific) per cell.
head(combined_df)

###################################################### Basic statistics: cell and region counts ######################################################
# Binarize combined_df: treat any non-zero MAPseq count as "yes, this cell projects to this grouped region".
# Number of grouped regions each cell projects to (presence/absence, not strength)
num_regions_per_cell <- rowSums(combined_df > 0)
# Number of cells projecting to each grouped region (presence/absence, not strength)
num_cells_per_region <- colSums(combined_df > 0)

# Plot bar plot for number of CELLS per grouped region
# This shows which regions receive projections from the largest number of cells
par(mar = c(14, 4, 4, 2)) 
barplot(num_cells_per_region, las = 2, col = "lightblue", 
        xlab = "", ylab = "Number of cells", main = "Number of Cells per Region")
mtext("Region", side = 1, line = 12) 
par(mar = c(5, 4, 4, 2))  
dev.copy(pdf, "Absolute number of cells with non-zero projections to ROI.pdf", width = 10, height = 8)
dev.off()

###################################################### Co-innervation matrix (Jaccard index) → overlap among cells that hit at least one of the pair ######################################################
# Ignore projection strength and work with binary presence/absence:
# a grouped region is "innervated" by a cell if the SUMMED raw count for that group is > 0.
# For each pair of grouped regions (i, j), we compute the Jaccard index across cells:
#   J(i, j) = (# cells that project to BOTH i and j) / (# cells that project to i OR j) - Fraction of cells projecting to both regions, relative to those projecting to either
# This measures how often two grouped regions are co-innervated by the same cells.
co_innervation_matrix <- matrix(0, nrow = ncol(combined_df), ncol = ncol(combined_df), 
                                dimnames = list(colnames(combined_df), colnames(combined_df)))
for (i in seq_along(colnames(combined_df))) {
  for (j in seq_along(colnames(combined_df))) {
    proj_i <- combined_df[, i] > 0
    proj_j <- combined_df[, j] > 0
    co_innervation_matrix[i, j] <- sum(proj_i & proj_j) / sum(proj_i | proj_j)
  }
}
#  Fraction of cells that project to both the corresponding pair of regions, identify regions that are commonly co-innervated by the same cells
co_innervation_matrix[is.nan(co_innervation_matrix)] <- 0
co_innervation_matrix_masked <- co_innervation_matrix
diag(co_innervation_matrix_masked) <- NA
heatmap.2(
  co_innervation_matrix_masked,
  col = viridis(256),
  na.color = "grey",  # NA values will appear grey
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "Normalized Co-Innervation Fraction (Jaccard)",
  margins = c(10, 10),
  cexRow = 0.7,
  cexCol = 0.7,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  keysize = 1.0,
  key.title = "Jaccard Index",
  key.xlab = "Value",
  density.info = "none",
  lmat = rbind(c(4, 3), c(2, 1)),
  lhei = c(1, 4),
  lwid = c(2, 4)
)
dev.copy(pdf, "Co-innervation fraction.pdf", width = 12, height = 10)
dev.off()

###################################################### Projection fraction per region  ######################################################
# Restrict to multi-projecting cells: cells that project to >1 grouped region.
# Binarize combined_df (combined_df > 0), focus on presence/absence of projections, not projection strength.
# For each grouped region compute: fraction of multi-projecting cells that project to that region at all, otherwise asking which regions are most targeted by multi-projecting cells?
# Number of grouped regions each cell projects to (binary: >0 counts)
num_regions_per_cell <- rowSums(combined_df > 0)
# Indices of multi-projecting cells (cells with projections to >1 grouped region)
multi_proj_cells <- which(num_regions_per_cell > 1)
# For each grouped region, count how many multi-projecting cells project there (binary presence)
num_multi_proj_per_region <- colSums((combined_df > 0)[multi_proj_cells, , drop = FALSE])
# Fraction of multi-projecting cells that target each grouped region
projection_fraction <- num_multi_proj_per_region / length(multi_proj_cells)
# Create a data frame for ggplot
projection_fraction_df <- data.frame(
  Region = names(projection_fraction),
  ProjectionFraction = projection_fraction
)
# Ensure Region retains the order from combined_df
projection_fraction_df$Region <- factor(projection_fraction_df$Region, levels = names(projection_fraction))
# This plot shows, for each region, the fraction of multi-projecting cells that target it
ggplot(projection_fraction_df, aes(x = Region, y = ProjectionFraction)) +
  geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7, color = "black") +
  labs(title = "Fraction of Multi-Projecting Cells Targeting Each Region",
       x = "Region",
       y = "Fraction of Multi-Projecting Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)))
dev.copy(pdf, "Projection fraction of multi-projecting cells.pdf", width = 10, height = 8)
dev.off()

###################################################### Conditional probability matrix ######################################################
# Goal: For each pair of grouped regions (A,B), compute
# P(cell projects to B | cell projects to A)
# using CELL-BASED, BINARY data (presence/absence of any projection), not projection strength.
# Fraction of cells that project to other regions given that they project to a specific region - If a cell projects to A, how likely to also project to B?
# Convert combined_df to binary 
df_bin <- (combined_df > 0) * 1
# Initialize probability matrix
# rows = "given" region A, columns = "also in" region B
prob_other_given_here <- matrix(0, nrow = ncol(df_bin), ncol = ncol(df_bin), 
                                dimnames = list(colnames(df_bin), colnames(df_bin)))
# Number of cells projecting to each region (binary, not strength)
cells_in_this_region <- colSums(df_bin)

# Compute conditional probabilities P(B | A)
for (region_A in colnames(df_bin)) {
  # Subset to cells that project to region_A
  df_temp <- df_bin[df_bin[, region_A] == 1, , drop = FALSE]
  
  # If no cells project to region_A, skip
  if (nrow(df_temp) == 0) next
  
  # Among region_A-positive cells, count how many also project to each OTHER region
  df_other <- colSums(df_temp[, colnames(df_temp) != region_A, drop = FALSE])
  
  # Normalize by the total number of region_A-positive cells:
  #   P(B | A) = (# cells with A AND B) / (# cells with A)
  prob_other_given_here[region_A, names(df_other)] <- df_other / cells_in_this_region[region_A]
}

# Clean up numeric mode (df_bin was already numeric, but this is harmless)
prob_other_given_here <- as.matrix(prob_other_given_here)
mode(prob_other_given_here) <- "numeric"
prob_other_given_here[is.nan(prob_other_given_here)] <- 0
prob_other_given_here[is.infinite(prob_other_given_here)] <- 0

# Inspect a few rows + denominators
print(prob_other_given_here[1:3, ])
print(cells_in_this_region)

# Plot:
# Each row A shows, for all B, the conditional probability that a cell projecting to A also projects to B.
par(mar = c(5, 8, 4, 2))
prob_other_given_here_masked <- prob_other_given_here
diag(prob_other_given_here_masked) <- NA  # we don't plot P(A | A) here

heatmap.2(
  prob_other_given_here_masked,
  col = viridis(256),
  na.color = "grey",
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "Conditional Probability Heatmap P(B | A)",
  xlab = "Cell is also in these regions (B)",
  ylab = "Cell is given to be in this region (A)",
  margins = c(10, 10),
  cexRow = 0.7,
  cexCol = 0.8,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  keysize = 1.0,
  key.title = "Probability",
  key.xlab = "Probability",
  density.info = "none",
  lmat = rbind(c(4, 3), c(2, 1)),
  lhei = c(1, 4),
  lwid = c(2, 4)
)
dev.copy(pdf, "Conditional probability between regions.pdf", width = 12, height = 10)
dev.off()

###################################################### Co-occurrence probability matrix (intersection/all) → absolute prevalence of a co-projection pattern in the whole population ######################################################
# Goal: For each pair of grouped regions (i, j), compute
#   P(cell projects to BOTH i AND j)
# i.e. (# cells that project to both i and j) / (# all cells).
# This is a CELL-BASED, BINARY measure (presence/absence),
# and it treats each hemisphere-specific grouped region (e.g. striatum_contra, striatum_ipsi)
# as a separate node.

# Binary presence/absence for ALL grouped regions (including spinal cord)
brain_regions <- colnames(combined_df)  
brain_binary <- (combined_df[, brain_regions] > 0) * 1

# Initialize co-occurrence probability matrix
co_occur_prob <- matrix(0, nrow = length(brain_regions), ncol = length(brain_regions),
                        dimnames = list(brain_regions, brain_regions))
# Compute P(projects to i AND j) for each pair, normalized by total # cells
for (i in seq_along(brain_regions)) {
  for (j in seq_along(brain_regions)) {
    # Fraction of all cells that project to BOTH regions i and j
    both_regions <- sum(brain_binary[, i] & brain_binary[, j])
    co_occur_prob[i, j] <- both_regions / nrow(brain_binary)
  }
}

# Plot the co-occurrence matrix with biological ordering (not clustering)
# Define the biological order: contra regions, then ipsi regions, then spinal cord
biological_order <- c(
  "olf_bulb_contra", "AON_contra", "frontal_ctx_contra", "dorsal_ctx_contra", 
  "middle_ctx_contra", "ventral_ctx_contra", "corpus_callosum_contra", 
  "striatum_contra", "NAc_contra", "septum_BNST_contra", "hippocampus_contra", 
  "amygdala_contra", "thalamus_contra", "hypothalamus_contra", "midbrain_contra", 
  "hindbrain_contra", "cerebellum_contra", "medulla_contra",
  "olf_bulb_ipsi", "AON_ipsi", "frontal_ctx_ipsi", "dorsal_ctx_ipsi", 
  "middle_ctx_ipsi", "ventral_ctx_ipsi", "corpus_callosum_ipsi", 
  "striatum_ipsi", "NAc_ipsi", "septum_BNST_ipsi", "hippocampus_ipsi", 
  "amygdala_ipsi", "thalamus_ipsi", "hypothalamus_ipsi", "midbrain_ipsi", 
  "hindbrain_ipsi", "cerebellum_ipsi", "medulla_ipsi",
  "spinal_cord_cervical", "spinal_cord_thoracic"
)
# (Optional safety) only keep regions that actually exist in combined_df
brain_regions_ordered <- intersect(biological_order, brain_regions)
# Reorder the matrix according to biological organization
co_occur_prob_ordered <- co_occur_prob[biological_order, biological_order]

co_occur_prob_ordered_masked <- co_occur_prob_ordered
diag(co_occur_prob_ordered_masked) <- NA # hide P(i AND i) = P(i)
heatmap.2(
  co_occur_prob_ordered_masked,
  col = viridis(256),
  na.color = "grey",
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "Co-Occurrence Fraction (Intersection / All Cells)",
  margins = c(12, 12),
  cexRow = 0.6,
  cexCol = 0.6,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  key.title = "Fraction of All Cells",
  key.xlab = "Fraction",
  density.info = "none"
)
dev.copy(pdf, "Co-occurence matrix of relationships between regions.pdf", width = 12, height = 10)
dev.off()

############################# Compare Jaccard co-innervation vs co-occurrence (intersection / all cells) ###############################################################
# Ensure both matrices have the same order
common_regions <- intersect(colnames(co_innervation_matrix), colnames(co_occur_prob))
co_inn_mat <- co_innervation_matrix[common_regions, common_regions]
co_occur_mat <- co_occur_prob[common_regions, common_regions]

# Difference: Jaccard minus co-occurrence
#   diff_mat[i,j] = J(i,j) - C(i,j)
#   where J(i,j) = (# cells A∧B)/(# cells A∨B),
#         C(i,j) = (# cells A∧B)/(# all cells)
diff_mat <- co_inn_mat - co_occur_mat
# Ratio (log2 scale): Jaccard / co-occurrence
#   log2_ratio_mat[i,j] = log2( J(i,j) / C(i,j) )
# Add a small pseudocount to avoid division by zero when C(i,j) = 0
pseudocount <- 1e-6
log2_ratio_mat <- log2((co_inn_mat + pseudocount) / (co_occur_mat + pseudocount))

# Difference heatmap
# Shows where normalized overlap among A∨B cells (Jaccard) is larger or smaller
# than the absolute prevalence among all cells (co-occurrence).
heatmap.2(
  diff_mat,
  col = viridis(256),
  na.color = "grey",
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "Jaccard Co-Innervation minus Co-Occurrence",
  margins = c(10, 10),
  cexRow = 0.7,
  cexCol = 0.7,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  key.title = "Difference",
  key.xlab = "Value",
  density.info = "none"
)
dev.copy(pdf, "Jaccard_co-innervation_minus_co-occurrence.pdf", width = 12, height = 10)
dev.off()

# Ratio heatmap (log2 scale)
# Shows fold-change of Jaccard (intersection/union) relative to co-occurrence (intersection/all cells):
#   >0  → Jaccard > co-occurrence
#   <0  → Jaccard < co-occurrence
heatmap.2(
  log2_ratio_mat,
  col = viridis(256),
  na.color = "grey",
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "log2(Jaccard Co-Innervation / Co-Occurrence)",
  margins = c(10, 10),
  cexRow = 0.7,
  cexCol = 0.7,
  Rowv = FALSE,
  Colv = FALSE,
  key = TRUE,
  key.title = "log2(J / C)",
  key.xlab = "log2 Ratio",
  density.info = "none"
)
dev.copy(pdf, "log2_Jaccard_Co-innervation_over_Co-occurrence.pdf", width = 12, height = 10)
dev.off()
