##########################################################################################################################################
source("~/capsule/code/02_prepare_brain3_4_combined_inputs.R")

OUT_DIR <- "/results/BARseq_780345-780346_combined/"

prep <- prepare_brain3_4_inputs(
  loaders_path = "~/capsule/code/01_loaders_brain3-4_combined.R",
  out_dir = OUT_DIR,
  verbose = TRUE,
  return_intermediates = TRUE,
  restore_wd = FALSE
)

setwd(OUT_DIR)  # explicit, guarantees downstream saves go here

combined_norm <- prep$combined_norm
combined_metadata <- prep$combined_metadata     
combined_inRH_lookup <- prep$combined_inRH_lookup

combined_norm_df <- as.data.frame(combined_norm)

# Create a custom color palette for the inRH column (red for RH, blue for LH)
inRH_colors <- ifelse(combined_inRH_lookup$inRH == 1, "red", "blue")

# Plot the combined heatmap (similar to the old sanity check)
numeric_matrix <- data.matrix(combined_norm_df)
heatmap.2(
  numeric_matrix, 
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
  col = viridis(256),  # Assuming viridis is loaded
  RowSideColors = inRH_colors,
  main = "Combined Normalized Ipsi-Contra Projections (Brain3 + Brain4)"
)
# Add legend for RowSideColors (RH/LH soma)
legend("topright", 
       legend = c("RH Soma", "LH Soma"), 
       fill = c("red", "blue"), 
       title = "Soma Hemisphere", 
       cex = 0.8,  # Adjust size if needed
       bty = "n")  # No box around legend
dev.copy(pdf, "Combined_heatmap_with_soma_loc_legend.pdf", width = 12, height = 10)
dev.off()

# Sort by top projection (using combined_norm_df as ipsi_contra_log equivalent)
# Assign a numeric order to each region (column)
region_order <- setNames(seq_along(colnames(combined_norm_df)), colnames(combined_norm_df))
# Find the top projection region for each cell
top_proj_region <- colnames(combined_norm_df)[apply(combined_norm_df, 1, which.max)]
# Map each cell's top region to its numeric order
top_proj_order <- region_order[top_proj_region]
# Add as new columns
combined_norm_df$top_proj_order <- top_proj_order
combined_norm_df$top_proj_region <- top_proj_region
# Sort the data frame by top_proj_order
combined_norm_sorted <- combined_norm_df[order(combined_norm_df$top_proj_order), ]
projection_cols <- setdiff(colnames(combined_norm_sorted), c("top_proj_order", "top_proj_region"))
heatmap_matrix <- as.matrix(combined_norm_sorted[, projection_cols])
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) <- rownames(combined_norm_sorted)
n_cells <- nrow(heatmap_matrix)
main_title <- paste0(
  "Combined Cell-by-Region Projections\n(Cells Sorted by Top Projection)\n",
  "n = ", n_cells, " cells"
)
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
  key.title = "Normalized Projection Strength"
)
dev.copy(pdf, "Combined_ipsi-contra_projections_heatmap_top_region_sorted.pdf", width = 14, height = 10)
dev.off()

# Select ipsi columns and spinal cord columns from combined_norm_df
ipsi_cols <- grep("-ipsi", colnames(combined_norm_df), value = TRUE)
spinal_cols <- grep("_SP", colnames(combined_norm_df), value = TRUE)
cols_to_plot <- c(ipsi_cols, spinal_cols[spinal_cols %in% colnames(combined_norm_df)])
ipsi_matrix <- combined_norm_df[, cols_to_plot, drop = FALSE]
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
  "Combined Cell-by-Region Projections\n(Cells Sorted by Top Ipsi Projection)\n",
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
  key.title = "Normalized Projection Strength"
)
dev.copy(pdf, "Combined_ipsi-only_projections_heatmap_top_region_sorted.pdf", width = 10, height = 8)
dev.off()

################################## Regional grouping and probability analysis #########################################
source_matrix <- combined_norm  # This is cell × aggregated region with normalized projection strengths
head(source_matrix)

###################################################### Basic statistics: cell and region counts ######################################################
# Binarize combined_df: treat any non-zero normalized projection as "yes, this cell projects to this aggregated region".
# Number of aggregated regions each cell projects to (presence/absence, not strength)
num_regions_per_cell <- rowSums(source_matrix > 0)

# Number of cells projecting to each aggregated region (presence/absence, not strength)
num_cells_per_region <- colSums(source_matrix > 0)

# Plot bar plot for number of CELLS per aggregated region
# This shows which regions receive projections from the largest number of cells
par(mar = c(14, 4, 4, 2)) 
barplot(num_cells_per_region, las = 2, col = "lightblue", 
        xlab = "", ylab = "Number of cells", main = "Number of Cells per Aggregated Region")
mtext("Aggregated Region", side = 1, line = 12) 
par(mar = c(5, 4, 4, 2))  
dev.copy(pdf, "Combined_absolute_number_of_cells_with_non-zero_projections_to_aggregated_ROI.pdf", width = 12, height = 7)
dev.off()

# Region level stats for results section
region_min   <- min(num_cells_per_region)
region_max   <- max(num_cells_per_region)
region_median<- median(num_cells_per_region)

c(region_min = region_min, region_max = region_max, region_median = region_median)

# Cell levels stas for results section
cell_median <- median(num_regions_per_cell)
cell_iqr    <- quantile(num_regions_per_cell, probs = c(0.25, 0.75), names = FALSE)
cell_range  <- range(num_regions_per_cell)

list(
  median = cell_median,
  IQR_25 = cell_iqr[1],
  IQR_75 = cell_iqr[2],
  min = cell_range[1],
  max = cell_range[2]
)

###################################################### Co-innervation matrix (Jaccard index) → overlap among cells that hit at least one of the pair ######################################################
# Ignore projection strength and work with binary presence/absence:
# a grouped region is "innervated" by a cell if the SUMMED raw count for that group is > 0.
# For each pair of grouped regions (i, j), we compute the Jaccard index across cells:
#   J(i, j) = (# cells that project to BOTH i and j) / (# cells that project to i OR j) - Fraction of cells projecting to both regions, relative to those projecting to either
# This measures how often two grouped regions are co-innervated by the same cells.
co_innervation_matrix <- matrix(0, nrow = ncol(source_matrix), ncol = ncol(source_matrix), 
                                dimnames = list(colnames(source_matrix), colnames(source_matrix)))
for (i in seq_along(colnames(source_matrix))) {
  for (j in seq_along(colnames(source_matrix))) {
    proj_i <- source_matrix[, i] > 0
    proj_j <- source_matrix[, j] > 0
    co_innervation_matrix[i, j] <- sum(proj_i & proj_j) / sum(proj_i | proj_j)
  }
}
# Handle NaN (e.g., if no cells project to either region)
co_innervation_matrix[is.nan(co_innervation_matrix)] <- 0
co_innervation_matrix_masked <- co_innervation_matrix
diag(co_innervation_matrix_masked) <- NA  # Mask diagonal for visualization

# Heatmap of co-innervation (Jaccard index)
heatmap.2(
  co_innervation_matrix_masked,
  col = viridis(256),
  na.color = "grey",  # NA values (diagonal) appear grey
  trace = "none",
  dendrogram = "none",
  scale = "none",
  main = "Normalized Co-Innervation Fraction",
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
dev.copy(pdf, "Combined_co-innervation_fraction.pdf", width = 12, height = 10)
dev.off()

# Statistical summaries for the methods section
## Basic counts
n_cells   <- nrow(source_matrix)
n_targets <- ncol(source_matrix)

## How many cells target each target (binary presence/absence)
cells_per_target <- colSums(source_matrix > 0)

## How many targets each cell hits (binary)
targets_per_cell <- rowSums(source_matrix > 0)
n_multi_proj     <- sum(targets_per_cell > 1)
frac_multi_proj  <- n_multi_proj / n_cells

## Quick summaries 
list(
  n_cells = n_cells,
  n_targets = n_targets,
  n_multi_proj = n_multi_proj,
  frac_multi_proj = frac_multi_proj,
  cells_per_target_summary = summary(cells_per_target),
  targets_per_cell_summary = summary(targets_per_cell)
)

get_top_jaccard_pairs <- function(source_matrix, co_innervation_matrix, top_n = 20) {
  stopifnot(identical(colnames(source_matrix), colnames(co_innervation_matrix)))
  m <- co_innervation_matrix
  diag(m) <- NA
  
  ut <- which(upper.tri(m) & !is.na(m), arr.ind = TRUE)
  vals <- m[ut]
  ord <- order(vals, decreasing = TRUE)
  ut <- ut[ord, , drop = FALSE]
  vals <- vals[ord]
  
  top <- head(seq_along(vals), top_n)
  
  out <- lapply(top, function(k) {
    i <- ut[k, 1]; j <- ut[k, 2]
    ti <- colnames(m)[i]; tj <- colnames(m)[j]
    a <- source_matrix[, i] > 0
    b <- source_matrix[, j] > 0
    inter <- sum(a & b)
    uni   <- sum(a | b)
    data.frame(
      target_i = ti,
      target_j = tj,
      jaccard  = vals[k],
      n_intersection = inter,
      n_union        = uni,
      n_i = sum(a),
      n_j = sum(b)
    )
  })
  
  do.call(rbind, out)
}

top_pairs <- get_top_jaccard_pairs(source_matrix, co_innervation_matrix, top_n = 20)
top_pairs
write.csv(top_pairs, "top_jaccard_pairs_with_counts.csv", row.names = FALSE)

# Off-diagonal distribution summary
m <- co_innervation_matrix
diag(m) <- NA
offdiag <- m[upper.tri(m)]
offdiag <- offdiag[!is.na(offdiag) & is.finite(offdiag)]

offdiag_summary <- list(
  n_targets = ncol(source_matrix),
  n_pairs = length(offdiag),
  summary = summary(offdiag), # Min/1Q/Median/Mean/3Q/Max
  quantiles = quantile(offdiag, c(.05,.1,.25,.5,.75,.9,.95), na.rm=TRUE),
  frac_ge_05 = mean(offdiag >= 0.5),
  frac_ge_07 = mean(offdiag >= 0.7),
  frac_ge_08 = mean(offdiag >= 0.8)
)
offdiag_summary

# Blockwise summaries (supports “ipsi–ipsi blocks” etc.)
# Blockwise summaries (ipsi/contra only; excludes "other" targets such as SP)
type <- ifelse(grepl("-ipsi$", colnames(source_matrix)), "ipsi",
               ifelse(grepl("-contra$", colnames(source_matrix)), "contra", "other"))

# Keep only ipsi/contra targets for hemisphere-class summaries
keep <- type %in% c("ipsi","contra")
m2 <- m[keep, keep, drop = FALSE]
type2 <- type[keep]

# Helper that always returns UNIQUE pairs (upper triangle), for all block types
get_block_vals_ut <- function(m2, type2, block = c("ipsi-ipsi","contra-contra","ipsi-contra")) {
  block <- match.arg(block)
  ut <- upper.tri(m2) & is.finite(m2) & !is.na(m2)
  
  if (block == "ipsi-ipsi") {
    sel <- outer(type2=="ipsi", type2=="ipsi", "&")
  } else if (block == "contra-contra") {
    sel <- outer(type2=="contra", type2=="contra", "&")
  } else { # ipsi-contra
    sel <- outer(type2=="ipsi", type2=="contra", "&") | outer(type2=="contra", type2=="ipsi", "&")
  }
  
  v <- m2[ut & sel]
  v[is.finite(v) & !is.na(v)]
}

vals_ipsi_ipsi     <- get_block_vals_ut(m2, type2, "ipsi-ipsi")
vals_contra_contra <- get_block_vals_ut(m2, type2, "contra-contra")
vals_ipsi_contra   <- get_block_vals_ut(m2, type2, "ipsi-contra")

# Quick sanity checks: unique pair counts
c(
  n_ipsi_ipsi = length(vals_ipsi_ipsi),
  n_ipsi_contra = length(vals_ipsi_contra),
  n_contra_contra = length(vals_contra_contra),
  n_total_pairs = length(vals_ipsi_ipsi) + length(vals_ipsi_contra) + length(vals_contra_contra)
)

block_df <- rbind(
  data.frame(block="ipsi-ipsi", t(summary(vals_ipsi_ipsi))),
  data.frame(block="contra-contra", t(summary(vals_contra_contra))),
  data.frame(block="ipsi-contra", t(summary(vals_ipsi_contra)))
)
block_df

#Report top pairs in a results-friendly way (you already have these, but here’s a formatting helper)
format_top_pairs <- function(top_pairs, n = 5, digits = 3) {
  tp <- head(top_pairs, n)
  apply(tp, 1, function(r) {
    sprintf("%s–%s (J=%.3f; shared %s/%s; n_i=%s, n_j=%s)",
            r["target_i"], r["target_j"],
            as.numeric(r["jaccard"]),
            r["n_intersection"], r["n_union"],
            r["n_i"], r["n_j"])
  })
}
format_top_pairs(top_pairs, n = 6)

hist(offdiag, breaks=50)
c(median=median(offdiag), mean=mean(offdiag), p90=quantile(offdiag,0.9), p10=quantile(offdiag,0.1))

#Exact IQR and range (if you want them explicitly in text)
offdiag_iqr <- quantile(offdiag, c(0.25, 0.75), na.rm=TRUE)
offdiag_range <- range(offdiag, na.rm=TRUE)
offdiag_iqr
offdiag_range

#Add “how many pairs” exceed thresholds (counts, not just fractions)
c(
  n_pairs = length(offdiag),
  n_ge_05 = sum(offdiag >= 0.5),
  n_ge_07 = sum(offdiag >= 0.7),
  n_ge_08 = sum(offdiag >= 0.8)
)

#If you want to report blockwise counts (how many values in each block)
c(
  n_ipsi_ipsi = length(vals_ipsi_ipsi),
  n_ipsi_contra = length(vals_ipsi_contra),
  n_contra_contra = length(vals_contra_contra)
)

#If you want to cite blockwise 90th percentiles (often nice for “upper tail” comparisons)
c(
  ipsi_ipsi_p90 = unname(quantile(vals_ipsi_ipsi, 0.90)),
  ipsi_contra_p90 = unname(quantile(vals_ipsi_contra, 0.90)),
  contra_contra_p90 = unname(quantile(vals_contra_contra, 0.90))
)

#Optional null model scaffold
#This preserves per-cell projection breadth (row sums) while breaking target identity; you can compare observed off-diagonal median or the fraction ≥0.7 to the null distribution.
set.seed(1)
df_bin <- (source_matrix > 0) * 1

permute_within_rows <- function(mat_bin) {
  t(apply(mat_bin, 1, sample))
}

compute_offdiag_jaccard <- function(mat_bin) {
  p <- ncol(mat_bin)
  m <- matrix(0, p, p)
  for(i in 1:p) for(j in 1:p) {
    a <- mat_bin[, i] == 1
    b <- mat_bin[, j] == 1
    u <- sum(a | b)
    m[i,j] <- ifelse(u==0, NA, sum(a & b)/u)
  }
  diag(m) <- NA
  v <- m[upper.tri(m)]
  v[!is.na(v)]
}

obs_offdiag <- compute_offdiag_jaccard(df_bin)
null_stats <- replicate(200, {
  perm <- permute_within_rows(df_bin)
  v <- compute_offdiag_jaccard(perm)
  c(median=median(v), frac_ge_07=mean(v>=0.7))
})
apply(null_stats, 1, quantile, c(.05,.5,.95))

null_df <- as.data.frame(t(null_stats))
head(null_df)
summary(null_df)

# Observed values
obs_median     <- median(obs_offdiag, na.rm = TRUE)
obs_frac_ge_07 <- mean(obs_offdiag >= 0.7, na.rm = TRUE)

# Null quantiles (e.g., 5/50/95%)
null_median_q <- quantile(null_stats["median", ], c(0.05, 0.5, 0.95), na.rm = TRUE)
null_frac07_q <- quantile(null_stats["frac_ge_07", ], c(0.05, 0.5, 0.95), na.rm = TRUE)

list(
  observed = c(median = obs_median, frac_ge_07 = obs_frac_ge_07),
  null_median_quantiles = null_median_q,
  null_frac_ge_07_quantiles = null_frac07_q
)


###################################################### Projection fraction per region  ######################################################
# Restrict to multi-projecting cells: cells that project to >1 grouped region.
# Binarize combined_df (combined_df > 0), focus on presence/absence of projections, not projection strength.
# For each grouped region compute: fraction of multi-projecting cells that project to that region at all, otherwise asking which regions are most targeted by multi-projecting cells?
# Number of aggregated regions each cell projects to (binary: >0 counts)
num_regions_per_cell <- rowSums(source_matrix > 0)
# Indices of multi-projecting cells (cells with projections to >1 aggregated region)
multi_proj_cells <- which(num_regions_per_cell > 1)
# For each aggregated region, count how many multi-projecting cells project there (binary presence)
num_multi_proj_per_region <- colSums((source_matrix > 0)[multi_proj_cells, , drop = FALSE])
# Fraction of multi-projecting cells that target each aggregated region
projection_fraction <- num_multi_proj_per_region / length(multi_proj_cells)
# Create a data frame for ggplot
projection_fraction_df <- data.frame(
  Region = names(projection_fraction),
  ProjectionFraction = projection_fraction
)
# Ensure Region retains the order from source_matrix
projection_fraction_df$Region <- factor(projection_fraction_df$Region, levels = names(projection_fraction))

# This plot shows, for each region, the fraction of multi-projecting cells that target it
ggplot(projection_fraction_df, aes(x = Region, y = ProjectionFraction)) +
  geom_bar(stat = "identity", fill = "lightblue", alpha = 0.7, color = "black") +
  labs(title = "Fraction of Multi-Projecting Cells Targeting Each Aggregated Region",
       x = "Aggregated Region",
       y = "Fraction of Multi-Projecting Cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_text(margin = margin(t = 10)))

dev.copy(pdf, "Combined_projection_fraction_of_multi-projecting_cells.pdf", width = 10, height = 8)
dev.off()

###################################################### Conditional probability matrix ######################################################
# Goal: For each pair of grouped regions (A,B), compute
# P(cell projects to B | cell projects to A)
# using CELL-BASED, BINARY data (presence/absence of any projection), not projection strength.
# Fraction of cells that project to other regions given that they project to a specific region - If a cell projects to A, how likely to also project to B?
# Convert source_matrix to binary 
df_bin <- (source_matrix > 0) * 1

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

# Clean up
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
dev.copy(pdf, "Combined_conditional_probability_between_aggregated_regions.pdf", width = 12, height = 10)
dev.off()

##############################################
# Stats calculations for conditional probabilities (RESULTS)
# Uses directed A→B pairs (A ≠ B), plus support-filtered top pairs,
# per-target summaries, and lift (base-rate corrected enrichment).
##############################################

# --- Global summaries of conditional probabilities (ALL off-diagonal, directed) ---
P <- prob_other_given_here
stopifnot(is.matrix(P))
stopifnot(identical(rownames(P), colnames(P)))

diag(P) <- NA

# Vectorize ALL A→B entries excluding diagonal
p_off <- as.vector(P)
p_off <- p_off[!is.na(p_off) & is.finite(p_off)]

# Sanity check: with 45 targets, should be 45*44 = 1980 directed pairs
stopifnot(length(p_off) == ncol(P) * (ncol(P) - 1))

condprob_summary <- list(
  n_targets = ncol(P),
  n_pairs   = length(p_off),
  summary   = summary(p_off),  # Min/1Q/Median/Mean/3Q/Max
  quantiles = quantile(p_off, c(.05,.10,.25,.50,.75,.90,.95), na.rm=TRUE),
  range     = range(p_off, na.rm=TRUE),
  n_ge_05   = sum(p_off >= 0.5, na.rm=TRUE),
  n_ge_07   = sum(p_off >= 0.7, na.rm=TRUE),
  n_ge_08   = sum(p_off >= 0.8, na.rm=TRUE),
  frac_ge_05 = mean(p_off >= 0.5, na.rm=TRUE),
  frac_ge_07 = mean(p_off >= 0.7, na.rm=TRUE),
  frac_ge_08 = mean(p_off >= 0.8, na.rm=TRUE)
)
condprob_summary


# --- Top directional conditional-probability pairs with support ---
get_top_conditional_pairs <- function(df_bin, P, top_n = 20,
                                      min_nA = 30, min_nAB = 20) {
  stopifnot(identical(colnames(df_bin), colnames(P)))
  stopifnot(identical(rownames(P), colnames(P)))
  
  # counts
  nA  <- colSums(df_bin)          # #cells with A
  nAB <- t(df_bin) %*% df_bin     # joint counts (#cells with A and B)
  
  M <- P
  diag(M) <- NA
  ij <- which(!is.na(M) & is.finite(M), arr.ind = TRUE)
  
  out <- data.frame(
    A = rownames(M)[ij[,1]],
    B = colnames(M)[ij[,2]],
    P_B_given_A = M[ij],
    n_A  = nA[rownames(M)[ij[,1]]],
    n_B  = nA[colnames(M)[ij[,2]]],
    n_AB = as.integer(nAB[cbind(ij[,1], ij[,2])]),
    stringsAsFactors = FALSE
  )
  
  # support filters
  out <- out[out$n_A >= min_nA & out$n_AB >= min_nAB, , drop = FALSE]
  out <- out[order(out$P_B_given_A, decreasing = TRUE), , drop = FALSE]
  head(out, top_n)
}

top_cond_pairs <- get_top_conditional_pairs(df_bin, P,
                                            top_n = 25, min_nA = 30, min_nAB = 20)
top_cond_pairs

format_top_cond_pairs <- function(df, n = 8, digits = 3) {
  df <- head(df, n)
  apply(df, 1, function(r) {
    sprintf("%s→%s: P=%.3f (n_A=%s, n_AB=%s; n_B=%s)",
            r["A"], r["B"], as.numeric(r["P_B_given_A"]),
            r["n_A"], r["n_AB"], r["n_B"])
  })
}
format_top_cond_pairs(top_cond_pairs, n = 10)


# --- For each A, list its top-K partners B by P(B|A) ---
top_partners_by_A <- function(P, df_bin, k = 5, min_nA = 30) {
  nA  <- colSums(df_bin)
  nAB <- t(df_bin) %*% df_bin
  
  res <- lapply(rownames(P), function(A) {
    if (nA[A] < min_nA) return(NULL)
    v <- P[A, ]
    v[A] <- NA
    ok <- which(!is.na(v) & is.finite(v))
    if (!length(ok)) return(NULL)
    ord <- ok[order(v[ok], decreasing = TRUE)]
    ord <- head(ord, k)
    data.frame(
      A = A,
      n_A = nA[A],
      B = colnames(P)[ord],
      P_B_given_A = v[ord],
      n_AB = as.integer(nAB[A, colnames(P)[ord]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}

top_by_A <- top_partners_by_A(P, df_bin, k = 5, min_nA = 30)
head(top_by_A, 20)


# --- Per-target summaries: "does each A have at least one strong partner" and "how many moderate partners" ---
A_selectivity <- function(P, df_bin, min_nA = 30, thr = 0.5) {
  nA <- colSums(df_bin)
  out <- lapply(rownames(P), function(A) {
    if (nA[A] < min_nA) return(NULL)
    v <- P[A, ]
    v[A] <- NA
    v <- v[!is.na(v) & is.finite(v)]
    data.frame(
      A = A,
      n_A = nA[A],
      maxP = max(v),
      medianP = median(v),
      n_partners_ge_thr = sum(v >= thr),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

A_sel_df <- A_selectivity(P, df_bin, min_nA = 30, thr = 0.5)
A_sel_df[order(A_sel_df$maxP, decreasing = TRUE), ]

row_summary <- list(
  n_A_rows = nrow(A_sel_df),
  maxP_summary = summary(A_sel_df$maxP),
  medianP_summary = summary(A_sel_df$medianP),
  partners_ge_05_summary = summary(A_sel_df$n_partners_ge_thr)
)
row_summary


# --- Lift (base-rate corrected enrichment): lift(A→B) = P(B|A) / P(B) ---
n_cells <- nrow(df_bin)
pB <- colSums(df_bin) / n_cells
pB[pB == 0] <- NA  # safety; should rarely/never happen with real targets

lift_mat <- P
diag(lift_mat) <- NA
lift_mat <- sweep(lift_mat, 2, pB, FUN = "/")

get_top_lift_pairs <- function(df_bin, P, lift_mat, pB,
                               top_n = 20, min_nA = 30, min_nAB = 20) {
  nA  <- colSums(df_bin)
  nAB <- t(df_bin) %*% df_bin
  
  ij <- which(!is.na(lift_mat) & is.finite(lift_mat), arr.ind = TRUE)
  
  out <- data.frame(
    A = rownames(lift_mat)[ij[,1]],
    B = colnames(lift_mat)[ij[,2]],
    lift = lift_mat[ij],
    P_B_given_A = P[ij],
    p_B = pB[colnames(lift_mat)[ij[,2]]],
    n_A = nA[rownames(lift_mat)[ij[,1]]],
    n_AB = as.integer(nAB[cbind(ij[,1], ij[,2])]),
    stringsAsFactors = FALSE
  )
  
  out <- out[out$n_A >= min_nA & out$n_AB >= min_nAB, , drop = FALSE]
  out <- out[order(out$lift, decreasing = TRUE), , drop = FALSE]
  head(out, top_n)
}

top_lift_pairs <- get_top_lift_pairs(df_bin, P, lift_mat, pB,
                                     top_n = 20, min_nA = 30, min_nAB = 20)
top_lift_pairs

###################################################### Co-occurrence probability matrix (intersection/all) → absolute prevalence of a co-projection pattern in the whole population ######################################################
# Goal: For each pair of grouped regions (i, j), compute
#   P(cell projects to BOTH i AND j)
# i.e. (# cells that project to both i and j) / (# all cells).
# This is a CELL-BASED, BINARY measure (presence/absence),
# and it treats each hemisphere-specific grouped region (e.g. striatum_contra, striatum_ipsi)
# as a separate node.

# Binary presence/absence for ALL grouped regions (including spinal cord)
brain_regions <- colnames(source_matrix)  
brain_binary <- (source_matrix[, brain_regions] > 0) * 1

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

# Plot the co-occurrence matrix (assuming it's already in the desired order)
co_occur_prob_masked <- co_occur_prob
diag(co_occur_prob_masked) <- NA # hide P(i AND i) = P(i)
heatmap.2(
  co_occur_prob_masked,
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

