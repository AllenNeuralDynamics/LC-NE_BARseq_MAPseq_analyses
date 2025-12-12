#set working directory 
setwd('/scratch/BARseq_780345/')

################################################## load and pre-process relevant files ######################################################
# coordinates of all possible cells for plotting, only keep cells with at least one gene count
# load files with all cells and coordinates, exclude low quality cells with zero total gene counts
barseq <- readRDS("combined_neurons_clust_CCFv2_uid_cpm_log.rds")
counts_mat <- assay(barseq, "counts")
# Sum total reads per cell
total_reads_per_cell <- colSums(counts_mat)
reads_distribution <- table(total_reads_per_cell)
# Plot the distribution
barplot(reads_distribution, 
        main = "Distribution of total reads per cell", 
        xlab = "Total reads per cell", 
        ylab = "Number of cells")
# Identify cells with at least one read, subset the SingleCellExperiment object to keep only those cells
cells_to_keep <- total_reads_per_cell > 0
barseq <- barseq[, cells_to_keep]


# load transcriptomic data for LC-NE neurons and check gene counts between barcoded and non-barcoded cells
LCNE_barseq  <- readRDS("LCNE_clusters_filtered_coherence_filtered_cpm_log_clust.rds")
counts_data <- assay(LCNE_barseq, "counts")
total_counts <- Matrix::colSums(counts_data)  # sum counts for each cell (column)
barcode_info <- colData(LCNE_barseq)$barcode
barcoded_counts <- total_counts[barcode_info == 1]
non_barcoded_counts <- total_counts[barcode_info == 0]
hist(non_barcoded_counts, col = "palegreen3", breaks = 50,
     xlim = c(0,150),
     main = "Total Gene Counts per Cell", xlab = "Total Counts", ylab = "Number of Cells")
hist(barcoded_counts, col = "pink", breaks = 50, add = TRUE)
legend("topright", legend = c("Non-barcoded", "Barcoded"), fill = c("palegreen3", "pink"))
dev.copy(pdf, file = "BC_nonBC_cells_total_gene_counts.pdf", width = 8, height = 8)
dev.off()

##########################################################################################################################################
# Check PCA relationship between gene expression and cell position in space
# Prepare gene expression data and filter genes with zero variance
expr <- assay(LCNE_barseq, "logcounts")  # or "counts" 
gene_var <- apply(expr, 1, var)
expr_filtered <- expr[gene_var > 0, ]
expr_t_filtered <- t(as.matrix(expr_filtered))  # cells x genes

# PCA on gene expression
pca <- prcomp(expr_t_filtered, scale. = TRUE)

# Prepare spatial coordinates and remove NAs
ccf_coords <- colData(LCNE_barseq)[, c("CCF_ML", "CCF_DV", "CCF_AP")]
df <- data.frame(
  PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], PC4 = pca$x[,4], PC5 = pca$x[,5],PC6 = pca$x[,6],
  CCF_ML = ccf_coords$CCF_ML, CCF_DV = ccf_coords$CCF_DV, CCF_AP = ccf_coords$CCF_AP
)
df <- na.omit(df)  # Remove rows with any NA

# PCA on spatial coordinates (3D)
spatial_coords <- df[, c("CCF_ML", "CCF_AP", "CCF_DV")]
spatial_pca <- prcomp(spatial_coords, scale. = TRUE)
df$spatial_PC1 <- spatial_pca$x[, 1]

# Correlate gene PCs with spatial PC1 and plot
plots <- list()
for (i in 1:6) {
  pc_name <- paste0("PC", i)
  fit <- lm(df[[pc_name]] ~ df$spatial_PC1)
  r2 <- summary(fit)$r.squared
  p <- ggplot(df, aes_string(x = "spatial_PC1", y = pc_name)) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", col = "red") +
    labs(x = "Spatial PC1", y = pc_name,
         title = paste0(pc_name, " vs Spatial PC1\nR² = ", round(r2, 3))) +
    theme_minimal()
  plots[[i]] <- p
}
grid.arrange(grobs = plots, nrow = 2, ncol = 3)
dev.copy(pdf, file = "PCs_vs_spatialPC1_correlations.pdf", width = 10, height = 7)
dev.off()

# Plot cell positions colored by spatial PC1 (coronal and sagittal)
p1 <- ggplot(df, aes(x = CCF_ML, y = CCF_DV, color = spatial_PC1)) +
  geom_point(alpha = 0.7) +
  scale_y_reverse() +
  scale_color_viridis_c() +
  labs(x = "ML", y = "DV", color = "Spatial PC1", title = "Coronal (ML vs DV)") +
  theme_minimal()
p2 <- ggplot(df, aes(x = CCF_AP, y = CCF_DV, color = spatial_PC1)) +
  geom_point(alpha = 0.7) +
  scale_y_reverse() +
  scale_color_viridis_c() +
  labs(x = "AP", y = "DV", color = "Spatial PC1", title = "Sagittal (AP vs DV)") +
  theme_minimal()
grid.arrange(p1, p2, nrow = 1)
dev.copy(pdf, file = "spatialPC1_coronal_sagittal.pdf", width = 10, height = 7)
dev.off()

plots <- list()
ccf_names <- c("CCF_ML", "CCF_DV", "CCF_AP")
for (i in 1:3) {
  for (j in 1:5) {
    p <- ggplot(df, aes_string(x = paste0("PC", j), y = ccf_names[i])) +
      geom_point(alpha = 0.5, size = 0.5) +
      geom_smooth(method = "lm", col = "red") +
      labs(x = paste0("PC", j), y = ccf_names[i]) +
      theme_minimal()
    plots[[length(plots) + 1]] <- p
  }
}
# Arrange plots in 3x5 grid
grid.arrange(grobs = plots, nrow = 3, ncol = 5)
dev.copy(pdf, file = "PCs_vs_CCF_axes.pdf", width = 12, height = 8)
dev.off()

##########################################################################################################################################
NMF_types <- readr::read_csv("./NMF_ids_factors.csv", col_names = TRUE)

#subset LC_barseq data to only include barcoded LC neurons
matches <- match(NMF_types$split_cellID, colnames(LCNE_barseq))
LCNE_barcoded <- LCNE_barseq[, matches]

#add cluster assignment information but first check that order is the same
identical(colnames(LCNE_barcoded), as.character(NMF_types$split_cellID)) #must be TRUE
colData(LCNE_barcoded)$NMF_factor <- NMF_types$factor
colData(LCNE_barcoded)$proj_pattern <- NMF_types$proj_pattern
colData(LCNE_barcoded)$proj_target <- NMF_types$proj_target
table(LCNE_barcoded$NMF_factor)
table(LCNE_barcoded$proj_pattern)
table(LCNE_barcoded$proj_target)
table(LCNE_barcoded$louvain_cluster)

#plot soma locations to visualize NMF_type and cluster topography
# Convert to data frame
barseq_df <- as.data.frame(barseq@colData)
LCNE_barseq_df <- as.data.frame(LCNE_barseq@colData)
LCNE_barcoded_df <- as.data.frame(LCNE_barcoded@colData)

# Identify unique slices in the data
print(unique(barseq_df$slice))
print(unique(LCNE_barseq_df$slice))
print(unique(LCNE_barcoded_df$slice))

################################################## plot in NMF cells on slices ######################################################
outdir <- "NMF_plots"
if (!dir.exists(outdir)) dir.create(outdir)
# Plot slices containing only the barcoded cells
#plot soma locations with a figure legend
slices_barcoded <- unique(LCNE_barcoded_df$slice)
# Loop over each slice
for (s in slices_barcoded) {
  # Subset all neurons in this slice for background
  barseq_in_slice <- subset(barseq_df, slice == s)
  # Subset LCNE neurons in this slice for cluster plotting
  LCNE_barcoded_in_slice <- subset(LCNE_barcoded_df, slice == s)
  
  # Create the plot
  p <- ggplot() +
    # Plot all neurons in grey
    geom_point(data = barseq_in_slice, aes(x = CCF_ML, y = CCF_DV), color = "grey", size = 0.1) +
    # Plot LCNE neurons colored by louvain cluster
    geom_point(data = LCNE_barcoded_in_slice, aes(x = CCF_ML, y = CCF_DV, color = as.factor(NMF_factor)), size = 0.8, alpha = 0.9) +
    coord_fixed() +
    ggtitle("LC clusters") +
    theme_minimal() +
    scale_y_reverse() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    # Add a color legend
    scale_color_manual(values = c("#469990", "#f032e6", "#3cb44b", "#42d4f4", "#911eb4"),
                       name = "NMF Factor",
                       breaks = c(1, 2, 3, 4, 5), 
                       labels = c("F1", "F2", "F3", "F4", "F5")) 
  plot(p)
  # Save the plot
  ggsave(filename = file.path(outdir, paste0(s, "_NMF.pdf")), plot = p, width = 6, height = 6, units = "in")
}

################################################## center of mass for NMF groups ######################################################
# Define midline
midline_ML <- 228

# Reflect all cells to one hemisphere (e.g., reflect left to right)
LCNE_barcoded_reflected <- LCNE_barcoded_df %>%
  mutate(
    # Determine which hemisphere each cell is in
    hemisphere = ifelse(CCF_ML < midline_ML, "left", "right"),
    # Reflect left hemisphere to right hemisphere
    CCF_ML_reflected = ifelse(CCF_ML < midline_ML, 
                              2 * midline_ML - CCF_ML,  # reflect across midline
                              CCF_ML)  # keep right hemisphere as is
  )
# Calculate center of mass for each NMF factor using reflected coordinates
center_of_mass_reflected <- LCNE_barcoded_reflected %>%
  group_by(NMF_factor) %>%
  summarise(
    center_ML = mean(CCF_ML_reflected, na.rm = TRUE),
    center_DV = mean(CCF_DV, na.rm = TRUE),
    center_AP = mean(CCF_AP, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  )
print(center_of_mass_reflected)

# Coronal view (ML vs DV)
p_coronal <- ggplot() +
  # Plot all reflected cells
  geom_point(data = LCNE_barcoded_reflected, 
             aes(x = CCF_ML_reflected, y = CCF_DV, color = as.factor(NMF_factor)), 
             size = 0.7, alpha = 0.4) +
  # Plot centers of mass
  geom_point(data = center_of_mass_reflected, 
             aes(x = center_ML, y = center_DV, color = as.factor(NMF_factor)), 
             size = 4, alpha = 0.8, shape = 17) +
  # Add midline reference
  geom_vline(xintercept = midline_ML, linetype = "dashed", color = "black") +
  coord_fixed() +
  ggtitle("NMF Factor Centers of Mass - Coronal View (ML vs DV)") +
  theme_minimal() +
  scale_y_reverse() +
  scale_color_manual(values = c("#469990", "#f032e6", "#3cb44b", "#42d4f4", "#911eb4"),
                     name = "NMF Factor",
                     breaks = c(1, 2, 3, 4, 5), 
                     labels = c("F1", "F2", "F3", "F4", "F5")) +
  labs(x = "CCF_ML (reflected)", y = "CCF_DV")
print(p_coronal)
ggsave("NMF_centers_coronal_reflected.pdf", width = 8, height = 8)

# Sagittal view (AP vs DV)
p_sagittal <- ggplot() +
  # Plot all reflected cells
  geom_point(data = LCNE_barcoded_reflected, 
             aes(x = CCF_AP, y = CCF_DV, color = as.factor(NMF_factor)), 
             size = 0.5, alpha = 0.3) +
  # Plot centers of mass
  geom_point(data = center_of_mass_reflected, 
             aes(x = center_AP, y = center_DV, color = as.factor(NMF_factor)), 
             size = 4, shape = 17) +
  coord_fixed() +
  ggtitle("NMF Factor Centers of Mass - Sagittal View (AP vs DV)") +
  theme_minimal() +
  scale_y_reverse() +
  scale_color_manual(values = c("#469990", "#f032e6", "#3cb44b", "#42d4f4", "#911eb4"),
                     name = "NMF Factor",
                     breaks = c(1, 2, 3, 4, 5), 
                     labels = c("F1", "F2", "F3", "F4", "F5")) +
  labs(x = "CCF_AP", y = "CCF_DV")

print(p_sagittal)
ggsave("NMF_centers_sagittal_reflected.pdf", width = 8, height = 8)

# Calculate pairwise distances between centers using reflected coordinates
centers_coords_reflected <- center_of_mass_reflected[, c("center_ML", "center_DV", "center_AP")]
rownames(centers_coords_reflected) <- paste0("F", center_of_mass_reflected$NMF_factor)
distance_matrix_reflected <- as.matrix(dist(centers_coords_reflected))
print(round(distance_matrix_reflected, 2))

################################################## plot in transcriptomic cluster cells on slices ######################################################
outdir <- "LCNE_cluster_plots"
if (!dir.exists(outdir)) dir.create(outdir)
#plot soma locations with a figure legend
# Identify unique slices in the LCNE_barseq data
slices_LCNE <- unique(LCNE_barseq_df$slice)
# Loop over each slice for louvain cluster plotting
for (s in slices_LCNE) {
  # Subset all neurons in this slice for background
  barseq_in_slice <- subset(barseq_df, slice == s)
  # Subset LCNE neurons in this slice for cluster plotting
  LCNE_in_slice <- subset(LCNE_barseq_df, slice == s)
  
  # Create the plot
  p <- ggplot() +
    # Plot all neurons in grey
    geom_point(data = barseq_in_slice, aes(x = CCF_ML, y = CCF_DV), color = "grey", size = 0.1) +
    # Plot LCNE neurons colored by louvain cluster
    geom_point(data = LCNE_in_slice, aes(x = CCF_ML, y = CCF_DV, color = as.factor(louvain_cluster)), size = 0.8, alpha = 0.9) +
    coord_fixed() +
    ggtitle("LC clusters") +
    theme_minimal() +
    scale_y_reverse() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    scale_color_manual(
      values = c("#dcbeff","#aaffc3", "#fabed4",  "#e6194B"),
      name = "cluster",
      breaks = c("1", "2", "3", "4"),
      labels = c("LC_1", "LC_2", "LC_3", "LC_4")
    )
  plot(p)
  # Save the plot
  ggsave(filename = file.path(outdir, paste0(s, "_louvain_LCNE_filtered.pdf")), plot = p, width = 6, height = 6, units = "in")
}

################################################## center of mass for transcriptomic groups ######################################################
# Reflect all LCNE cells to one hemisphere (e.g., reflect left to right)
LCNE_barseq_reflected <- LCNE_barseq_df %>%
  mutate(
    # Determine which hemisphere each cell is in
    hemisphere = ifelse(CCF_ML < midline_ML, "left", "right"),
    # Reflect left hemisphere to right hemisphere
    CCF_ML_reflected = ifelse(CCF_ML < midline_ML, 
                              2 * midline_ML - CCF_ML,  # reflect across midline
                              CCF_ML)  # keep right hemisphere as is
  )

# Calculate center of mass for each Louvain cluster using reflected coordinates
center_of_mass_louvain <- LCNE_barseq_reflected %>%
  group_by(louvain_cluster) %>%
  summarise(
    center_ML = mean(CCF_ML_reflected, na.rm = TRUE),
    center_DV = mean(CCF_DV, na.rm = TRUE),
    center_AP = mean(CCF_AP, na.rm = TRUE),
    n_cells = n(),
    .groups = 'drop'
  )
print(center_of_mass_louvain)

# Coronal view (ML vs DV)
p_coronal_louvain <- ggplot() +
  # Plot all reflected cells
  geom_point(data = LCNE_barseq_reflected, 
             aes(x = CCF_ML_reflected, y = CCF_DV, color = as.factor(louvain_cluster)), 
             size = 0.7, alpha = 0.4) +
  # Plot centers of mass
  geom_point(data = center_of_mass_louvain, 
             aes(x = center_ML, y = center_DV, color = as.factor(louvain_cluster)), 
             size = 4, alpha = 0.8, shape = 17) +
  # Add midline reference
  geom_vline(xintercept = midline_ML, linetype = "dashed", color = "black") +
  coord_fixed() +
  ggtitle("Louvain Cluster Centers of Mass - Coronal View (ML vs DV)") +
  theme_minimal() +
  scale_y_reverse() +
  scale_color_manual(
    values = c("#dcbeff","#aaffc3", "#fabed4",  "#e6194B"),
    name = "cluster",
    breaks = c("1", "2", "3", "4"),
    labels = c("LC_1", "LC_2", "LC_3", "LC_4")
  ) +
  labs(x = "CCF_ML (reflected)", y = "CCF_DV")

print(p_coronal_louvain)
ggsave("Louvain_centers_coronal_reflected.pdf", width = 8, height = 8)

# Sagittal view (AP vs DV)
p_sagittal_louvain <- ggplot() +
  # Plot all reflected cells
  geom_point(data = LCNE_barseq_reflected, 
             aes(x = CCF_AP, y = CCF_DV, color = as.factor(louvain_cluster)), 
             size = 0.5, alpha = 0.3) +
  # Plot centers of mass
  geom_point(data = center_of_mass_louvain, 
             aes(x = center_AP, y = center_DV, color = as.factor(louvain_cluster)), 
             size = 4, shape = 17) +
  coord_fixed() +
  ggtitle("Louvain Cluster Centers of Mass - Sagittal View (AP vs DV)") +
  theme_minimal() +
  scale_y_reverse() +
  scale_color_manual(
    values = c("#dcbeff","#aaffc3", "#fabed4",  "#e6194B"),
    name = "cluster",
    breaks = c("1", "2", "3", "4"),
    labels = c("LC_1", "LC_2", "LC_3", "LC_4")
  ) +
  labs(x = "CCF_AP", y = "CCF_DV")

print(p_sagittal_louvain)
ggsave("Louvain_centers_sagittal_reflected.pdf", width = 8, height = 8)

# Calculate pairwise distances between centers using reflected coordinates
centers_coords_louvain <- center_of_mass_louvain[, c("center_ML", "center_DV", "center_AP")]
rownames(centers_coords_louvain) <- paste0("LC_", center_of_mass_louvain$louvain_cluster)
distance_matrix_louvain <- as.matrix(dist(centers_coords_louvain))
print(round(distance_matrix_louvain, 2))

################################################## plot in 3D using CCF coordinates ######################################################
# whole brain plot to check LC registration positioning
# Assign LCNE neurons to LC group, the rest are other group 
LC_uids <- colData(LCNE_barseq)$uid
colData(barseq)$cluster_type <- ifelse(colData(barseq)$uid %in% LC_uids, "LC", "other")
table(barseq$cluster_type)

# Extract the 3D coordinates and the characteristic
x <- colData(barseq)$CCF_ML
y <- colData(barseq)$CCF_DV
z <- colData(barseq)$CCF_AP
characteristic <- colData(barseq)$cluster_type

# Create a named vector of colors and color vector based on the characteristic
color_map <- c("LC" = "#00FF00", "other" = "#696969")
colors <- color_map[characteristic]
# Create a named vector of sizes and size vector based on the characteristic
size_map <- c("LC" = 3, "other" = 1)
sizes <- size_map[characteristic]
# Create a named vector of alphas and alpha vector based on the characteristic
alpha_map <- c("LC" = 0.9, "other" = 0.1)
alphas <- alpha_map[characteristic]
# Initialize an empty plotly object
p <- plot_ly()
# Add a trace for the "LC" cells
p <- add_trace(p, 
               x = ~x[characteristic == "LC"], 
               y = ~y[characteristic == "LC"], 
               z = ~z[characteristic == "LC"], 
               name = "LC", 
               type = "scatter3d", 
               mode = "markers",
               marker = list(size = sizes[characteristic == "LC"],
                             color = colors[characteristic == "LC"],
                             opacity = alphas[characteristic == "LC"],
                             line = list(color = colors[characteristic == "LC"], width = 2)))
# Add a trace for the "other" cells
p <- add_trace(p, 
               x = ~x[characteristic == "other"], 
               y = ~y[characteristic == "other"], 
               z = ~z[characteristic == "other"], 
               name = "other", 
               type = "scatter3d", 
               mode = "markers",
               marker = list(size = sizes[characteristic == "other"],
                             color = colors[characteristic == "other"],
                             opacity = alphas[characteristic == "other"],
                             line = list(color = colors[characteristic == "other"], width = 2)))
# Save the plot as an HTML file
htmlwidgets::saveWidget(p, "whole_brain_LCgroup_green.html")
# Open the HTML file in a web browser
browseURL("whole_brain_LCgroup_green.html")

################################################## LC with cell and NFM types plots ##################################################
# add NMF info directly to LCNE_barseq
# match NMF entries to LCNE_barseq columns using split_cellID and column names
matches <- match(NMF_types$split_cellID, colnames(LCNE_barseq))
valid <- !is.na(matches)
stopifnot(all(colnames(LCNE_barseq)[matches[valid]] == NMF_types$split_cellID[valid]))

# initialise everything as "other"
colData(LCNE_barseq)$NMF_type     <- "other"
colData(LCNE_barseq)$proj_pattern <- "other"
colData(LCNE_barseq)$proj_target  <- "other"

# only use valid matches
valid <- !is.na(matches)

colData(LCNE_barseq)$NMF_type[matches[valid]]     <- as.character(NMF_types$factor[valid])
colData(LCNE_barseq)$proj_pattern[matches[valid]] <- NMF_types$proj_pattern[valid]
colData(LCNE_barseq)$proj_target[matches[valid]]  <- NMF_types$proj_target[valid]

# sanity checks
table(LCNE_barseq$NMF_type)
table(LCNE_barseq$proj_pattern)
table(LCNE_barseq$proj_target)
table(LCNE_barseq$louvain_cluster)

# Extract the 3D coordinates and the characteristic
x <- colData(LCNE_barseq)$CCF_ML
y <- colData(LCNE_barseq)$CCF_DV
z <- colData(LCNE_barseq)$CCF_AP
NMF <- colData(LCNE_barseq)$NMF_type
proj_pattern <- colData(LCNE_barseq)$proj_pattern
proj_target <- colData(LCNE_barseq)$proj_target
cluster_type <- colData(LCNE_barseq)$louvain_cluster

################################################## all NFM factors ##################################################
# Create a named vector of colors and color vector based on the characteristic
color_map_NMF <- c("1"="#469990", "2"="#f032e6", "3"="#3cb44b","4"="#42d4f4","5"="#911eb4","other"="lightgrey")
colors <- color_map_NMF[NMF]
# Create a named vector of sizes and size vector based on the characteristic
size_map_NMF <- c("1"=6, "2"=6, "3"=6,"4"=6,"5"=6,"other"=5)
sizes <- size_map_NMF[NMF]
# Create a named vector of alphas and alpha vector based on the characteristic
alpha_map_NMF <- c("1"=0.9, "2"=0.9, "3"=0.9,"4"=0.9,"5"=0.9, "other"=0.6)
alphas <- alpha_map_NMF[NMF]

# Initialize an empty plotly object
p <- plot_ly()
# Add a trace for each NMF type
for (type in unique(NMF)) {
  indices <- NMF == type
  df <- data.frame(x = x[indices], 
                   y = y[indices], 
                   z = z[indices], 
                   size = sizes[indices],
                   color = colors[indices],
                   alpha = alphas[indices])
  p <- add_trace(p, 
                 data = df,
                 x = ~x, 
                 y = ~y, 
                 z = ~z, 
                 name = type, 
                 type = "scatter3d", 
                 mode = "markers",
                 marker = list(size = ~size,
                               color = ~color,
                               opacity = ~alpha,
                               line = list(width = 0)))
}
# Save the plot as an HTML file
htmlwidgets::saveWidget(p, "LC_all_NMF_factors.html")
# Open the HTML file in a web browser
browseURL("LC_all_NMF_factors.html")

###################### projection patterns A-P #################################
color_map_proj_pattern <- c("Anterior"="blue", "Posterior"="red", "other"="lightgrey")
colors <- color_map_proj_pattern[proj_pattern]
size_map_proj_pattern <- c("Anterior"=6, "Posterior"=6, "other"=5)
sizes <- size_map_proj_pattern[proj_pattern]
alpha_map_proj_pattern <- c("Anterior"=0.9, "Posterior"=0.9, "other"=0.6)
alphas <- alpha_map_proj_pattern[proj_pattern]

# Initialize an empty plotly object
p <- plot_ly()
# Add a trace for each proj_pattern type
for (type in unique(proj_pattern)) {
  indices <- proj_pattern == type
  df <- data.frame(x = x[indices], 
                   y = y[indices], 
                   z = z[indices], 
                   size = sizes[indices],
                   color = colors[indices],
                   alpha = alphas[indices])
  p <- add_trace(p, 
                 data = df,
                 x = ~x, 
                 y = ~y, 
                 z = ~z, 
                 name = type, 
                 type = "scatter3d", 
                 mode = "markers",
                 marker = list(size = ~size,
                               color = ~color,
                               opacity = ~alpha,
                               line = list(width = 0)))
}
# Save the plot as an HTML file
htmlwidgets::saveWidget(p, "LC_NMF_A-P.html")
# Open the HTML file in a web browser
browseURL("LC_NMF_A-P.html")

####################### all BarSeq_cluster types ########################################
# Create a named vector of colors and color vector based on the characteristic
color_map_louvain_cluster <- c("1"="#dcbeff", "2"="#aaffc3", "3"="#fabed4","4"="#e6194B", "other"="lightgrey")
colors <- color_map_louvain_cluster[cluster_type]
# Create a named vector of sizes and size vector based on the characteristic
size_map_louvain_cluster <- c("1"=5, "2"=5, "3"=5,"4"=5,"other"=5)
sizes <- size_map_louvain_cluster[cluster_type]
# Create a named vector of alphas and alpha vector based on the characteristic
alpha_map_louvain_cluster <- c("1"=0.8, "2"=0.8, "3"=0.8,"4"=0.8,"other"=0.6)
alphas <- alpha_map_louvain_cluster[cluster_type]

# Initialize an empty plotly object
p <- plot_ly()
# Add a trace for each louvain cluster
for (type in unique(cluster_type)) {
  indices <- cluster_type == type
  df <- data.frame(x = x[indices], 
                   y = y[indices], 
                   z = z[indices], 
                   size = sizes[indices],
                   color = colors[indices],
                   alpha = alphas[indices])
  p <- add_trace(p, 
                 data = df,
                 x = ~x, 
                 y = ~y, 
                 z = ~z, 
                 name = type, 
                 type = "scatter3d", 
                 mode = "markers",
                 marker = list(size = ~size,
                               color = ~color,
                               opacity = ~alpha,
                               line = list(width = 0)))
}
# Save the plot as an HTML file
htmlwidgets::saveWidget(p, "LC_all_BarSeq_clusters.html")
# Open the HTML file in a web browser
browseURL("LC_all_BarSeq_clusters.html")
