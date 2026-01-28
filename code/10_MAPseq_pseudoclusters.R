## ============================================================
## 00) Setup
## ============================================================
OUT_DIR <- "/scratch/BARseq_780345-780346_combined/"
setwd(OUT_DIR)

## Helper
strip_suffix <- function(x) sub("\\.[0-9]+.*$", "", as.character(x))

## ============================================================
## 01) Load inputs
## ============================================================
df_top <- read_csv("cell_top_projections_with_coords.csv", show_col_types = FALSE)
df_nmf <- read_csv("NMF_ids_factors.csv", show_col_types = FALSE)
pc3    <- read_csv("imputed_pseudoclusters_brain3.csv", show_col_types = FALSE)
pc4    <- read_csv("imputed_pseudoclusters_brain4.csv", show_col_types = FALSE)

## ============================================================
## 02) Standardize keys + build final_df (drop ambiguous in-core keys)
## ============================================================

# --- Top projection table (506 rows)
df_top2 <- df_top %>%
  mutate(uid_key = strip_suffix(cell_id)) %>%
  select(uid_key, cell_id, top_projection, top_projection_strength, CCF_DV, CCF_ML, CCF_AP)

# --- NMF table (506 rows)
df_nmf2 <- df_nmf %>%
  mutate(uid_key = as.character(split_cellID)) %>%
  select(uid_key, cellID, factor, proj_pattern, proj_target, split_cellID)

# --- Pseudocluster tables
pc3u <- pc3 %>%
  transmute(uid_key = as.character(uid),
            imputed_pseudoclusters = as.numeric(imputed_pseudoclusters),
            CI_width = as.numeric(CI_width),
            brain = "brain3")

pc4u <- pc4 %>%
  transmute(uid_key = as.character(uid),
            imputed_pseudoclusters = as.numeric(imputed_pseudoclusters),
            CI_width = as.numeric(CI_width),
            brain = "brain4")

pc_all <- bind_rows(pc3u, pc4u) %>% as_tibble()

# Core merge: should be 506
core_506 <- inner_join(df_top2, df_nmf2, by = "uid_key")
cat("core_506 rows:", nrow(core_506), "\n")

# Find ambiguous uid_key (present in BOTH pc3 and pc4)
ambiguous_keys <- pc_all %>%
  dplyr::group_by(uid_key) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(uid_key)

cat("Ambiguous uid_key (pc3 vs pc4 overlap):", length(ambiguous_keys), "\n")

# Only drop ambiguous keys that are actually in the 506-cell set
ambig_in_core <- intersect(ambiguous_keys, core_506$uid_key)
cat("Ambiguous keys in core_506:", length(ambig_in_core), "\n")
if (length(ambig_in_core) > 0) print(ambig_in_core)

core_no_ambig <- core_506 %>% filter(!uid_key %in% ambig_in_core)

# Keep only unique pseudocluster rows (uid_key appears exactly once across pc3+pc4)
pc_unique <- pc_all %>%
  group_by(uid_key) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(uid_key, imputed_pseudoclusters, CI_width)  # drop brain

# Final merge
final_df <- left_join(core_no_ambig, pc_unique, by = "uid_key")

cat("final_df rows:", nrow(final_df), "\n")
cat("duplicated uid_key:", sum(duplicated(final_df$uid_key)), "\n")
cat("missing pseudoclusters:", sum(is.na(final_df$imputed_pseudoclusters)), "\n")

write_csv(final_df, "pseudocluster_df_for_plotting_and_regression.csv")

## ============================================================
## 03) Define analysis datasets
##     df_stats: all groups (no dropping)
##     df_plot : min_n filter for visualization only
## ============================================================

df_stats <- final_df %>%
  mutate(
    top_projection = as.factor(top_projection),
    proj_target    = as.factor(proj_target),
    imputed_pseudoclusters = as.numeric(imputed_pseudoclusters),
    CI_width = as.numeric(CI_width)
  )

# Plotting subset only (do NOT use for stats)
min_n_plot <- 10

keep_top <- df_stats %>%
  dplyr::group_by(top_projection) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= min_n_plot) %>%
  dplyr::pull(top_projection)

keep_tgt <- df_stats %>%
  dplyr::group_by(proj_target) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n >= min_n_plot) %>%
  dplyr::pull(proj_target)

df_plot <- df_stats %>%
  filter(top_projection %in% keep_top,
         proj_target    %in% keep_tgt) %>%
  droplevels()

## ============================================================
## 04) Ridge plots (coarse groups colored; NMF grey)
## ============================================================

# Plasma palette (Python-equivalent) for coarse groups
top_projection_order <- c("OLF","Isocortex","HPF","CTXsp","CNU","TH","HY","MB","CB","P","MY","SP")

# If viridisLite is available, use it; otherwise fall back to base rainbow
if (requireNamespace("viridisLite", quietly = TRUE)) {
  topProj_colors_hex <- setNames(viridisLite::plasma(n = length(top_projection_order)),
                                 top_projection_order)
} else {
  topProj_colors_hex <- setNames(grDevices::rainbow(length(top_projection_order)),
                                 top_projection_order)
}

# enforce canonical order (even if some missing in df_plot)
df_plot <- df_plot %>%
  mutate(top_projection = factor(top_projection, levels = top_projection_order))

bw_numeric <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NA_real_)
  stats::density(x)$bw
}

order_by_median_str <- function(data, group_col, x_col) {
  # returns group levels ordered by median of x_col
  med_df <- data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(med = stats::median(.data[[x_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(med)
  as.character(med_df[[group_col]])
}

ridge_plot_str <- function(data, group_col, x_col, xlab, ylab,
                           alpha = 0.75, palette = NULL, levels = NULL) {
  
  data <- tibble::as_tibble(data) %>%
    dplyr::mutate(.x = as.numeric(.data[[x_col]])) %>%
    dplyr::filter(is.finite(.x)) %>%
    dplyr::mutate(.grp_raw = .data[[group_col]])
  
  bw <- bw_numeric(data$.x)
  if (!is.finite(bw)) stop("Bandwidth could not be estimated (too few finite x values).")
  
  levs <- if (!is.null(levels)) levels else order_by_median_str(data, ".grp_raw", ".x")
  data <- data %>% dplyr::mutate(.grp = factor(.grp_raw, levels = levs))
  
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .x, y = .grp))
  
  if (!is.null(palette)) {
    p <- ggplot2::ggplot(data, ggplot2::aes(x = .x, y = .grp, fill = .grp)) +
      ggplot2::scale_fill_manual(values = palette, guide = "none")
  }
  
  p +
    ggridges::geom_density_ridges(
      bandwidth = bw,
      scale = 1.2,
      rel_min_height = 0.01,
      alpha = alpha,
      color = "grey30",
      linewidth = 0.25
    ) +
    ggplot2::stat_summary(fun = stats::median, geom = "point", size = 1.5) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.y = ggplot2::element_blank()
    )
}
p1 <- ridge_plot_str(df_plot, "top_projection", "imputed_pseudoclusters",
                     "Imputed pseudocluster", "Top projection (coarse group)",
                     palette = topProj_colors_hex,
                     levels  = top_projection_order)

p2 <- ridge_plot_str(df_plot, "top_projection", "CI_width",
                     "CI width (uncertainty)", "Top projection (coarse group)",
                     palette = topProj_colors_hex,
                     levels  = top_projection_order)

p3 <- ridge_plot_str(df_plot, "proj_target", "imputed_pseudoclusters",
                     "Imputed pseudocluster", "NMF proj_target")

p4 <- ridge_plot_str(df_plot, "proj_target", "CI_width",
                     "CI width (uncertainty)", "NMF proj_target")

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)

## ============================================================
## 05) Stats (base R): use ALL groups (no dropping)
##     - omnibus ANOVA + Kruskal
##     - eta^2 from ANOVA SS
##     - TukeyHSD for ANOVA (only if >=2 groups)
##     - pairwise Wilcoxon with BH adjustment
## ============================================================

run_group_stats_base <- function(df, y, g, label = NULL) {
  if (is.null(label)) label <- paste(y, "~", g)
  stopifnot(all(c(y, g) %in% names(df)))
  
  yy <- as.numeric(df[[y]])
  gg <- as.factor(df[[g]])
  
  keep <- is.finite(yy) & !is.na(gg)
  yy <- yy[keep]
  gg <- droplevels(gg[keep])
  
  if (length(yy) == 0) stop("No finite values for: ", label)
  if (nlevels(gg) < 2) stop("Grouping var has <2 levels for: ", label)
  
  # group summary
  split_y <- split(yy, gg, drop = TRUE)
  n_by <- vapply(split_y, length, integer(1))
  mean_by <- vapply(split_y, mean, numeric(1), na.rm = TRUE)
  med_by  <- vapply(split_y, median, numeric(1), na.rm = TRUE)
  sd_by   <- vapply(split_y, sd, numeric(1), na.rm = TRUE)
  
  group_summary <- data.frame(
    group = names(n_by),
    n = as.integer(n_by),
    mean = as.numeric(mean_by),
    median = as.numeric(med_by),
    sd = as.numeric(sd_by),
    row.names = NULL
  )
  
  # ANOVA
  aov_fit <- aov(yy ~ gg)
  aov_tab <- summary(aov_fit)[[1]]
  aov_p <- aov_tab[["Pr(>F)"]][1]
  
  ss <- aov_tab[["Sum Sq"]]
  eta2 <- ss[1] / sum(ss)
  
  # Kruskal
  kw_p <- kruskal.test(yy ~ gg)$p.value
  
  omnibus <- data.frame(
    test = c("ANOVA", "KruskalWallis"),
    p_value = c(as.numeric(aov_p), as.numeric(kw_p)),
    eta2 = c(as.numeric(eta2), NA_real_),
    label = label,
    n_total = length(yy),
    n_groups = nlevels(gg),
    row.names = NULL
  )
  
  # Tukey (ANOVA post-hoc)
  tukey_df <- NULL
  if (nlevels(gg) >= 2) {
    tk <- TukeyHSD(aov_fit)
    tk_mat <- tk[[1]]
    tukey_df <- data.frame(
      contrast = rownames(tk_mat),
      diff = tk_mat[, "diff"],
      lwr  = tk_mat[, "lwr"],
      upr  = tk_mat[, "upr"],
      p_adj = tk_mat[, "p adj"],
      row.names = NULL
    )
  }
  
  # Pairwise Wilcoxon BH
  pw <- pairwise.wilcox.test(yy, gg, p.adjust.method = "BH")
  pw_mat <- pw$p.value  # matrix (upper triangle)
  
  list(
    group_summary = group_summary,
    omnibus = omnibus,
    tukey = tukey_df,
    pairwise_wilcox_p = pw_mat
  )
}

# Run the 4 conceptual tests (raw CI only; no log)
res_top_pc <- run_group_stats_base(df_stats, "imputed_pseudoclusters", "top_projection",
                                   label = "pseudocluster ~ top_projection (all groups)")
res_top_ci <- run_group_stats_base(df_stats, "CI_width", "top_projection",
                                   label = "CI_width ~ top_projection (all groups)")
res_tgt_pc <- run_group_stats_base(df_stats, "imputed_pseudoclusters", "proj_target",
                                   label = "pseudocluster ~ proj_target (all groups)")
res_tgt_ci <- run_group_stats_base(df_stats, "CI_width", "proj_target",
                                   label = "CI_width ~ proj_target (all groups)")

omnibus_summary <- rbind(res_top_pc$omnibus, res_top_ci$omnibus,
                         res_tgt_pc$omnibus, res_tgt_ci$omnibus)
print(omnibus_summary)

# Save stats outputs
write.csv(omnibus_summary, "stats_omnibus_summary_baseR_allGroups_rawCI.csv", row.names = FALSE)

write.csv(res_top_pc$group_summary, "group_summary_pseudocluster_by_top_projection_allGroups.csv", row.names = FALSE)
write.csv(res_tgt_pc$group_summary, "group_summary_pseudocluster_by_proj_target_allGroups.csv", row.names = FALSE)
write.csv(res_top_ci$group_summary, "group_summary_CIwidth_by_top_projection_allGroups.csv", row.names = FALSE)
write.csv(res_tgt_ci$group_summary, "group_summary_CIwidth_by_proj_target_allGroups.csv", row.names = FALSE)

if (!is.null(res_top_pc$tukey)) write.csv(res_top_pc$tukey, "tukey_pseudocluster_by_top_projection_allGroups.csv", row.names = FALSE)
if (!is.null(res_tgt_pc$tukey)) write.csv(res_tgt_pc$tukey, "tukey_pseudocluster_by_proj_target_allGroups.csv", row.names = FALSE)
if (!is.null(res_top_ci$tukey)) write.csv(res_top_ci$tukey, "tukey_CIwidth_by_top_projection_allGroups.csv", row.names = FALSE)
if (!is.null(res_tgt_ci$tukey)) write.csv(res_tgt_ci$tukey, "tukey_CIwidth_by_proj_target_allGroups.csv", row.names = FALSE)

write.csv(as.data.frame(res_top_pc$pairwise_wilcox_p), "wilcoxBH_pseudocluster_by_top_projection_allGroups.csv", row.names = TRUE)
write.csv(as.data.frame(res_tgt_pc$pairwise_wilcox_p), "wilcoxBH_pseudocluster_by_proj_target_allGroups.csv", row.names = TRUE)
write.csv(as.data.frame(res_top_ci$pairwise_wilcox_p), "wilcoxBH_CIwidth_by_top_projection_allGroups.csv", row.names = TRUE)
write.csv(as.data.frame(res_tgt_ci$pairwise_wilcox_p), "wilcoxBH_CIwidth_by_proj_target_allGroups.csv", row.names = TRUE)

## ============================================================
## 06) Regression: pseudocluster is the response (0–1 continuous)
##     CI_width is uncertainty/quality (NOT the outcome)
##     Two complementary analyses:
##       A) Unweighted: adjust for CI_width as a covariate
##       B) Weighted: downweight uncertain cells using weights = 1/(CI_width^2)
## ============================================================

df_reg <- df_stats %>%
  filter(is.finite(imputed_pseudoclusters),
         is.finite(CI_width),
         !is.na(top_projection),
         !is.na(proj_target)) %>%
  mutate(
    top_projection = droplevels(as.factor(top_projection)),
    proj_target    = droplevels(as.factor(proj_target))
  )

cat("Regression rows:", nrow(df_reg), "\n")
cat("Top projection levels:", nlevels(df_reg$top_projection), "\n")
cat("Proj_target levels:", nlevels(df_reg$proj_target), "\n")

## --- A) Covariate-adjusted models (unweighted) ---
m_top_adj <- lm(imputed_pseudoclusters ~ top_projection + CI_width, data = df_reg)
m_tgt_adj <- lm(imputed_pseudoclusters ~ proj_target + CI_width, data = df_reg)

cat("\n[Unweighted] pseudocluster ~ top_projection + CI_width\n")
print(summary(m_top_adj))
cat("\nANOVA (type I) for top_projection term:\n")
print(anova(m_top_adj))

cat("\n[Unweighted] pseudocluster ~ proj_target + CI_width\n")
print(summary(m_tgt_adj))
cat("\nANOVA (type I) for proj_target term:\n")
print(anova(m_tgt_adj))

## --- B) Weighted models (more weight to precise pseudoclusters) ---
## CI_width is in (0,1]; smaller = more confident. A standard choice:
## weights = 1 / (CI_width^2 + eps)
eps <- 1e-8
w <- 1 / (df_reg$CI_width^2 + eps)

m_top_w <- lm(imputed_pseudoclusters ~ top_projection, data = df_reg, weights = w)
m_tgt_w <- lm(imputed_pseudoclusters ~ proj_target, data = df_reg, weights = w)

cat("\n[Weighted] pseudocluster ~ top_projection  (weights = 1/CI_width^2)\n")
print(summary(m_top_w))
cat("\nWeighted ANOVA-like table:\n")
print(anova(m_top_w))

cat("\n[Weighted] pseudocluster ~ proj_target  (weights = 1/CI_width^2)\n")
print(summary(m_tgt_w))
cat("\nWeighted ANOVA-like table:\n")
print(anova(m_tgt_w))

## --- Effect size (base R): partial eta^2 for the categorical term ---
partial_eta2 <- function(lm_fit, term_name) {
  a <- anova(lm_fit)
  if (!(term_name %in% rownames(a))) return(NA_real_)
  ss_term <- a[term_name, "Sum Sq"]
  ss_res  <- a["Residuals", "Sum Sq"]
  ss_term / (ss_term + ss_res)
}

eta_top_adj <- partial_eta2(m_top_adj, "top_projection")
eta_tgt_adj <- partial_eta2(m_tgt_adj, "proj_target")
eta_top_w   <- partial_eta2(m_top_w, "top_projection")
eta_tgt_w   <- partial_eta2(m_tgt_w, "proj_target")

cat("\nPartial eta^2 (unweighted, adjusted):\n")
cat("  top_projection:", eta_top_adj, "\n")
cat("  proj_target   :", eta_tgt_adj, "\n")

cat("\nPartial eta^2 (weighted):\n")
cat("  top_projection:", eta_top_w, "\n")
cat("  proj_target   :", eta_tgt_w, "\n")

## --- Optional: robust check using Kruskal-Wallis (ignores CI) ---
cat("\nKruskal-Wallis (no CI adjustment; distribution-free):\n")
cat("  pseudocluster ~ top_projection p =",
    kruskal.test(imputed_pseudoclusters ~ top_projection, data = df_reg)$p.value, "\n")
cat("  pseudocluster ~ proj_target p =",
    kruskal.test(imputed_pseudoclusters ~ proj_target, data = df_reg)$p.value, "\n")

## Save a compact regression summary
reg_summary <- data.frame(
  model = c("top_adj", "tgt_adj", "top_weighted", "tgt_weighted"),
  n = c(nobs(m_top_adj), nobs(m_tgt_adj), nobs(m_top_w), nobs(m_tgt_w)),
  r2 = c(summary(m_top_adj)$r.squared,
         summary(m_tgt_adj)$r.squared,
         summary(m_top_w)$r.squared,
         summary(m_tgt_w)$r.squared),
  partial_eta2 = c(eta_top_adj, eta_tgt_adj, eta_top_w, eta_tgt_w),
  stringsAsFactors = FALSE
)
write.csv(reg_summary, "regression_summary_pseudocluster_response.csv", row.names = FALSE)

summary(w)
quantile(w, c(0.5, 0.9, 0.99, 0.999))
sum(w > quantile(w, 0.99))
plot(df_reg$CI_width, df_reg$imputed_pseudoclusters, pch=16, cex=0.6,
     xlab="CI_width (uncertainty)", ylab="imputed_pseudoclusters")
abline(lm(imputed_pseudoclusters ~ CI_width, data=df_reg), lwd=2)
cor(df_reg$CI_width, df_reg$imputed_pseudoclusters, method="spearman")
