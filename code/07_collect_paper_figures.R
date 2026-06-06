# Reorganize /results/ for publication.
#
# All analysis stages write where they always have (the three per-cohort
# folders under /results/). This final step does no analysis — it only
# relocates outputs so the published structure matches AIND publication
# standards (https://docs.allenneuraldynamics.org/en/latest/policies_practices/publication_standards.html):
#
#   results/
#     paper_figures/FigureS5/   canonically-named copies of the manuscript panels
#     other_figures/            everything else (per-cohort QC/exploratory output)
#
# Panel -> source-file mapping confirmed by the paper authors on issue #5
# (only panels e, f, g of Figure S5 come from this capsule; a/b are hand-drawn
# schematics, and c/d/h are not produced here).

RESULTS_DIR <- "/results"
COMBINED_DIR <- file.path(RESULTS_DIR, "BARseq_780345-780346_combined")
PAPER_FIG_DIR <- file.path(RESULTS_DIR, "paper_figures", "FigureS5")
OTHER_FIG_DIR <- file.path(RESULTS_DIR, "other_figures")

# source basename (in COMBINED_DIR) -> canonical published name
panel_map <- list(
  list(src = "Combined_ipsi-contra_projections_heatmap_top_region_sorted.pdf",
       dst = "FigureS5e_ipsi_contra_projection_heatmap.pdf"),
  list(src = "sorted_proj_heatmap_ipsi-contra.pdf",
       dst = "FigureS5f_sorted_projection_heatmap_ipsi_contra.pdf"),
  list(src = "rank_corr_ipsi-contra.pdf",
       dst = "FigureS5g_rank_correlation_ipsi_contra.pdf")
)

# --- 1. Copy the manuscript panels into paper_figures/FigureS5/ ---------------
dir.create(PAPER_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
for (panel in panel_map) {
  src <- file.path(COMBINED_DIR, panel$src)
  dst <- file.path(PAPER_FIG_DIR, panel$dst)
  if (!file.exists(src)) {
    stop(sprintf("Expected paper-figure source not found: %s", src))
  }
  if (!file.copy(src, dst, overwrite = TRUE)) {
    stop(sprintf("Failed to copy %s -> %s", src, dst))
  }
  cat(sprintf("paper_figures: %s -> %s\n", panel$src, panel$dst))
}

# --- 2. Move the per-cohort folders under other_figures/ ----------------------
dir.create(OTHER_FIG_DIR, recursive = TRUE, showWarnings = FALSE)
cohort_dirs <- c("BARseq_780345", "BARseq_780346", "BARseq_780345-780346_combined")
for (cohort in cohort_dirs) {
  from <- file.path(RESULTS_DIR, cohort)
  to   <- file.path(OTHER_FIG_DIR, cohort)
  if (!dir.exists(from)) {
    cat(sprintf("other_figures: %s not present, skipping\n", cohort))
    next
  }
  if (!file.rename(from, to)) {
    stop(sprintf("Failed to move %s -> %s", from, to))
  }
  cat(sprintf("other_figures: moved %s/\n", cohort))
}

cat("Paper-figure collection complete.\n")
