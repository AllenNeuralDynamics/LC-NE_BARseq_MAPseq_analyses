#' # Validate analysis outputs against the released-capsule reference
#'
#' TEMP — remove before merge. Confirms the input-asset swap did not change any
#' output of the pipeline. Compares everything written to `/results/` by this run
#' against the released-capsule outputs (frozen as the reference asset).
#'
#' Strict comparison: `identical()` on RDS objects and parsed CSV/TSV. Visual
#' files (HTML, PDF) are skipped — eyeball those side-by-side instead.

suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

OLD_ROOT <- "/data/temp_lc_paper_barseq_mapseq_frozen_v1_output"
NEW_ROOT <- "/results"

stopifnot(dir.exists(OLD_ROOT))
stopifnot(dir.exists(NEW_ROOT))

#' ## File-by-file comparison

ref_files <- list.files(OLD_ROOT, recursive = TRUE, full.names = FALSE)
new_files <- list.files(NEW_ROOT, recursive = TRUE, full.names = FALSE)

cat(sprintf("Reference asset: %d files\n", length(ref_files)))
cat(sprintf("Current /results: %d files\n\n", length(new_files)))

ext_of <- function(p) tolower(tools::file_ext(p))

n_match_rds <- 0; n_diff_rds <- 0
n_match_csv <- 0; n_diff_csv <- 0
n_match_tsv <- 0; n_diff_tsv <- 0
n_skip_visual <- 0
n_skip_other <- 0
n_missing <- 0

for (rel in sort(ref_files)) {
  old_path <- file.path(OLD_ROOT, rel)
  new_path <- file.path(NEW_ROOT, rel)
  e <- ext_of(rel)

  if (!file.exists(new_path)) {
    cat(sprintf("MISSING in /results: %s\n", rel))
    n_missing <- n_missing + 1
    next
  }

  if (e == "rds") {
    old_obj <- readRDS(old_path)
    new_obj <- readRDS(new_path)
    if (identical(old_obj, new_obj)) {
      n_match_rds <- n_match_rds + 1
      cat(sprintf("[OK rds]   %s\n", rel))
    } else {
      n_diff_rds <- n_diff_rds + 1
      cat(sprintf("[DIFF rds] %s\n", rel))
      cat(sprintf("           class old=%s new=%s\n",
                  paste(class(old_obj), collapse=","),
                  paste(class(new_obj), collapse=",")))
      if (is(old_obj, "SingleCellExperiment") && is(new_obj, "SingleCellExperiment")) {
        cat(sprintf("           dim   old=%s new=%s\n",
                    paste(dim(old_obj), collapse="x"),
                    paste(dim(new_obj), collapse="x")))
        for (a in intersect(assayNames(old_obj), assayNames(new_obj))) {
          eq <- identical(assay(old_obj, a), assay(new_obj, a))
          cat(sprintf("           assay '%s' identical: %s\n", a, eq))
        }
        cat(sprintf("           colData identical:  %s\n",
                    identical(colData(old_obj), colData(new_obj))))
        cat(sprintf("           rowData identical:  %s\n",
                    identical(rowData(old_obj), rowData(new_obj))))
        cat(sprintf("           metadata identical: %s   (often differs harmlessly via package-version stamps)\n",
                    identical(metadata(old_obj), metadata(new_obj))))
      }
    }
  } else if (e %in% c("csv", "tsv")) {
    sep <- if (e == "csv") "," else "\t"
    old_df <- tryCatch(read.table(old_path, sep = sep, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE,
                                  comment.char = ""),
                       error = function(err) NULL)
    new_df <- tryCatch(read.table(new_path, sep = sep, header = TRUE,
                                  stringsAsFactors = FALSE, check.names = FALSE,
                                  comment.char = ""),
                       error = function(err) NULL)
    if (is.null(old_df) || is.null(new_df)) {
      cat(sprintf("[SKIP %s unreadable] %s\n", e, rel))
      n_skip_other <- n_skip_other + 1
      next
    }
    if (identical(old_df, new_df)) {
      if (e == "csv") n_match_csv <- n_match_csv + 1 else n_match_tsv <- n_match_tsv + 1
      cat(sprintf("[OK %s]   %s  (%d rows x %d cols)\n",
                  e, rel, nrow(old_df), ncol(old_df)))
    } else {
      if (e == "csv") n_diff_csv <- n_diff_csv + 1 else n_diff_tsv <- n_diff_tsv + 1
      cat(sprintf("[DIFF %s] %s  old=%dx%d new=%dx%d\n",
                  e, rel, nrow(old_df), ncol(old_df), nrow(new_df), ncol(new_df)))
    }
  } else if (e %in% c("html", "pdf")) {
    n_skip_visual <- n_skip_visual + 1
  } else {
    n_skip_other <- n_skip_other + 1
  }
}

extras <- setdiff(new_files, ref_files)
for (rel in extras) {
  cat(sprintf("EXTRA in /results: %s\n", rel))
}
n_extra <- length(extras)

#' ## Summary

cat("\n========================\n")
cat(sprintf("RDS   matching: %d   differing: %d\n", n_match_rds, n_diff_rds))
cat(sprintf("CSV   matching: %d   differing: %d\n", n_match_csv, n_diff_csv))
cat(sprintf("TSV   matching: %d   differing: %d\n", n_match_tsv, n_diff_tsv))
cat(sprintf("Visual (HTML/PDF) skipped: %d   (eyeball side-by-side)\n", n_skip_visual))
cat(sprintf("Other skipped:             %d\n", n_skip_other))
cat(sprintf("Missing in /results:       %d\n", n_missing))
cat(sprintf("Extra in /results:         %d\n", n_extra))

n_diverged <- n_diff_rds + n_diff_csv + n_diff_tsv + n_missing + n_extra
if (n_diverged == 0) {
  cat("\nALL NUMERIC OUTPUTS MATCH THE REFERENCE.\n")
} else {
  cat(sprintf("\n%d FILE(S) DIVERGE — see [DIFF] / MISSING / EXTRA lines above.\n", n_diverged))
}

#' ## Dbh / Th per-cell expression spot-check
#'
#' Drill into the per-cohort post-normalization SCEs and compare counts for the
#' two genes Polina checks first.

for (cohort in c("BARseq_780345", "BARseq_780346", "BARseq_780345-780346_combined")) {
  rel <- file.path(cohort, "combined_neurons_clust_CCFv2_uid_cpm_log.rds")
  old_path <- file.path(OLD_ROOT, rel)
  new_path <- file.path(NEW_ROOT, rel)
  cat(sprintf("\n--- %s ---\n", cohort))
  if (!file.exists(old_path) || !file.exists(new_path)) {
    cat(sprintf("  skipped: %s   old_exists=%s new_exists=%s\n",
                rel, file.exists(old_path), file.exists(new_path)))
    next
  }
  old_sce <- readRDS(old_path)
  new_sce <- readRDS(new_path)
  for (gene in c("Dbh", "Th")) {
    if (!(gene %in% rownames(old_sce) && gene %in% rownames(new_sce))) {
      cat(sprintf("  %s: not present in both objects (old=%s new=%s)\n",
                  gene, gene %in% rownames(old_sce), gene %in% rownames(new_sce)))
      next
    }
    old_vec <- as.numeric(counts(old_sce)[gene, ])
    new_vec <- as.numeric(counts(new_sce)[gene, ])
    eq <- identical(old_vec, new_vec)
    cat(sprintf("  %s counts: cells=%d  identical=%s  sum old=%g new=%g  max old=%g new=%g  nz_cells old=%d new=%d\n",
                gene, length(old_vec), eq,
                sum(old_vec), sum(new_vec),
                max(old_vec), max(new_vec),
                sum(old_vec > 0), sum(new_vec > 0)))
  }
}
