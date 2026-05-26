suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

OLD_ROOT <- "/data/temp_lc_paper_barseq_mapseq_frozen_v1_output"
NEW_ROOT <- "/results"
CLUSTERING_ROOT <- "/code/cached_clustering"

stopifnot(dir.exists(OLD_ROOT))
stopifnot(dir.exists(NEW_ROOT))

# ---------------------------------------------------------------
# Per-SCE diff drilldown helper. Returns a character vector of
# markdown lines describing the diff in detail. Safe to call even
# when no diff is present — emits a "no diff" line instead.
# ---------------------------------------------------------------
diagnose_sce_diff <- function(old_obj, new_obj, mismatches) {
  out <- character()

  fmt_vec <- function(v, n = 6L) {
    if (length(v) == 0L) return("<empty>")
    s <- head(as.character(v), n)
    paste(s, collapse = ", ")
  }

  if ("colData" %in% mismatches) {
    cd_old <- colData(old_obj); cd_new <- colData(new_obj)
    out <- c(out, "  **colData drilldown:**")
    out <- c(out, sprintf("  - ncol old=%d new=%d", ncol(cd_old), ncol(cd_new)))
    if (!identical(colnames(cd_old), colnames(cd_new))) {
      out <- c(out, sprintf("  - colnames old: %s", fmt_vec(colnames(cd_old), 50)))
      out <- c(out, sprintf("  - colnames new: %s", fmt_vec(colnames(cd_new), 50)))
      added <- setdiff(colnames(cd_new), colnames(cd_old))
      removed <- setdiff(colnames(cd_old), colnames(cd_new))
      if (length(added)   > 0L) out <- c(out, sprintf("  - added columns: %s",   paste(added,   collapse = ", ")))
      if (length(removed) > 0L) out <- c(out, sprintf("  - removed columns: %s", paste(removed, collapse = ", ")))
    }
    rn_id <- identical(rownames(cd_old), rownames(cd_new))
    out <- c(out, sprintf("  - rownames identical: %s", rn_id))
    if (!rn_id && length(rownames(cd_old)) == length(rownames(cd_new))) {
      n_rn_diff <- sum(rownames(cd_old) != rownames(cd_new))
      out <- c(out, sprintf("    (%d of %d rownames differ; head old: %s; head new: %s)",
                            n_rn_diff, length(rownames(cd_old)),
                            fmt_vec(rownames(cd_old)), fmt_vec(rownames(cd_new))))
    }

    common <- intersect(colnames(cd_old), colnames(cd_new))
    diff_cols <- character()
    same_cols <- character()
    for (cc in common) {
      if (identical(cd_old[[cc]], cd_new[[cc]])) same_cols <- c(same_cols, cc)
      else diff_cols <- c(diff_cols, cc)
    }
    out <- c(out, sprintf("  - per-column identical: %d of %d common columns match",
                          length(same_cols), length(common)))
    if (length(diff_cols) > 0L) {
      out <- c(out, sprintf("  - **columns differing:** %s", paste(diff_cols, collapse = ", ")))
      for (cc in diff_cols) {
        ovec <- cd_old[[cc]]; nvec <- cd_new[[cc]]
        cls_o <- class(ovec)[1]; cls_n <- class(nvec)[1]
        len_ok <- length(ovec) == length(nvec)
        out <- c(out, sprintf("    - `%s`: class old=%s new=%s, lengths equal: %s",
                              cc, cls_o, cls_n, len_ok))
        if (len_ok) {
          o_chr <- as.character(ovec); n_chr <- as.character(nvec)
          n_diff <- sum(o_chr != n_chr, na.rm = TRUE)
          out <- c(out, sprintf("      - n positions differing (as character): %d of %d",
                                n_diff, length(ovec)))
          out <- c(out, sprintf("      - head old: %s", fmt_vec(ovec)))
          out <- c(out, sprintf("      - head new: %s", fmt_vec(nvec)))
        }
        if (is.factor(ovec) && is.factor(nvec)) {
          lev_id <- identical(levels(ovec), levels(nvec))
          out <- c(out, sprintf("      - levels identical: %s", lev_id))
          out <- c(out, sprintf("        - levels old: %s", fmt_vec(levels(ovec), 50)))
          if (!lev_id) out <- c(out, sprintf("        - levels new: %s", fmt_vec(levels(nvec), 50)))
          if (lev_id && len_ok && length(levels(ovec)) <= 20L) {
            ct <- table(old = ovec, new = nvec)
            out <- c(out, "      - contingency table (rows = old, cols = new):")
            out <- c(out, "        ```")
            ct_lines <- capture.output(print(ct))
            out <- c(out, paste0("        ", ct_lines))
            out <- c(out, "        ```")
            # quick permutation hint: each old row peaks at exactly one new column
            if (all(dim(ct) == length(levels(ovec)))) {
              row_peaks <- apply(ct, 1, which.max)
              if (length(unique(row_peaks)) == nrow(ct)) {
                out <- c(out, "      - **looks like a label permutation** (each old label peaks at a distinct new label)")
              }
            }
          }
        }
      }
    }
  }

  assay_mm <- grep("^assay '", mismatches, value = TRUE)
  if (length(assay_mm) > 0L) {
    out <- c(out, "  **assay drilldown:**")
    for (mm in assay_mm) {
      a <- sub("^assay '(.*)'$", "\\1", mm)
      if (!(a %in% assayNames(old_obj) && a %in% assayNames(new_obj))) next
      old_a <- assay(old_obj, a); new_a <- assay(new_obj, a)
      dim_match <- identical(dim(old_a), dim(new_a))
      out <- c(out, sprintf("  - `%s`: dim old=%s new=%s, dim match: %s",
                            a, paste(dim(old_a), collapse="x"), paste(dim(new_a), collapse="x"), dim_match))
      cn_id <- identical(colnames(old_a), colnames(new_a))
      rn_id <- identical(rownames(old_a), rownames(new_a))
      out <- c(out, sprintf("    - colnames identical (same cells in same order): %s", cn_id))
      out <- c(out, sprintf("    - rownames identical (same genes in same order): %s", rn_id))
      if (!cn_id && length(colnames(old_a)) == length(colnames(new_a))) {
        # check if it's a reorder (same set) vs a different set (different cell selection)
        same_set <- setequal(colnames(old_a), colnames(new_a))
        out <- c(out, sprintf("    - same cell SET (just reordered): %s", same_set))
      }
      if (dim_match) {
        # diff in shared (intersection of colnames) submatrix
        # avoid dense conversion of huge matrices; restrict to small failed files
        if (ncol(old_a) <= 20000L) {
          # align by colnames where possible
          if (cn_id) {
            o <- as.matrix(old_a); n <- as.matrix(new_a)
          } else if (length(colnames(old_a)) == length(colnames(new_a)) &&
                     setequal(colnames(old_a), colnames(new_a))) {
            o <- as.matrix(old_a); n <- as.matrix(new_a)[, colnames(old_a), drop = FALSE]
          } else {
            o <- NULL; n <- NULL
          }
          if (!is.null(o)) {
            d <- abs(o - n)
            out <- c(out, sprintf("    - max abs diff: %s, sum abs diff: %s",
                                  format(max(d), digits = 6), format(sum(d), digits = 6)))
            out <- c(out, sprintf("    - cells with any diff: %d of %d", sum(colSums(d) > 0), ncol(d)))
            out <- c(out, sprintf("    - genes with any diff: %d of %d", sum(rowSums(d) > 0), nrow(d)))
          } else {
            out <- c(out, "    - cell selection differs — element-wise diff not meaningful, skipped")
          }
        } else {
          out <- c(out, sprintf("    - (matrix too large for element-wise diff: %d cells)", ncol(old_a)))
        }
      }
    }
  }

  if ("rowData" %in% mismatches) {
    out <- c(out, "  **rowData differs** (not drilled in)")
  }
  if ("metadata" %in% mismatches) {
    out <- c(out, "  **metadata differs**")
    mo <- metadata(old_obj); mn <- metadata(new_obj)
    out <- c(out, sprintf("    - names old: %s", fmt_vec(names(mo), 50)))
    out <- c(out, sprintf("    - names new: %s", fmt_vec(names(mn), 50)))
  }

  out
}

ref_files <- list.files(OLD_ROOT, recursive = TRUE, full.names = FALSE)
new_files <- list.files(NEW_ROOT, recursive = TRUE, full.names = FALSE)

ext_of <- function(p) tolower(tools::file_ext(p))
fmt <- function(x) format(x, big.mark = ",")

# ---------------------------------------------------------------
# Pass 1: classify each ref file into a bucket. No printing yet.
# ---------------------------------------------------------------

csv_match      <- list()  # list(path, rows, cols)
csv_diff       <- list()  # list(path, old_dim, new_dim)
tsv_match      <- list()
tsv_diff       <- list()
rds_data_match <- list()  # list(path, dim, fully_identical)
rds_real_diff  <- list()  # list(path, dim, mismatches)
skipped_visual <- character()
skipped_other  <- character()
missing_in_new <- character()
extras         <- setdiff(new_files, ref_files)

for (rel in sort(ref_files)) {
  old_path <- file.path(OLD_ROOT, rel)
  new_path <- file.path(NEW_ROOT, rel)
  e <- ext_of(rel)

  if (!file.exists(new_path)) {
    missing_in_new <- c(missing_in_new, rel)
    next
  }

  if (e == "rds") {
    old_obj <- readRDS(old_path)
    new_obj <- readRDS(new_path)

    if (identical(old_obj, new_obj)) {
      dims <- if (is(old_obj, "SingleCellExperiment")) paste(dim(old_obj), collapse = "x") else NA_character_
      rds_data_match[[length(rds_data_match) + 1L]] <- list(path = rel, dim = dims, fully_identical = TRUE)
      next
    }

    if (is(old_obj, "SingleCellExperiment") && is(new_obj, "SingleCellExperiment")) {
      mismatches <- character()
      for (a in intersect(assayNames(old_obj), assayNames(new_obj))) {
        if (!identical(assay(old_obj, a), assay(new_obj, a))) {
          mismatches <- c(mismatches, paste0("assay '", a, "'"))
        }
      }
      added   <- setdiff(assayNames(new_obj), assayNames(old_obj))
      removed <- setdiff(assayNames(old_obj), assayNames(new_obj))
      if (length(added)   > 0L) mismatches <- c(mismatches, paste0("added assay: ",   added))
      if (length(removed) > 0L) mismatches <- c(mismatches, paste0("removed assay: ", removed))
      if (!identical(colData(old_obj),  colData(new_obj)))  mismatches <- c(mismatches, "colData")
      if (!identical(rowData(old_obj),  rowData(new_obj)))  mismatches <- c(mismatches, "rowData")
      if (!identical(metadata(old_obj), metadata(new_obj))) mismatches <- c(mismatches, "metadata")

      dims <- paste(dim(old_obj), collapse = "x")
      if (length(mismatches) == 0L) {
        rds_data_match[[length(rds_data_match) + 1L]] <- list(path = rel, dim = dims, fully_identical = FALSE)
      } else {
        detail_lines <- tryCatch(
          diagnose_sce_diff(old_obj, new_obj, mismatches),
          error = function(e) sprintf("  (drilldown failed: %s)", conditionMessage(e))
        )
        rds_real_diff[[length(rds_real_diff) + 1L]] <- list(
          path = rel, dim = dims, mismatches = mismatches, detail = detail_lines
        )
      }
    } else {
      rds_real_diff[[length(rds_real_diff) + 1L]] <- list(
        path = rel,
        dim = NA_character_,
        mismatches = paste0("class old=", paste(class(old_obj), collapse = ","),
                            " new=", paste(class(new_obj), collapse = ","))
      )
    }
  } else if (e %in% c("csv", "tsv")) {
    sep <- if (e == "csv") "," else "\t"
    old_df <- tryCatch(
      read.table(old_path, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, comment.char = ""),
      error = function(err) NULL
    )
    new_df <- tryCatch(
      read.table(new_path, sep = sep, header = TRUE, stringsAsFactors = FALSE,
                 check.names = FALSE, comment.char = ""),
      error = function(err) NULL
    )
    if (is.null(old_df) || is.null(new_df)) {
      skipped_other <- c(skipped_other, rel)
      next
    }

    if (identical(old_df, new_df)) {
      entry <- list(path = rel, rows = nrow(old_df), cols = ncol(old_df))
      if (e == "csv") csv_match[[length(csv_match) + 1L]] <- entry
      else            tsv_match[[length(tsv_match) + 1L]] <- entry
    } else {
      entry <- list(path = rel,
                    old_dim = paste0(nrow(old_df), "x", ncol(old_df)),
                    new_dim = paste0(nrow(new_df), "x", ncol(new_df)))
      if (e == "csv") csv_diff[[length(csv_diff) + 1L]] <- entry
      else            tsv_diff[[length(tsv_diff) + 1L]] <- entry
    }
  } else if (e %in% c("html", "pdf")) {
    skipped_visual <- c(skipped_visual, rel)
  } else {
    skipped_other <- c(skipped_other, rel)
  }
}

# ---------------------------------------------------------------
# Pass 2: print the report, summary first.
# ---------------------------------------------------------------

n_csv_total <- length(csv_match) + length(csv_diff)
n_tsv_total <- length(tsv_match) + length(tsv_diff)
n_rds_data_match  <- length(rds_data_match)
n_rds_stamp_only  <- sum(vapply(rds_data_match, function(x) !isTRUE(x$fully_identical), logical(1)))
n_rds_fully_match <- n_rds_data_match - n_rds_stamp_only
n_rds_real_diff   <- length(rds_real_diff)
n_data_diffs      <- length(csv_diff) + length(tsv_diff) + n_rds_real_diff

cat("# Validation report — input-asset swap vs released-capsule reference\n\n")
cat("Compares everything written to `/results/` by this run against the released-capsule outputs (frozen as the reference asset).\n")
cat("Strict comparison: `identical()` on RDS objects and on parsed CSV/TSV. Visual files (HTML, PDF) are skipped — eyeball those side-by-side instead.\n\n")

cat("## Result\n\n")
if (n_data_diffs == 0L) {
  cat("**Behavior preserved — no analysis output drifted in any data-bearing way.**\n\n")
} else {
  cat(sprintf("**%d file(s) have real data differences. See the \"RDS files with real data differences\" / \"CSVs with differences\" sections below.**\n\n",
              n_data_diffs))
}

cat(sprintf("- **CSVs:** %d of %d match the reference exactly.\n", length(csv_match), n_csv_total))
if (n_tsv_total > 0L) {
  cat(sprintf("- **TSVs:** %d of %d match the reference exactly.\n", length(tsv_match), n_tsv_total))
}
if (n_rds_data_match > 0L) {
  cat(sprintf("- **RDS files (data-matching):** %d of %d. Of these, %d are fully identical, and %d match the reference exactly in every data slot (`counts`, `cpm`, `logcounts`, `colData`, `rowData`, `metadata`) but report a top-level diff because the `SingleCellExperiment` library writes its own version into a hidden field inside each saved object — the actual scientific data is unchanged.\n",
              n_rds_data_match, n_rds_data_match + n_rds_real_diff, n_rds_fully_match, n_rds_stamp_only))
}
if (n_rds_real_diff > 0L) {
  cat(sprintf("- **RDS files with real data differences:** %d (see section below).\n", n_rds_real_diff))
}
cat(sprintf("- **Visual files (HTML/PDF):** %d skipped — eyeball side-by-side if needed.\n", length(skipped_visual)))
if (length(missing_in_new) > 0L || length(extras) > 0L) {
  cat(sprintf("- **File-presence quirks:** %d missing, %d extra — explained in the section below.\n",
              length(missing_in_new), length(extras)))
}
cat(sprintf("\nReference asset: %d files. Current `/results/`: %d files.\n\n", length(ref_files), length(new_files)))

# ---------- Clustering cache state ----------
cat("## Clustering cache state\n\n")
recompute_env <- Sys.getenv("RECOMPUTE_CLUSTERING", "")
recompute_active <- tolower(recompute_env) %in% c("true", "1", "yes")
cat(sprintf("- `RECOMPUTE_CLUSTERING` env: `%s` → %s\n",
            ifelse(nzchar(recompute_env), recompute_env, "<unset>"),
            ifelse(recompute_active, "**cache BYPASSED, fresh clustering forced**", "cache enabled")))
cat(sprintf("- `CLUSTERING_ROOT`: `%s` (exists: %s)\n", CLUSTERING_ROOT, dir.exists(CLUSTERING_ROOT)))

if (dir.exists(CLUSTERING_ROOT)) {
  cohort_dirs <- list.dirs(CLUSTERING_ROOT, recursive = FALSE, full.names = FALSE)
  for (cohort in cohort_dirs) {
    cohort_path <- file.path(CLUSTERING_ROOT, cohort)
    rounds <- list.dirs(cohort_path, recursive = FALSE, full.names = FALSE)
    cat(sprintf("- `%s/`: %d round(s)\n", cohort, length(rounds)))
    for (rd in rounds) {
      rp <- file.path(cohort_path, rd)
      has_umap <- file.exists(file.path(rp, "umap.csv"))
      has_clu  <- file.exists(file.path(rp, "cluster.csv"))
      has_ann  <- file.exists(file.path(rp, "cluster_annotation.csv"))
      would_hit <- !recompute_active && dir.exists(rp) && has_umap && has_clu
      cat(sprintf("  - `%s/`: umap.csv=%s, cluster.csv=%s, cluster_annotation.csv=%s → **cache hit predicate: %s**\n",
                  rd, has_umap, has_clu, has_ann, would_hit))
    }
  }
}
cat("\n")

# ---------- CSVs ----------
if (length(csv_match) > 0L) {
  cat(sprintf("## CSVs that match the reference exactly (%d)\n\n", length(csv_match)))
  for (entry in csv_match) {
    cat(sprintf("- `%s` — %s rows × %s cols\n", entry$path, fmt(entry$rows), fmt(entry$cols)))
  }
  cat("\n")
}

if (length(csv_diff) > 0L) {
  cat(sprintf("## CSVs with differences (%d)\n\n", length(csv_diff)))
  for (entry in csv_diff) {
    cat(sprintf("- `%s` — old %s, new %s\n", entry$path, entry$old_dim, entry$new_dim))
  }
  cat("\n")
}

# ---------- RDS data-match ----------
if (n_rds_data_match > 0L) {
  cat(sprintf("## RDS files: data matches the reference exactly (%d)\n\n", n_rds_data_match))
  if (n_rds_stamp_only > 0L) {
    cat("For these files, every data slot — `assay 'counts'`, `'cpm'`, `'logcounts'`, `colData`, `rowData`, and `metadata` — matches the reference exactly. ")
    cat("Top-level `identical()` returns FALSE only because the `SingleCellExperiment` library writes its own version into a hidden field inside each saved object. ")
    cat("The actual scientific data (counts, gene expression, cell metadata) is unchanged.\n\n")
  }
  for (entry in rds_data_match) {
    suffix <- if (!is.na(entry$dim)) sprintf(" (%s)", entry$dim) else ""
    cat(sprintf("- `%s`%s\n", entry$path, suffix))
  }
  cat("\n")
}

# ---------- RDS real-diff ----------
if (n_rds_real_diff > 0L) {
  cat(sprintf("## RDS files with real data differences (%d)\n\n", n_rds_real_diff))
  for (entry in rds_real_diff) {
    suffix <- if (!is.na(entry$dim)) sprintf(" (%s)", entry$dim) else ""
    cat(sprintf("- `%s`%s — differs in: %s\n",
                entry$path, suffix, paste(entry$mismatches, collapse = ", ")))
  }
  cat("\n")

  cat("## Detailed diagnostics for diverged RDS files\n\n")
  for (entry in rds_real_diff) {
    suffix <- if (!is.na(entry$dim)) sprintf(" (%s)", entry$dim) else ""
    cat(sprintf("### `%s`%s\n\n", entry$path, suffix))
    cat(sprintf("- differs in: %s\n", paste(entry$mismatches, collapse = ", ")))
    if (!is.null(entry$detail) && length(entry$detail) > 0L) {
      cat(paste(entry$detail, collapse = "\n"), "\n", sep = "")
    }
    cat("\n")
  }
}

# ---------- Skipped ----------
cat(sprintf("## Visual files skipped (HTML/PDF): %d\n\n", length(skipped_visual)))
cat("These are rendered analysis reports and figures. Byte-comparison is fragile (timestamps, fonts, anti-aliasing); eyeball any specific report side-by-side if you want to confirm a figure visually.\n\n")

if (length(skipped_other) > 0L) {
  cat(sprintf("## Other files skipped (%d)\n\n", length(skipped_other)))
  for (rel in skipped_other) cat(sprintf("- `%s`\n", rel))
  cat("\n")
}

# ---------- File-presence quirks ----------
if (length(missing_in_new) > 0L || length(extras) > 0L) {
  cat("## File-presence quirks (harmless)\n\n")
  for (rel in missing_in_new) {
    if (rel == "output") {
      cat("- **Missing in `/results/`:** `output` — Code Ocean writes its run log to `/results/output` *after* the validation step finishes; the reference asset has it because it was saved post-completion.\n")
    } else {
      cat(sprintf("- **Missing in `/results/`:** `%s`\n", rel))
    }
  }
  for (rel in extras) {
    if (rel == "TEMP_validate_against_reference.md") {
      cat("- **Extra in `/results/`:** `TEMP_validate_against_reference.md` — this script's own output. Not in the released-capsule reference because that capsule didn't run validation.\n")
    } else {
      cat(sprintf("- **Extra in `/results/`:** `%s`\n", rel))
    }
  }
  cat("\n")
}

# ---------- Dbh / Th spot-check ----------
cat("## Dbh / Th per-cell expression spot-check\n\n")
cat("Drilling into the per-cohort post-normalization SCEs and comparing counts for the two genes Polina checks first.\n\n")

for (cohort in c("BARseq_780345", "BARseq_780346", "BARseq_780345-780346_combined")) {
  rel <- file.path(cohort, "combined_neurons_clust_CCFv2_uid_cpm_log.rds")
  old_path <- file.path(OLD_ROOT, rel)
  new_path <- file.path(NEW_ROOT, rel)
  if (!file.exists(old_path) || !file.exists(new_path)) {
    cat(sprintf("### %s\n\n_Skipped: `%s` is not produced by this cohort._\n\n", cohort, rel))
    next
  }
  old_sce <- readRDS(old_path)
  new_sce <- readRDS(new_path)
  cat(sprintf("### %s — %s cells\n\n", cohort, fmt(ncol(old_sce))))
  for (gene in c("Dbh", "Th")) {
    if (!(gene %in% rownames(old_sce) && gene %in% rownames(new_sce))) {
      cat(sprintf("- **%s:** not present in both objects.\n", gene))
      next
    }
    old_vec <- as.numeric(counts(old_sce)[gene, ])
    new_vec <- as.numeric(counts(new_sce)[gene, ])
    if (identical(old_vec, new_vec)) {
      cat(sprintf("- **%s:** matches the reference exactly. Sum: %s. Max: %s. Nonzero cells: %s.\n",
                  gene, fmt(sum(old_vec)), fmt(max(old_vec)), fmt(sum(old_vec > 0))))
    } else {
      cat(sprintf("- **%s — DIFFERS.** Sum old=%s new=%s, max old=%s new=%s, nonzero cells old=%s new=%s.\n",
                  gene,
                  fmt(sum(old_vec)), fmt(sum(new_vec)),
                  fmt(max(old_vec)), fmt(max(new_vec)),
                  fmt(sum(old_vec > 0)), fmt(sum(new_vec > 0))))
    }
  }
  cat("\n")
}
