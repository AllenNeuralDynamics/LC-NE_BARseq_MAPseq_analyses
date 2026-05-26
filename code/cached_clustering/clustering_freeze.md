# Frozen clustering snapshot

The CSVs in this directory are frozen UMAP layouts + Louvain cluster assignments produced by `analyze_barseq()` in `code/1_BARseq_analyses_functions_*.R`. They live here in the repo (rather than being recomputed at runtime) so the analysis pipeline produces the same figures across rebuilds.

The three subdirectories correspond to the three analysis cohorts: `BARseq_780345/` (brain3), `BARseq_780346/` (brain4), and `BARseq_780345-780346_combined/` (combined).

## Where they came from

All three sets of CSVs were produced by a single computation in a separate test-build capsule, then extracted from the resulting Code Ocean data assets and committed here. Shared provenance:

- **Capsule:** `PK_BARseq_MAPseq_LC-NE_test_build`
- **Code version:** `e5f6e52`
- **Computation ID:** `b9ee59ed-ef4b-496e-8654-a68e8cdb64ca`
- **Frozen on:** 2026-04-03 by Polina Kosillo

The three Code Ocean data assets the CSVs were extracted from:

| Subdirectory | CO asset name (display) | CO asset ID |
|---|---|---|
| `BARseq_780345/` | `BARseq_780345_UMAPclustering` | `a78668f6-3451-42d4-9a83-d88f6884e83d` |
| `BARseq_780346/` | `BARseq_780346_UMAPclustering` | `0e3f1b4f-2f63-4692-acfe-1411d3f44e25` |
| `BARseq_780345-780346_combined/` | `BARseq_780345-780346_combined_UMAPclustering` | `546c49f5-56c0-45fa-9787-6582332482a5` |

Only the `umap.csv`, `cluster.csv`, and (for the combined cohort) `cluster_annotation.csv` files are committed — the analysis code reads only those. Other folders that exist in the source data assets (e.g. `Dbh_corrected_test/`, `fromLCNE_combined/`) are not referenced by any analysis script and are not copied here.

## Why these are committed rather than regenerated at runtime

The first round of clustering operates on ~3M QCed cells and takes many hours of compute on a standard CO machine. Committing the resulting CSVs avoids that delay on every Reproducible Run and guarantees that the pipeline produces consistent figures across rebuilds since every run uses the same cluster IDs and UMAP coordinates as the original analysis.

## Forcing a fresh clustering run

To bypass the cached CSVs and recompute clustering from scratch, set the `RECOMPUTE_CLUSTERING` environment variable to `true` (or `1`, or `yes`) before triggering a Reproducible Run. On Code Ocean this is set in the capsule's Environment Variables (gear icon → Environment Variables). The flag is read once at the top of `code/setup.R` and short-circuits every cache-check `if` block in the three `code/2_BARseq_norm_cluster_analyze_*.R` scripts.

Plan for a long run: clustering all ~3M QCed cells in the first round takes many hours of compute.
