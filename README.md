# LC-NE BARseq and MAPseq Analyses

Code for analyzing BARseq and MAPseq projection data from locus coeruleus norepinephrine (LC-NE) neurons, as described in:

> Su, Kosillo, Jung, Chen *et al.* (2026). Topographic structure and function of locus coeruleus norepinephrine neurons. [bioRxiv 2026.04.10.717727](https://www.biorxiv.org/content/10.64898/2026.04.10.717727v1)

The capsule processes BARseq gene-expression barcoding and MAPseq projection-mapping data from two specimens (780345 and 780346) to identify LC-NE neuron subtypes and characterize their projection patterns. Outputs feed **Figure S5** of the manuscript.

**GitHub:** https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAPseq_analyses
**Code Ocean:** https://codeocean.allenneuraldynamics.org/capsule/2195789/tree
**Collection:** https://codeocean.allenneuraldynamics.org/collections/9cf044ce-93c7-4c7e-bfa1-5d8c37aa42ec

## Running

Click **Reproducible Run** in Code Ocean. The `run` script renders each numbered analysis stage to a self-contained HTML report (~1–2 hours on a large instance).

## Code

- `setup.R`, `00_env_lib_loading.R` — load R libraries
- `01_loaders_*.R`, `02_prepare_brain3_4_combined_inputs.R` — per-brain data loaders + combined-brain prep
- `1_BARseq_analyses_functions_*.R` — shared functions (normalization, clustering, spatial coherence)
- `2_BARseq_norm_cluster_analyze_*.R` — normalize counts, cluster, identify LC-NE cells
- `3_MAPseq_match_BARseq_*.R` — match BARseq barcodes to MAPseq barcodes (Hamming distance)
- `4_MAPseq_Klebschull_replicate_CTX_proj_*.R` — replicate Bhatt/Kebschull et al. (2022) cortical projection analysis
- `5_MAPseq_probability_*.R` — projection probabilities, heatmaps, co-innervation
- `6_MAPseq_ExA-SPIM_*.R` — comparison with ExA-SPIM single-neuron morphology

Each numbered stage has three variants: `_brain3.R`, `_brain4.R`, `_brain3-4_combined.R`.

Clustering UMAPs + cluster-label CSVs are committed under `code/cached_clustering/` and reloaded by default — see `code/cached_clustering/clustering_freeze.md` for provenance. Set `RECOMPUTE_CLUSTERING=true` in the capsule's environment variables to recompute clustering from scratch.

## Inputs

Four data assets are attached, one BARseq + one MAPseq per specimen. Files within each are loaded directly from `/data/<asset_mount>/...` by the analysis scripts.

### BARseq per-specimen assets (`780345_..._processed-MAT2RDS_...`, `780346_..._processed-MAT2RDS_...`)

Outputs of the `LC-NE_BARseq_MAT-RDS_conversion` capsule, which converts the upstream MATLAB BARseq pipeline outputs into R-friendly formats ([GitHub](https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAT-RDS_conversion), [Code Ocean](https://codeocean.allenneuraldynamics.org/capsule/3953531)). All files live under `BARseq/` inside each asset.

| File | Description |
|---|---|
| `combined_neurons_clust_CCFv2_uid.rds` | `SingleCellExperiment` of all QCed BARseq cells for the specimen (~300–500 K cells × 103 genes). Contains the raw count matrix plus per-cell `colData` columns: CCF coordinates (`CCF_AP`, `CCF_DV`, `CCF_ML`, `CCFano`), slice index, imaging-FOV coordinates, somatic-barcode index, batch, and a unique cell id (`uid`). Loaded at the top of stage 2; the entry point for the whole pipeline. |
| `combined_neurons_clust_CCFv2.rds` | Predecessor of the above without `uid` assignment. Not used; superseded. |
| `DBHfilteredneurons_clust_CCFv2_uid.rds` | An earlier `Dbh`-positive-only subset, kept for historical context. Not used by this pipeline. |
| `barcodes_BC_qc_<subject>.csv` | Per-cell BARseq somatic barcode sequences (15 nt) for cells that passed barcode QC. Joined to MAPseq projection barcodes via Hamming-distance matching in stage 3. |
| `LC_visualQC_barcoded_cells_<subject>.csv` | Manual visual-QC annotations of barcoded LC-NE cells (`uid` + QC flags). Used in stages 3 and 4 to restrict barcode matching and projection analyses to cells that passed visual QC. |

### MAPseq per-specimen assets (`mapseq_780345_2025-03-24_12-00-00`, `mapseq_780346_2025-07-23_12-00-00`)

Raw MAPseq projection-barcode counts and dissection metadata. The relevant files live in two places inside each asset:

| File | Location | Description |
|---|---|---|
| `<subject>.nbcm.tsv` | `MAPseq/M<run>_<date>_USEthis/` | Filtered (background-subtracted, spike-in-normalized) MAPseq UMI count matrix — rows = projection barcodes, columns = ROIs (`BC*`) and a soma column. The primary MAPseq input for downstream matching. |
| `<subject>.rbcm.tsv` | same | Raw MAPseq UMI count matrix (pre-filter). Used in stage 3 QC checks only. |
| `<subject>.sbcm.tsv` | same | Spike-in barcode counts. Used in stage 3 QC checks only. |
| `M<run>sampleinfo.tsv` or `M<run>_<date>.sampleinfo.xlsx` | same | Per-tube experiment metadata (tube number, dissection labels, processing notes). The brain3 asset has both; brain4 has only the xlsx. |
| `sampleinfo_<subject>.tsv` (brain3) / `sampleinfo_<subject>.xlsx` (brain4) | `MAPseq/` (asset root) | Curated lookup table mapping MAPseq sample-tube numbers (`BC*`) to CCF brain-region names + dissection metadata. Used in stages 1/4 to label projection columns with region names. |

## Outputs

After all analysis stages run, a final step (`code/07_collect_paper_figures.R`) reorganizes `/results/` for publication:

```
results/
  <stage>.html              # rendered analysis report, one per stage
  paper_figures/FigureS5/   # the manuscript panels, named by figure
  other_results/            # per-cohort data, QC plots, intermediate CSVs
    BARseq_780345/                  (brain3, specimen 780345)
    BARseq_780346/                  (brain4, specimen 780346)
    BARseq_780345-780346_combined/  (combined cross-brain analyses)
  data_description.json     # AIND derived-asset metadata (code/08_generate_metadata.py)
  processing.json
```

### Figure S5: LC-NE projections measured with MAPseq and BARseq

This capsule produces three panels of manuscript Figure S5 (confirmed by the authors). All come from the combined cross-brain analyses and are copied into `paper_figures/FigureS5/` under canonical names:

| Panel | Published file | Source file (in `other_results/BARseq_780345-780346_combined/`) | Stage |
|-------|----------------|------------------------------------------------------------------|-------|
| S5e | `FigureS5e_ipsi_contra_projection_heatmap.pdf` | `Combined_ipsi-contra_projections_heatmap_top_region_sorted.pdf` | 5 |
| S5f | `FigureS5f_sorted_projection_heatmap_ipsi_contra.pdf` | `sorted_proj_heatmap_ipsi-contra.pdf` | 6 |
| S5g | `FigureS5g_rank_correlation_ipsi_contra.pdf` | `rank_corr_ipsi-contra.pdf` | 6 |

Other Figure S5 panels are not produced here: S5a (schema) and S5b (dissection boundaries) are hand-drawn schematics, and the remaining panels come from other capsules. Everything in `other_results/` is exploratory / QC output and intermediate data, not manuscript figures.

### Publishing the results as a data asset

`code/08_generate_metadata.py` writes `data_description.json` and `processing.json` into `/results/` so the run output can be saved as a DERIVED asset in `aind-open-data` (a downstream capsule consumes it as input). Provenance (capsule URL + release version) is pulled from the Code Ocean API at runtime, which requires the **"Code Ocean API Credentials" Secret** attached to the capsule (Capsule Settings → Credentials). A start-of-run preflight (`code/00_check_credentials.py`) verifies these credentials before the long pipeline runs.

## Environment

R 4.3.0 with Seurat, SingleCellExperiment, scater, scran, MetaNeighbor, and ~30 additional packages. Pinned versions in `environment/barseq-r4.yml`; consumed by `environment/Dockerfile`.

## License

MIT — see [LICENSE](LICENSE).
