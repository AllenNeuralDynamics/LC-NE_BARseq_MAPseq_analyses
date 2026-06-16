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

This capsule expects **four data assets** — one BARseq and one MAPseq asset per specimen (780345 = "brain 3", 780346 = "brain 4"). The assets are detached in `.codeocean/datasets.json`; attach them in Code Ocean before a run. The analysis scripts hard-code the `/data/<mount>/...` paths, so each asset must be mounted under the exact name below:

| Asset (mount name) | Modality | Specimen | Source |
|---|---|---|---|
| `780345_2025-02-24_12-00-00_processed-MAT2RDS_2026-06-16_00-49-36` | BARseq | 780345 (brain 3) | derived |
| `780346_2025-06-13_12-00-00_processed-MAT2RDS_2026-06-16_00-51-19` | BARseq | 780346 (brain 4) | derived |
| `780345_2025-03-24_12-00-00` | MAPseq | 780345 (brain 3) | raw |
| `780346_2025-07-23_12-00-00` | MAPseq | 780346 (brain 4) | raw |

The two BARseq assets are outputs of the [`LC-NE_BARseq_MAT-RDS_conversion`](https://github.com/AllenNeuralDynamics/LC-NE_BARseq_MAT-RDS_conversion) capsule ([Code Ocean](https://codeocean.allenneuraldynamics.org/capsule/3953531)), which converts the upstream MATLAB BARseq pipeline outputs into R-friendly formats. Their mount names embed the conversion run's creation timestamp, so they change if the conversion capsule is re-run. The two MAPseq assets are raw projection-barcode counts and dissection metadata.

The files actually consumed by the pipeline, listed under the asset each comes from:

### `780345_..._processed-MAT2RDS_...` — BARseq, brain 3

Paths relative to `BARseq/` within the asset.

| File | Description |
|---|---|
| `combined_neurons_clust_CCFv2_uid.rds` | `SingleCellExperiment` of all QCed BARseq cells for the specimen (~300–500 K cells × 103 genes). Raw count matrix plus per-cell `colData`: CCF coordinates (`CCF_AP`, `CCF_DV`, `CCF_ML`, `CCFano`), slice index, imaging-FOV coordinates, somatic-barcode index, batch, and a unique cell id (`uid`). Loaded at the top of stage 2 — the entry point for the whole pipeline. |
| `barcodes_BC_qc_780345.csv` | Per-cell BARseq somatic barcode sequences (15 nt) for cells that passed barcode QC. Joined to MAPseq projection barcodes via Hamming-distance matching in stage 3. |
| `LC_visualQC_barcoded_cells_780345.csv` | Manual visual-QC annotations of barcoded LC-NE cells (`uid` + QC flags). Used in stages 3 and 4 to restrict matching and projection analyses to cells that passed visual QC. |

The asset also contains `combined_neurons_clust_CCFv2.rds` (same without `uid`, superseded) and `DBHfilteredneurons_clust_CCFv2_uid.rds` (an earlier `Dbh`-positive-only subset). Neither is used by this pipeline.

### `780346_..._processed-MAT2RDS_...` — BARseq, brain 4

Same layout as the brain-3 BARseq asset, with `780346` in place of `780345`. Paths relative to `BARseq/`.

| File | Description |
|---|---|
| `combined_neurons_clust_CCFv2_uid.rds` | As above, for specimen 780346. Entry point for stage 2. |
| `barcodes_BC_qc_780346.csv` | As above, for specimen 780346. |
| `LC_visualQC_barcoded_cells_780346.csv` | As above, for specimen 780346. |

Also contains the unused `combined_neurons_clust_CCFv2.rds` and `DBHfilteredneurons_clust_CCFv2_uid.rds`.

### `780345_2025-03-24_12-00-00` — MAPseq, brain 3

| File | Description |
|---|---|
| `MAPseq/M295_20250729_USEthis/780345.nbcm.tsv` | Filtered (background-subtracted, spike-in-normalized) MAPseq UMI count matrix — rows = projection barcodes, columns = ROIs (`BC*`) and a soma column. The primary MAPseq input for downstream matching (stage 3). |
| `MAPseq/M295_20250729_USEthis/780345.rbcm.tsv` | Raw MAPseq UMI count matrix (pre-filter). Stage 3 QC checks only. |
| `MAPseq/M295_20250729_USEthis/780345.sbcm.tsv` | Spike-in barcode counts. Stage 3 QC checks only. |
| `MAPseq/M295_20250729_USEthis/M295_20250721.sampleinfo.xlsx` | Per-tube experiment metadata (tube number, dissection labels, processing notes). Read by the stage-1 loaders. |
| `MAPseq/sampleinfo_780345.tsv` | Curated lookup mapping sample-tube numbers (`BC*`) to CCF brain-region names + dissection metadata. Labels projection columns with region names (stages 1/4). |

### `780346_2025-07-23_12-00-00` — MAPseq, brain 4

Same structure as the brain-3 MAPseq asset, but the count-matrix files carry a `1025` suffix and the per-tube metadata is a `.tsv`.

| File | Description |
|---|---|
| `MAPseq/M305_20251030_USEthis/780346.nbcm1025.tsv` | Filtered MAPseq UMI count matrix (primary input, stage 3). |
| `MAPseq/M305_20251030_USEthis/780346.rbcm1025.tsv` | Raw MAPseq UMI count matrix (pre-filter). Stage 3 QC checks only. |
| `MAPseq/M305_20251030_USEthis/780346.sbcm1025.tsv` | Spike-in barcode counts. Stage 3 QC checks only. |
| `MAPseq/M305_20251030_USEthis/M305sampleinfo.tsv` | Per-tube experiment metadata. Read by the stage-1 loaders. |
| `MAPseq/sampleinfo_780346.xlsx` | Curated tube → CCF-region lookup + dissection metadata (stages 1/4). |

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
