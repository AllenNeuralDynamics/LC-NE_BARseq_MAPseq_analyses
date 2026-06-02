# Setup: loaded at the start of every rendered script via a hidden knitr chunk.
# Do NOT load this in individual scripts — it is injected by the run script.

.libPaths("/opt/conda/envs/barseq-r4//lib/R/library/")

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(SingleCellExperiment)
  library(Seurat)
  library(MetaNeighbor)
  library(gridExtra)
  library(sctransform)
  library(clustree)
  library(viridis)
  library(dplyr)
  library(gplots)
  library(readr)
  library(pheatmap)
  library(tidyr)
  library(FNN)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  library(stringr)
  library(readxl)
  library(digest)
  library(proxy)
  library(matrixStats)
  library(grid)
  library(doParallel)
  library(cowplot)
  library(ggrepel)
  library(limma)
  library(edgeR)
  library(tibble)
})

# Set the RECOMPUTE_CLUSTERING env var to "true" (or "1", "yes") to bypass the
# committed clustering cache in code/cached_clustering/ and run clustering
# fresh. Default uses the cache. See code/cached_clustering/clustering_freeze.md.
RECOMPUTE_CLUSTERING <- tolower(Sys.getenv("RECOMPUTE_CLUSTERING", "false")) %in% c("true", "1", "yes")