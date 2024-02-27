#use barseq-r4
.libPaths("/opt/conda/envs/barseq-r4//lib/R/library/")
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(Seurat)
library(MetaNeighbor)

library(googlesheets4)
library(sctransform)
library(clustree)

#example - reinstall a package under barseq-r4
install.packages("reticulate")
library(reticulate)
#use envs 
conda_list()
system("ldconfig /usr/local/lib")
use_python("/Python-3.7.12/python")
use_condaenv("rstudio-r4-base",required = FALSE)
use_condaenv("barseq-r4",required = FALSE)
