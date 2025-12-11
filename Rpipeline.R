suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(tidyverse)
  library(magick)      # For image processing
  library(jsonlite)    # For JSON handling
  library(cowplot)     # For plot arrangement
  library(here)        # For path handling
  library(SpatialExperiment)
  library(SummarizedExperiment)
  library(spacexr)     # Main RCTD package
  library(Matrix)
})


## Load reference scRNA-seq data
reference <- readRDS("./data/scRNA-seq_Data_post_qc.rds")

print("Reference Celltype Annotations...")
print(table(reference$celltypes))

## Extract counts/celltypes
counts <- GetAssayData(reference, slot = "counts")
cell_types <- reference$celltypes
names(cell_types) <- colnames(reference)

## Create SummarizedExperiment for RCTD reference
reference_se <- SummarizedExperiment(
  assays = list(counts = counts),
  colData = DataFrame(cell_type = cell_types, nUMI = reference$nCount_RNA)
)

print("Reference SummarizedExperiment:")
print(dim(reference_se))
print(table(colData(reference_se)$cell_type))

# Load spatial data from pipeline.py output
print("Loading ST data...")
counts_matrix <- readMM("results/rctd_input/spatial_counts_raw.npz")
genes <- read.csv("results/rctd_input/genes.csv", header = FALSE)[,1]
spot_metadata <- read.csv("results/rctd_input/spot_metadata.csv", row.names = 1)

rownames(counts_matrix) <- genes
colnames(counts_matrix) <- rownames(spot_metadata)

spatial_spe <- SpatialExperiment(assays = list(counts = counts_matrix))

sprintf("Loaded: %d genes x %d spots\n\n", nrow(spatial_spe), ncol(spatial_spe))

## Preprocess data for RCTD
print("Preprocessing data for RCTD...")
rctd_data <- createRctd(spatial_spe, reference_se)

## Run RCTD
print("Running RCTD...")
results_spe <- runRctd(rctd_data, rctd_mode = "doublet", max_cores = 1)

print("RCTD completed successfully!")
print("Results available in results_spe object")
