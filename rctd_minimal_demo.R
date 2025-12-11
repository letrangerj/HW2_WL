#!/usr/bin/env Rscript
#' Minimal RCTD Demo with Placeholder Reference
#' This shows the minimal code needed to run RCTD once you have a reference

cat("Minimal RCTD Demo\n")
cat("==================\n\n")

# Load libraries
suppressMessages(library(spacexr))
suppressMessages(library(SpatialExperiment))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(Matrix))

# Load spatial data
cat("Loading spatial transcriptomics data...\n")
counts_matrix <- readMM("results/rctd_input/spatial_counts_raw.npz")
genes <- read.csv("results/rctd_input/genes.csv", header = FALSE)[,1]
spot_metadata <- read.csv("results/rctd_input/spot_metadata.csv", row.names = 1)

rownames(counts_matrix) <- genes
colnames(counts_matrix) <- rownames(spot_metadata)

spatial_spe <- SpatialExperiment(assays = list(counts = counts_matrix))

cat(sprintf("Loaded: %d genes x %d spots\n\n", nrow(spatial_spe), ncol(spatial_spe)))

# --- IMPORTANT: Replace with your actual reference! ---
cat("⚠️  USING PLACEHOLDER REFERENCE - REPLACE WITH YOUR DATA!\n\n")

# Simulate small reference
set.seed(42)
gene_subset <- sample(genes, 300)
reference_counts <- matrix(rpois(300 * 50, lambda = 5), nrow = 300, 
                          dimnames = list(gene_subset, paste0("cell_", 1:50)))

reference_cell_types <- factor(sample(c("T_cell", "B_cell", "Macrophage", "Tumor"), 50, replace = TRUE))

reference_se <- SummarizedExperiment(
  assays = list(counts = reference_counts),
  colData = DataFrame(cell_type = reference_cell_types)
)

cat("Reference dimensions:", dim(reference_se), "\n")
cat("Cell types:", paste(unique(reference_cell_types), collapse = ", "), "\n\n")
# --------------------------------------------------------

# Run RCTD
cat("Running RCTD (this may take a few minutes)...\n\n")

# Preprocess
rctd_data <- createRctd(spatial_spe, reference_se, max_cores = 2)

# Deconvolve
results_spe <- runRctd(rctd_data, rctd_mode = "doublet", max_cores = 2)

# Save results
dir.create("results/rctd_output", showWarnings = FALSE)
saveRDS(results_spe, "results/rctd_output/rctd_demo_results.rds")

# Show results summary
cat("Results summary:\n")
cat(sprintf("Spots analyzed: %d\n", ncol(results_spe)))
cat("Classifications:\n")
print(table(results_spe$spot_class))

cat("\nCell types identified:\n")
print(sort(rowSums(assay(results_spe, "weights_full")), decreasing = TRUE)[1:5])

cat("\n✓ Demo complete! Results saved to results/rctd_output/\n")
cat("\nNext: Replace placeholder reference with your actual scRNA-seq data!\n")
