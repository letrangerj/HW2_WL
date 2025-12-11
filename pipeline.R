# Load required libraries
suppressMessages(suppressWarnings(suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
    library(jsonlite)
    library(png)
    library(dplyr)
    library(Matrix)
    library(hdf5r)
})))

print("Loading Data...")

## Helper function to read h5ad file manually using hdf5r
read_h5ad <- function(h5ad_path) {
    # Open h5ad file
    h5 <- H5File$new(h5ad_path, mode = "r")

    # Read count matrix (X)
    # h5ad files store data in /X
    if (h5$exists("X")) {
        X_group <- h5[["X"]]

        # Check if sparse or dense matrix
        if (h5$exists("X/data")) {
            # Sparse matrix (CSR format)
            print("Processing sparse matrix...")
            data <- X_group[["data"]]$read()
            indices <- X_group[["indices"]]$read()
            indptr <- X_group[["indptr"]]$read()

            # Get dimensions from shape attribute
            n_obs <- length(indptr) - 1
            shape <- X_group$attr_open("shape")$read()
            n_vars <- shape[2]

            # Convert to dgCMatrix (CSC format that Seurat expects)
            # First convert CSR to COO, then to CSC
            counts <- sparseMatrix(
                i = rep(1:n_obs, diff(indptr)),
                j = indices + 1,  # R is 1-indexed
                x = data,
                dims = c(n_obs, n_vars)
            )
            counts <- t(counts)  # Transpose to genes x cells

        } else {
            # Dense matrix
            print("Processing dense matrix...")
            counts <- t(X_group$read())
            counts <- as(counts, "sparseMatrix")
        }
    } else {
        stop("No X (count matrix) found in h5ad file")
    }

    # Read obs (cell/spot metadata)
    obs <- list()
    if (h5$exists("obs")) {
        obs_group <- h5[["obs"]]
        # Get obs index (barcodes)
        if (h5$exists("obs/_index")) {
            barcodes <- obs_group[["_index"]]$read()
        } else if (h5$exists("obs/index")) {
            barcodes <- obs_group[["index"]]$read()
        } else {
            barcodes <- paste0("cell_", 1:ncol(counts))
        }

        # Read other obs columns
        obs_names <- names(obs_group)
        obs_names <- obs_names[!obs_names %in% c("_index", "index", "__categories")]
        for (col_name in obs_names) {
            tryCatch({
                col_data <- obs_group[[col_name]]$read()
                # Handle categorical data - check if categories exist
                if (h5$exists(paste0("obs/__categories/", col_name))) {
                    categories <- h5[[paste0("obs/__categories/", col_name)]]$read()
                    col_data <- categories[col_data + 1]  # R is 1-indexed
                }
                # Only add if length matches barcodes
                if (length(col_data) == length(barcodes)) {
                    obs[[col_name]] <- col_data
                }
            }, error = function(e) {
                message(paste("Warning: Could not read column", col_name, ":", e$message))
            })
        }

        # Create data frame with barcodes as row names
        if (length(obs) > 0) {
            obs <- as.data.frame(obs, stringsAsFactors = FALSE)
            rownames(obs) <- barcodes
        } else {
            obs <- data.frame(row.names = barcodes)
        }
    } else {
        barcodes <- paste0("cell_", 1:ncol(counts))
        obs <- data.frame(row.names = barcodes)
    }

    # Read var (gene metadata)
    var <- list()
    if (h5$exists("var")) {
        var_group <- h5[["var"]]
        # Get var index (gene names)
        if (h5$exists("var/_index")) {
            genes <- var_group[["_index"]]$read()
        } else if (h5$exists("var/index")) {
            genes <- var_group[["index"]]$read()
        } else {
            genes <- paste0("gene_", 1:nrow(counts))
        }

        # Read other var columns
        var_names <- names(var_group)
        var_names <- var_names[!var_names %in% c("_index", "index", "__categories")]
        for (col_name in var_names) {
            tryCatch({
                col_data <- var_group[[col_name]]$read()
                # Handle categorical data - check if categories exist
                if (h5$exists(paste0("var/__categories/", col_name))) {
                    categories <- h5[[paste0("var/__categories/", col_name)]]$read()
                    col_data <- categories[col_data + 1]  # R is 1-indexed
                }
                # Only add if length matches genes
                if (length(col_data) == length(genes)) {
                    var[[col_name]] <- col_data
                }
            }, error = function(e) {
                message(paste("Warning: Could not read column", col_name, ":", e$message))
            })
        }

        # Create data frame with genes as row names
        if (length(var) > 0) {
            var <- as.data.frame(var, stringsAsFactors = FALSE)
            rownames(var) <- genes
        } else {
            var <- data.frame(row.names = genes)
        }
    } else {
        genes <- paste0("gene_", 1:nrow(counts))
        var <- data.frame(row.names = genes)
    }

    # Close h5 file
    h5$close()

    # Set dimnames for count matrix
    rownames(counts) <- genes
    colnames(counts) <- barcodes

    # Make gene names unique
    rownames(counts) <- make.unique(rownames(counts))

    # Create Seurat object
    seurat_obj <- CreateSeuratObject(
        counts = counts,
        meta.data = obs,
        project = "h5ad"
    )

    return(seurat_obj)
}

## Load Data
load_sample <- function(sample_id, base_dir = 'data/GSE272362/HM_files/') {
    data_dir <- file.path(base_dir, sample_id)

    # Read h5ad file using custom function
    h5ad_path <- file.path(data_dir, "expression.h5ad")

    # Load the h5ad file
    adata <- read_h5ad(h5ad_path)

    # Make gene names unique
    rownames(adata) <- make.unique(rownames(adata))

    ## Images/coordinates
    spatial_dir <- file.path(data_dir, "spatial")
    tissue_positions_path <- list.files(spatial_dir, pattern = "*_tissue_positions.csv", full.names = TRUE)[1]
    scalefactors_path <- list.files(spatial_dir, pattern = "*_scalefactors_json.json", full.names = TRUE)[1]
    image_path <- list.files(spatial_dir, pattern = "*_tissue_hires_image.png", full.names = TRUE)[1]

    # Read positions
    positions <- read.csv(tissue_positions_path)
    rownames(positions) <- positions$barcode

    # Merge spatial info into metadata (left join to keep all spots)
    adata@meta.data <- merge(adata@meta.data, positions,
                             by.x = "row.names", by.y = "barcode",
                             all.x = TRUE)
    rownames(adata@meta.data) <- adata@meta.data$Row.names
    adata@meta.data$Row.names <- NULL

    # Read scalefactors
    scalefactors <- fromJSON(scalefactors_path)

    # Read image
    image <- readPNG(image_path)

    # Store spatial information
    # Create spatial coordinates
    spatial_coords <- as.matrix(adata@meta.data[, c("pxl_row_in_fullres", "pxl_col_in_fullres")])
    rownames(spatial_coords) <- rownames(adata@meta.data)

    # Add to Seurat object
    adata[["spatial"]] <- CreateDimReducObject(
        embeddings = spatial_coords,
        key = "spatial_",
        assay = DefaultAssay(adata)
    )

    # Store image and scalefactors in misc slot
    adata@misc[["spatial"]] <- list()
    adata@misc[["spatial"]][[sample_id]] <- list(
        images = list(hires = image),
        scalefactors = scalefactors
    )

    # Add sample_id to metadata
    adata@meta.data$sample_id <- sample_id

    return(adata)
}

sample_ids <- paste0("HM33", 32:61)
adatas <- lapply(sample_ids, load_sample)

# Merge all samples
adata <- merge(adatas[[1]], y = adatas[-1],
               add.cell.ids = sample_ids,
               project = "spatial")

# Add batch information
adata@meta.data$batch <- sapply(strsplit(rownames(adata@meta.data), "_"), function(x) x[1])

print("Load data complete")
print(paste0("Loaded: ", ncol(adata), " spots, ", nrow(adata), " genes"))

# QC

## Output directories
figures_dir <- "./figures"
h5ad_dir <- "./results"

# Create directories if they don't exist
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(h5ad_dir, showWarnings = FALSE, recursive = TRUE)

# Calculate QC metrics
adata@meta.data$mt <- grepl("^MT-", rownames(adata))
adata <- PercentageFeatureSet(adata, pattern = "^MT-", col.name = "pct_counts_mt")

# Calculate total counts and n_genes_by_counts
adata@meta.data$total_counts <- colSums(GetAssayData(adata, slot = "counts"))
adata@meta.data$n_genes_by_counts <- colSums(GetAssayData(adata, slot = "counts") > 0)

## Plot processing
plot_qc <- function(adata, group) {
    # Histograms for QC metrics distribution
    png(paste0("./figures/", group, "_qc_histograms.png"),
        width = 3600, height = 1200, res = 300)
    par(mfrow = c(1, 3))

    # Total counts histogram
    hist(adata@meta.data$total_counts,
         breaks = 50,
         col = rgb(70/255, 130/255, 180/255, 0.7),
         border = "black",
         main = "Total counts per spot",
         xlab = "Total counts",
         ylab = "Frequency")
    abline(v = 3, col = "red", lty = 2, lwd = 2)
    legend("topright", legend = "Filter threshold", col = "red", lty = 2, lwd = 2)

    # Genes per spot histogram
    hist(adata@meta.data$n_genes_by_counts,
         breaks = 50,
         col = rgb(46/255, 139/255, 87/255, 0.7),
         border = "black",
         main = "Genes per spot",
         xlab = "Number of genes",
         ylab = "Frequency")
    abline(v = 3, col = "red", lty = 2, lwd = 2)
    legend("topright", legend = "Filter threshold", col = "red", lty = 2, lwd = 2)

    # Mitochondrial percentage histogram
    hist(adata@meta.data$pct_counts_mt,
         breaks = 50,
         col = rgb(255/255, 140/255, 0/255, 0.7),
         border = "black",
         main = "Mitochondrial percentage",
         xlab = "% MT counts",
         ylab = "Frequency")

    dev.off()

    # Create violin plots 
    png(paste0("./figures/", group, "_QC.png"),
        width = 4500, height = 1500, res = 300)
    par(mfrow = c(1, 3))

    # Violin plot for total_counts
    vioplot::vioplot(adata@meta.data$total_counts,
                     col = rgb(70/255, 130/255, 180/255, 0.7),
                     main = "total_counts",
                     ylab = "Counts")

    # Violin plot for n_genes_by_counts
    vioplot::vioplot(adata@meta.data$n_genes_by_counts,
                     col = rgb(46/255, 139/255, 87/255, 0.7),
                     main = "n_genes_by_counts",
                     ylab = "Number of genes")

    # Violin plot for pct_counts_mt
    vioplot::vioplot(adata@meta.data$pct_counts_mt,
                     col = rgb(255/255, 140/255, 0/255, 0.7),
                     main = "pct_counts_mt",
                     ylab = "% MT counts")

    dev.off()

    return(NULL)
}

## Preprocessing function
process_qc <- function(adata) {
    print("Processing QC metrics...")

    # Filter genes and cells
    # Filter genes: keep genes with min_counts >= 3
    gene_counts <- rowSums(GetAssayData(adata, slot = "counts"))
    adata <- adata[gene_counts >= 3, ]

    # Filter cells: min_counts >= 500
    adata <- subset(adata, subset = total_counts >= 500)

    # Filter out extreme outliers: n_genes_by_counts > 300
    adata <- subset(adata, subset = n_genes_by_counts > 300)

    # Normalization (target_sum = 1e6, equivalent to CPM)
    adata <- NormalizeData(adata,
                          normalization.method = "RC",
                          scale.factor = 1e6,
                          verbose = FALSE)

    # Log transformation (already done by NormalizeData, but we ensure log1p)
    # Seurat's NormalizeData already applies log1p, so we're good

    # PCA
    adata <- FindVariableFeatures(adata, verbose = FALSE)
    adata <- ScaleData(adata, verbose = FALSE)
    adata <- RunPCA(adata, verbose = FALSE)

    return(adata)
}

## Plot before/after processing
plot_qc(adata, "1b_before")
adata <- process_qc(adata)
plot_qc(adata, "1b_after")

# Save processed data
print("QC processing complete.")
