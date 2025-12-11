# RCTD-Compatible Data Export for Spatial Transcriptomics

## Overview

This document describes how to prepare and analyze spatial transcriptomics data from Python/Scanpy using RCTD (Robust Cell Type Decomposition) in R.

## Workflow Summary

```
Python (Scanpy) → Raw Counts → RCTD (R) → Deconvolution Results
```

## Key Requirements

RCTD requires **untransformed integer counts** for both:
1. **Spatial transcriptomics data** (pixels/spots)
2. **Reference data** (single-cell/nuclei RNA-seq with cell type annotations)

⚠️ **Important**: Normalized or log-transformed data cannot be used with RCTD!

## What the Pipeline Does

The modified `pipeline.py` now:

1. **Loads data** from multiple Visium samples
2. **Saves raw counts** before normalization (for RCTD compatibility)
3. **Normalizes data** for standard Scanpy analysis
4. **Exports multiple formats** for different analysis workflows

## Output Files

### Standard Analysis (Normalized Data)

Located in `./results/`:

- `processed_data.h5ad` - Processed AnnData (normalized, log-transformed)
- `expression_matrix.npz` - Normalized expression matrix (sparse format)
- `cell_metadata.csv` - Cell/spot metadata
- `gene_metadata.csv` - Gene information
- `spatial_coordinates.csv` - Spatial coordinates
- `pca_coords.csv` - PCA coordinates

### RCTD Analysis (Raw Counts)

Located in `./results/rctd_input/`:

- `spatial_counts_raw.npz` - **Raw count matrix** (genes × spots) - sparse Matrix Market format
- `genes.csv` - Gene names (rownames for matrix)
- `spot_metadata.csv` - Spot metadata with barcodes (colnames for matrix)
- `spatial_coords_rctd.csv` - Spatial coordinates

## Running the Analysis

### Step 1: Run Python Pipeline

```bash
python pipeline.py
```

This will generate all the files needed for both Scanpy and RCTD analysis.

### Step 2: Run RCTD in R

```bash
Rscript rctd_analysis.R
```

Or run interactively in R:

```r
source("rctd_analysis.R")
```

### Step 3: Review Results

RCTD results are saved in `./results/rctd_output/`:

- `rctd_weights.mtx` - Cell type proportions (sparse matrix)
- `rctd_weights_full.mtx` - Unrestricted proportions
- `rctd_classifications.csv` - Spot classifications (singlet/doublet)
- `rctd_results.rds` - Full SpatialExperiment object
- `rctd_all_weights.pdf` - Visualization (pie charts)
- `rctd_{celltype}_density.pdf` - Individual cell type maps

## Reference Data

RCTD requires a **single-cell/nuclei reference dataset** with:

1. **Raw count matrix** (genes × cells)
2. **Cell type annotations** (factor vector)

### Example: Creating Reference in R

```r
# Load reference counts
reference_counts <- readMM("reference_counts.mtx")

# Load annotations
reference_cell_types <- read.csv("cell_annotations.csv")

# Create SummarizedExperiment
reference_se <- SummarizedExperiment(
  assays = list(counts = reference_counts),
  colData = DataFrame(cell_type = reference_cell_types$cell_type)
)
```

### Common Reference Sources

- **Healthy tissue atlases** (e.g., Human Cell Atlas)
- **Tumor microenvironment datasets**
- **Patient-matched scRNA-seq**
- **Public databases** (e.g., GEO, ArrayExpress)

## Understanding RCTD Modes

RCTD has three modes for different spatial technologies:

| Mode | Best For | Description |
|------|----------|-------------|
| `"doublet"` | **Visium** (55μm spots) | 1-2 cell types per spot |
| `"multi"` | Lower resolution | Up to `max_multi_types` cell types per spot |
| `"full"` | High resolution | No restrictions |

For 10x Visium data, `"doublet"` mode is **recommended**.

## Troubleshooting

### Common Issues

1. **"Reference and spatial data have no overlapping genes"**
   - Ensure gene names match between reference and spatial data
   - Check for consistent gene nomenclature (Ensembl, HGNC symbols)

2. **"Counts are not integers"**
   - Verify you're using raw counts, not normalized data
   - Use files from `./results/rctd_input/` (not `./results/`)

3. **Memory errors**
   - Reduce `max_cores` parameter in `runRctd()`
   - Try a smaller reference dataset
   - Filter spots by UMI count before running RCTD

4. **Poor deconvolution performance**
   - Increase reference dataset size
   - Ensure reference includes all expected cell types
   - Check that reference and spatial data are from similar tissues/conditions

### Adjusting Parameters

In `rctd_analysis.R`, adjust these parameters based on your data:

```r
# In createRctd()
gene_cutoff = 0.00005,     # Lower = more genes included
fc_cutoff = 1.5,           # Fold change cutoff for DE genes
UMI_min = 100,             # Minimum UMI per spot
UMI_max = 50000            # Maximum UMI per spot

# In runRctd()
rctd_mode = "doublet",     # Or "multi", "full"
max_cores = 4,             # Number of cores for parallel processing
doublet_threshold = 0.7    # Confidence threshold
```

## Integration with Scanpy

To use RCTD results in Python/Scanpy:

```python
import scanpy as sc
import pandas as pd
from scipy.io import mmread

# Load RCTD results
weights = mmread("results/rctd_output/rctd_weights_full.mtx")
classifications = pd.read_csv("results/rctd_output/rctd_classifications.csv", index_col=0)

# Load processed AnnData
adata = sc.read_h5ad("results/processed_data.h5ad")

# Add RCTD results to AnnData
adata.obsm["rctd_weights"] = weights.T  # transpose to spots × cell_types
adata.obs["rctd_class"] = classifications["spot_class"]
adata.obs["rctd_first_type"] = classifications["first_type"]

# Visualize
sc.pl.spatial(adata, color=["rctd_class", "rctd_first_type"])
```

## Performance Considerations

### Dataset Size

| Dataset Size | Approx. Runtime | Memory Required |
|--------------|----------------|----------------|
| 1,000 spots | 5-10 minutes | 4-8 GB |
| 10,000 spots | 1-2 hours | 16-32 GB |
| 50,000+ spots | 4+ hours | 64+ GB |

### Tips for Large Datasets

1. **Subset spots** by tissue region or UMI count
2. **Subset genes** to highly variable genes
3. **Use a representative reference** (not all cells needed)
4. **Increase cores** if memory allows (`max_cores = 8` or higher)
5. **Run on compute cluster** for very large datasets

## Additional Resources

- **RCTD Tutorial**: https://www.bioconductor.org/packages/spacexr
- **RCTD GitHub**: https://github.com/dmcable/spacexr
- **RCTD Paper**: Cable et al., Nature (2022)
- **Load Seurat Script**: `./load_seurat_simple.R`
- **Load Seurat Guide**: `./load_seurat.R` (comprehensive)

## Citation

If you use RCTD in your research, please cite:

```
Cable, D.M., Murray, E., Zou, L.S. et al. Robust decomposition of cell type
mixtures in spatial transcriptomics. Nat Biotechnol 40, 517–526 (2022).
https://doi.org/10.1038/s41587-021-00830-w
```

## Support

For issues with:
- **Python pipeline**: Check pipeline.py comments or file an issue
- **RCTD deconvolution**: See RCTD documentation or GitHub issues
- **Data preparation**: Compare your data format to this template

## Quick Reference Checklist

- [ ] Run `python pipeline.py` ✓
- [ ] Verify files in `./results/rctd_input/` ✓
- [ ] Prepare reference dataset ✓
- [ ] Update `rctd_analysis.R` with your reference ✓
- [ ] Run `Rscript rctd_analysis.R` ✓
- [ ] Review results in `./results/rctd_output/` ✓
- [ ] Visualize and interpret ✓

---

**Version**: 1.0  
**Last Updated**: 2025-12-11  
**Compatible with**: spacexr ≥ 1.2.0, Scanpy ≥ 1.9.0
