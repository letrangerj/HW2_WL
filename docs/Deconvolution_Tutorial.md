# Cell Type Deconvolution Tutorial for Spatial Transcriptomics

## Table of Contents
1. [Overview](#overview)
2. [Understanding Deconvolution](#understanding-deconvolution)
3. [RCTD Method (Paper's Approach)](#rctd-method)
4. [Python Alternatives](#python-alternatives)
5. [Complete Workflows](#complete-workflows)
6. [Visualization](#visualization)
7. [Common Issues](#common-issues)

---

## Overview

### What is Deconvolution?

In Visium spatial transcriptomics:
- **Problem**: Each spot (~55μm diameter) contains multiple cells (1-10+ cells)
- **Spot-by-gene matrix**: Expression represents mixed signal from all cells in that spot
- **Goal**: Decompose the mixed signal to estimate cell-type proportions

**Analogy**: Like listening to an orchestra and trying to identify which instruments are playing and their relative volumes.

### Why Do We Need It?

Without deconvolution:
- Can't determine cell-type composition
- Can't link spatial patterns to specific cell types
- Can't study cell-cell interactions at cellular resolution

---

## Understanding Deconvolution

### Basic Principle

```
Observed expression in spot = Σ (cell_type_fraction × cell_type_expression_profile)
```

**Example:**
```
Spot Expression:
  Gene A: 100 counts
  Gene B: 50 counts
  Gene C: 200 counts

Possible Decomposition:
  30% Tumor cells (high Gene A, B, C)
  50% Fibroblasts (high Gene C, low A, B)
  20% Immune cells (high Gene B, moderate A, C)
```

### Types of Deconvolution

**1. Reference-Based (Supervised)**
- **Input**: scRNA-seq reference + spatial data
- **Method**: Use known cell-type profiles to decompose spots
- **Advantages**: More accurate, cell-type specific
- **Examples**: RCTD, Cell2location, Tangram, SPOTlight

**2. Reference-Free (Unsupervised)**
- **Input**: Spatial data only
- **Method**: Discover cell states from spatial patterns
- **Advantages**: No reference needed
- **Examples**: NMF, SpatialDE, STdeconvolve

**For this assignment, use reference-based methods (RCTD as in the paper).**

---

## RCTD Method

### Overview

**RCTD (Robust Cell Type Decomposition)** was used in the original paper.

**Key Features:**
- Developed specifically for spatial transcriptomics
- Handles platform differences (10x scRNA-seq vs Visium)
- Accounts for doublets (multiple cell types per spot)
- Provides confidence scores

**Citation**: Cable et al., Nature Biotechnology 2022

### RCTD Modes

1. **"doublet" mode**: Assigns 1-2 dominant cell types per spot
2. **"full" mode**: Provides fractional abundances for all cell types (used in paper)
3. **"multi" mode**: Allows multiple cell types with weights

**The paper used "full mode"** to get complete cell-type fractions.

### Reference Data Required

From the paper (Methods section):
- **PDAC samples**: Peng et al. 2019 (Cell Research)
- **Metastatic samples**: Raghavan et al. 2021 (Cell)
- **Normal liver**: MacParland et al. 2018 (Nature Communications)
- **Lymph node**: Abe et al. 2022 (Nature Cell Biology)

**Download link** (from assignment): https://zenodo.org/records/10712047/files/scRNA-seq_Data_post_qc.rds

---

## Complete Workflows

### Option 1: RCTD in R (Exact Paper Method)

#### Step 1: Install Required Packages

```r
# Install RCTD
if (!require("devtools")) install.packages("devtools")
devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)

# Other packages
install.packages(c("Seurat", "Matrix", "ggplot2"))
```

#### Step 2: Prepare Reference scRNA-seq Data

```r
library(Seurat)
library(spacexr)
library(Matrix)

# Load reference scRNA-seq data (downloaded from Zenodo)
reference <- readRDS("path/to/scRNA-seq_Data_post_qc.rds")

# The reference should be a Seurat object with cell-type annotations
# Check structure
print(reference)
print(table(reference$cell_type))  # Check cell-type annotations

# Extract counts and cell types for RCTD
counts <- GetAssayData(reference, slot = "counts")
cell_types <- reference$cell_type  # Or whatever column contains annotations
names(cell_types) <- colnames(reference)

# Create RCTD reference object
reference_rctd <- Reference(
  counts = counts,
  cell_types = cell_types,
  nUMI = reference$nCount_RNA  # Total UMI counts per cell
)

print(reference_rctd)
```

#### Step 3: Prepare Spatial Data

```r
# Load your processed spatial data
adata_path <- "results/processed_data.h5ad"

# Option A: If using Python h5ad, convert to Seurat
library(reticulate)
sc <- import("scanpy")
adata <- sc$read_h5ad(adata_path)

# Convert to Seurat (or load directly if already in Seurat format)
# This is a simplified example - you may need to adjust based on your data structure

# Option B: Load directly if you saved as Seurat object
# spatial_data <- readRDS("path/to/spatial_seurat.rds")

# For RCTD, we need:
# 1. Counts matrix
# 2. Spatial coordinates
# 3. Number of UMIs per spot

# Extract from your spatial data
spot_counts <- GetAssayData(spatial_data, slot = "counts")
spot_coords <- GetTissueCoordinates(spatial_data)
spot_nUMI <- colSums(spot_counts)

# Create RCTD SpatialRNA object
spatial_rctd <- SpatialRNA(
  coords = spot_coords,
  counts = spot_counts,
  nUMI = spot_nUMI
)

print(spatial_rctd)
```

#### Step 4: Run RCTD Deconvolution

```r
# Create RCTD object
rctd_obj <- create.RCTD(
  spatialRNA = spatial_rctd,
  reference = reference_rctd,
  max_cores = 4  # Adjust based on your system
)

# Run deconvolution in FULL mode (as in paper)
rctd_obj <- run.RCTD(
  RCTD = rctd_obj,
  doublet_mode = "full"  # Get all cell-type fractions
)

# Extract results
results <- rctd_obj@results

# Get cell-type proportions
weights <- results$weights  # Cell-type fractions for each spot

# Get per-cell-type weights normalized to sum to 1
norm_weights <- normalize_weights(weights)

print(dim(norm_weights))
print(head(norm_weights))
```

#### Step 5: Visualize Results

```r
library(ggplot2)

# Add deconvolution results to spatial object
spatial_data@meta.data <- cbind(
  spatial_data@meta.data,
  norm_weights[rownames(spatial_data@meta.data), ]
)

# Plot individual cell types
cell_types_to_plot <- colnames(norm_weights)

for (ct in cell_types_to_plot) {
  p <- SpatialFeaturePlot(
    spatial_data,
    features = ct,
    pt.size.factor = 1.6,
    alpha = c(0.1, 1)
  ) +
    scale_fill_gradientn(
      colors = c("lightgray", "blue", "red"),
      na.value = "transparent"
    ) +
    ggtitle(paste("Proportion of", ct))

  print(p)
  ggsave(paste0("figures/deconv_", ct, ".png"), p, width = 8, height = 6)
}

# Plot dominant cell type per spot
spatial_data$dominant_celltype <- colnames(norm_weights)[
  apply(norm_weights, 1, which.max)
]

p <- SpatialDimPlot(
  spatial_data,
  group.by = "dominant_celltype",
  pt.size.factor = 1.6
) +
  ggtitle("Dominant Cell Type per Spot")

print(p)
ggsave("figures/deconv_dominant_celltype.png", p, width = 10, height = 8)
```

#### Step 6: Validate with Marker Genes

```r
# Plot marker genes alongside deconvolution results
marker_genes <- list(
  "Tumor" = c("EPCAM", "KRT8", "KRT18"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN"),
  "Immune" = c("PTPRC", "CD3D", "CD68"),
  "Endothelial" = c("PECAM1", "VWF", "CLDN5")
)

for (ct in names(marker_genes)) {
  genes <- marker_genes[[ct]]

  # Plot markers
  p1 <- SpatialFeaturePlot(
    spatial_data,
    features = genes,
    ncol = length(genes),
    pt.size.factor = 1.6
  )

  # Plot deconvolution proportion
  p2 <- SpatialFeaturePlot(
    spatial_data,
    features = paste0(ct, "_proportion"),  # Adjust to match your column names
    pt.size.factor = 1.6
  )

  # Combine plots
  combined <- p1 / p2
  ggsave(
    paste0("figures/validation_", ct, ".png"),
    combined,
    width = 12,
    height = 10
  )
}
```

---

### Option 2: Python Alternatives

If you prefer Python (assignment allows using alternative algorithms):

#### Option A: Cell2location

**Advantages:**
- Python-native
- Probabilistic framework
- Handles technical artifacts well

```python
import scanpy as sc
import cell2location
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load reference scRNA-seq data
# Note: You'll need to convert the .rds file to h5ad first
adata_ref = sc.read_h5ad("reference_scrna.h5ad")

# Load spatial data
adata_spatial = sc.read_h5ad("results/processed_data.h5ad")

# Prepare reference (similar to RCTD)
cell2location.models.RegressionModel.setup_anndata(
    adata=adata_ref,
    labels_key='cell_type'  # Column with cell-type annotations
)

# Train reference model
ref_model = cell2location.models.RegressionModel(adata_ref)
ref_model.train(max_epochs=250)

# Export estimated expression per cell type
adata_ref = ref_model.export_posterior(
    adata_ref,
    sample_kwargs={'num_samples': 1000}
)

# Prepare spatial data
cell2location.models.Cell2location.setup_anndata(adata=adata_spatial)

# Create cell2location model
model = cell2location.models.Cell2location(
    adata_spatial,
    cell_state_df=adata_ref.varm['means_per_cluster_mu_fg'],
    N_cells_per_location=30,  # Expected cells per spot
    detection_alpha=20
)

# Train
model.train(max_epochs=30000, batch_size=None, train_size=1)

# Export results
adata_spatial = model.export_posterior(
    adata_spatial,
    sample_kwargs={'num_samples': 1000}
)

# Extract cell-type proportions
cell_type_props = adata_spatial.obsm['q05_cell_abundance_w_sf']

# Add to adata
for i, ct in enumerate(adata_ref.obs['cell_type'].cat.categories):
    adata_spatial.obs[f'{ct}_proportion'] = cell_type_props[:, i]

# Save results
adata_spatial.write("results/adata_deconvolved.h5ad")

# Visualize
import squidpy as sq

for ct in adata_ref.obs['cell_type'].cat.categories:
    sq.pl.spatial_scatter(
        adata_spatial,
        color=f'{ct}_proportion',
        spot_size=150,
        save=f"deconv_{ct}.png"
    )
```

#### Option B: SpatialDecon (R, alternative)

```r
library(SpatialDecon)

# This is lighter weight than RCTD but less robust
# Good for quick analysis

# Prepare reference
ref_matrix <- GetAssayData(reference, slot = "data")  # Log-normalized
cell_types <- reference$cell_type

# Create cell profile matrix
cell_profiles <- create_profile_matrix(ref_matrix, cell_types)

# Run deconvolution
results <- spatialdecon(
  norm = GetAssayData(spatial_data, slot = "data"),
  bg = NULL,
  X = cell_profiles,
  align_genes = TRUE
)

# Extract proportions
props <- results$prop_of_all

# Add to Seurat object
spatial_data@meta.data <- cbind(
  spatial_data@meta.data,
  t(props)
)
```

#### Option C: Tangram (Python, mapping-based)

```python
import tangram as tg
import scanpy as sc

# Load data
adata_ref = sc.read_h5ad("reference_scrna.h5ad")
adata_spatial = sc.read_h5ad("results/processed_data.h5ad")

# Map cells to space
tg.pp_adatas(adata_ref, adata_spatial, genes=None)

# Train mapping
ad_map = tg.map_cells_to_space(
    adata_sc=adata_ref,
    adata_sp=adata_spatial,
    mode='cells',  # Map individual cells
    device='cpu'
)

# Project cell types to space
tg.project_cell_annotations(
    ad_map,
    adata_spatial,
    annotation='cell_type'
)

# The result is in adata_spatial.obsm['tangram_ct_pred']
# Contains predicted cell-type fractions per spot

# Convert to proportions
cell_type_props = adata_spatial.obsm['tangram_ct_pred']
cell_type_props = cell_type_props / cell_type_props.sum(axis=1, keepdims=True)

# Add to obs
for i, ct in enumerate(adata_ref.obs['cell_type'].cat.categories):
    adata_spatial.obs[f'{ct}_proportion'] = cell_type_props[:, i]
```

---

## Visualization

### Essential Visualizations for Assignment Task 2b

#### 1. Cell-Type Marker Expression (Fig. 1c from paper)

```python
import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq

adata = sc.read_h5ad("results/adata_deconvolved.h5ad")

# Define markers
markers = {
    'Tumor': ['EPCAM', 'KRT8', 'KRT18'],
    'Fibroblasts': ['COL1A1', 'COL1A2', 'DCN'],
    'Immune_T': ['CD3D', 'CD4', 'CD8A'],
    'Immune_B': ['CD79A', 'MS4A1'],
    'Macrophages': ['CD68', 'CD163', 'LYZ'],
    'Endothelial': ['PECAM1', 'VWF', 'CLDN5']
}

# Plot markers spatially
fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.flatten()

for idx, (celltype, genes) in enumerate(markers.items()):
    # Average expression of markers
    marker_expr = adata[:, genes].X.mean(axis=1).A1

    sq.pl.spatial_scatter(
        adata,
        color=marker_expr,
        title=f'{celltype} markers',
        spot_size=150,
        ax=axes[idx],
        colorbar_loc=None
    )

plt.tight_layout()
plt.savefig("figures/marker_expression_spatial.png", dpi=300, bbox_inches='tight')
plt.close()
```

#### 2. Deconvolution Results (Fig. 1d from paper)

```python
# Plot cell-type proportions
cell_types = ['Tumor_epithelial_cells', 'myCAF', 'iCAF', 'B_cells',
              'CD4_cells', 'CD8_NK_cells', 'FCN1_TAM', 'C1Q_TAM',
              'SPP1_TAM', 'Endothelial_cells', 'DCs']

fig, axes = plt.subplots(3, 4, figsize=(16, 12))
axes = axes.flatten()

for idx, ct in enumerate(cell_types):
    if ct in adata.obs.columns:
        sq.pl.spatial_scatter(
            adata,
            color=ct,
            title=ct.replace('_', ' '),
            spot_size=150,
            cmap='viridis',
            vmin=0,
            vmax=1,
            ax=axes[idx]
        )

plt.tight_layout()
plt.savefig("figures/celltype_proportions.png", dpi=300, bbox_inches='tight')
plt.close()
```

#### 3. Gene Set Expression In Situ (Fig. 1e from paper)

```python
# Define functional gene sets (from paper methods)
gene_sets = {
    'CAFs': ['COL1A1', 'COL1A2', 'COL3A1', 'FAP', 'PDGFRA'],
    'Tumor_Proliferation': ['MKI67', 'TOP2A', 'PCNA', 'CCNB1'],
    'Hypoxia': ['VEGFA', 'HIF1A', 'LDHA', 'PDK1', 'SLC2A1'],
    'EMT': ['VIM', 'TWIST1', 'SNAI1', 'SNAI2', 'ZEB1'],
    'Immune_Suppression': ['FOXP3', 'CTLA4', 'PDCD1', 'LAG3'],
    'Angiogenesis': ['VEGFA', 'ANGPT2', 'FLT1', 'KDR']
}

# Calculate signature scores
for signature, genes in gene_sets.items():
    # Filter genes present in dataset
    genes_present = [g for g in genes if g in adata.var_names]

    if len(genes_present) > 0:
        sc.tl.score_genes(
            adata,
            genes_present,
            score_name=f'{signature}_score',
            use_raw=False
        )

# Plot gene set scores
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for idx, (signature, genes) in enumerate(gene_sets.items()):
    if f'{signature}_score' in adata.obs.columns:
        sq.pl.spatial_scatter(
            adata,
            color=f'{signature}_score',
            title=signature.replace('_', ' '),
            spot_size=150,
            cmap='RdYlBu_r',
            ax=axes[idx]
        )

plt.tight_layout()
plt.savefig("figures/geneset_scores_spatial.png", dpi=300, bbox_inches='tight')
plt.close()
```

#### 4. Combined Visualization

```python
# Create a comprehensive figure matching paper style
fig = plt.figure(figsize=(20, 12))

# Subplot 1: H&E image
ax1 = plt.subplot(2, 4, 1)
if 'hires' in adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['images']:
    img = adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['images']['hires']
    ax1.imshow(img)
    ax1.set_title('H&E Image')
    ax1.axis('off')

# Subplots 2-4: Key markers
markers_to_plot = [
    ('COL1A1', 'Fibroblasts'),
    ('KRT8', 'Tumor'),
    ('CD79A', 'B cells')
]

for idx, (gene, label) in enumerate(markers_to_plot):
    ax = plt.subplot(2, 4, idx + 2)
    sq.pl.spatial_scatter(
        adata,
        color=gene,
        title=f'{gene} ({label})',
        spot_size=150,
        ax=ax,
        show=False
    )

# Subplots 5-8: Deconvolution results
cell_types_key = ['myCAF', 'Tumor_epithelial_cells', 'B_cells', 'CAFs']

for idx, ct in enumerate(cell_types_key):
    ax = plt.subplot(2, 4, idx + 5)
    if ct in adata.obs.columns:
        sq.pl.spatial_scatter(
            adata,
            color=ct,
            title=ct.replace('_', ' '),
            spot_size=150,
            ax=ax,
            show=False
        )

plt.tight_layout()
plt.savefig("figures/deconvolution_overview.png", dpi=300, bbox_inches='tight')
plt.close()
```

---

## Common Issues and Solutions

### Issue 1: Reference Data Format

**Problem**: Downloaded .rds file, but need h5ad for Python

**Solution**:
```r
# In R: Convert RDS to CSV
library(Seurat)
reference <- readRDS("scRNA-seq_Data_post_qc.rds")

# Export counts
write.csv(
  as.matrix(GetAssayData(reference, slot = "counts")),
  "reference_counts.csv"
)

# Export metadata
write.csv(
  reference@meta.data,
  "reference_metadata.csv"
)
```

```python
# In Python: Reconstruct h5ad
import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix

# Read CSV files
counts = pd.read_csv("reference_counts.csv", index_col=0)
metadata = pd.read_csv("reference_metadata.csv", index_col=0)

# Create AnnData
adata_ref = sc.AnnData(
    X=csr_matrix(counts.T.values),
    obs=metadata,
    var=pd.DataFrame(index=counts.index)
)

# Save
adata_ref.write("reference_scrna.h5ad")
```

### Issue 2: Gene Name Mismatch

**Problem**: Reference and spatial data have different gene names

**Solution**:
```python
# Find common genes
common_genes = list(set(adata_ref.var_names) & set(adata_spatial.var_names))
print(f"Common genes: {len(common_genes)} / {adata_spatial.n_vars}")

# Subset both datasets
adata_ref = adata_ref[:, common_genes]
adata_spatial = adata_spatial[:, common_genes]
```

### Issue 3: Different Platforms

**Problem**: scRNA-seq (10x) vs Spatial (Visium) have different characteristics

**Solution**:
- RCTD accounts for this automatically
- For other methods, normalize carefully:

```python
# Normalize both datasets similarly
sc.pp.normalize_total(adata_ref, target_sum=1e4)
sc.pp.log1p(adata_ref)

sc.pp.normalize_total(adata_spatial, target_sum=1e4)
sc.pp.log1p(adata_spatial)
```

### Issue 4: Memory Issues

**Problem**: Large datasets cause out-of-memory errors

**Solutions**:
```python
# 1. Process samples individually
sample_ids = adata.obs['sample_id'].unique()

for sample in sample_ids:
    adata_sample = adata[adata.obs['sample_id'] == sample].copy()
    # Run deconvolution
    # ...

# 2. Use subset of reference
# Take top N cells per cell type
from sklearn.utils import resample

ref_subset = []
for ct in adata_ref.obs['cell_type'].unique():
    cells = adata_ref[adata_ref.obs['cell_type'] == ct]
    if cells.n_obs > 500:
        cells = resample(cells, n_samples=500, random_state=42)
    ref_subset.append(cells)

adata_ref = sc.concat(ref_subset)
```

### Issue 5: Low-Quality Results

**Problem**: Deconvolution gives unrealistic proportions

**Diagnosis**:
```python
# Check proportion distributions
import seaborn as sns

cell_types = [col for col in adata.obs.columns if '_proportion' in col]

fig, axes = plt.subplots(1, len(cell_types), figsize=(20, 4))

for idx, ct in enumerate(cell_types):
    sns.histplot(adata.obs[ct], ax=axes[idx], bins=50)
    axes[idx].set_title(ct)
    axes[idx].set_xlabel('Proportion')

plt.tight_layout()
plt.savefig("figures/proportion_distributions.png")
```

**Solutions**:
1. **Check reference quality**: Ensure cell-type annotations are accurate
2. **Filter low-quality spots**: Remove spots with very low UMI counts
3. **Try different methods**: Compare RCTD, Cell2location, and others
4. **Adjust parameters**: Try different modes or settings

---

## Complete Pipeline Example

Here's a complete end-to-end example:

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import squidpy as sq
from pathlib import Path

# Configuration
results_dir = Path("results")
figures_dir = Path("figures")
figures_dir.mkdir(exist_ok=True, parents=True)

print("Step 1: Load processed spatial data")
adata_spatial = sc.read_h5ad(results_dir / "processed_data.h5ad")

print("Step 2: Load reference scRNA-seq data")
# Assume you've converted the RDS file to h5ad
adata_ref = sc.read_h5ad("data/reference_scrna.h5ad")

print("Step 3: Align gene names")
common_genes = list(set(adata_ref.var_names) & set(adata_spatial.var_names))
print(f"  Common genes: {len(common_genes)}")

adata_ref = adata_ref[:, common_genes].copy()
adata_spatial = adata_spatial[:, common_genes].copy()

print("Step 4: Run deconvolution")
# Here you would use RCTD in R, or Cell2location/Tangram in Python
# For this example, assume we've run RCTD and imported results

# Import RCTD results (saved from R)
rctd_results = pd.read_csv("results/rctd_proportions.csv", index_col=0)

# Add to spatial data
for col in rctd_results.columns:
    adata_spatial.obs[col] = rctd_results[col]

print("Step 5: Calculate gene set scores")
gene_sets = {
    'CAF_signature': ['COL1A1', 'COL1A2', 'COL3A1', 'FAP'],
    'Tumor_prolif': ['MKI67', 'TOP2A', 'PCNA'],
    'Immune_supp': ['FOXP3', 'CTLA4', 'PDCD1']
}

for signature, genes in gene_sets.items():
    genes_present = [g for g in genes if g in adata_spatial.var_names]
    if genes_present:
        sc.tl.score_genes(
            adata_spatial,
            genes_present,
            score_name=signature,
            use_raw=False
        )

print("Step 6: Visualization")
# Cell-type proportions
cell_types = [col for col in adata_spatial.obs.columns
              if any(x in col for x in ['CAF', 'Tumor', 'immune', 'cells'])]

for ct in cell_types[:4]:  # Plot first 4
    fig, ax = plt.subplots(figsize=(8, 6))
    sq.pl.spatial_scatter(
        adata_spatial,
        color=ct,
        title=ct.replace('_', ' '),
        spot_size=150,
        ax=ax,
        show=False
    )
    plt.savefig(figures_dir / f"deconv_{ct}.png", dpi=300, bbox_inches='tight')
    plt.close()

print("Step 7: Save results")
adata_spatial.write(results_dir / "adata_deconvolved.h5ad")

print("Deconvolution complete!")
print(f"Results saved to: {results_dir / 'adata_deconvolved.h5ad'}")
print(f"Figures saved to: {figures_dir}")
```

---

## Validation Checklist

Before finalizing your deconvolution results, verify:

- [ ] Cell-type proportions sum to ~1 for each spot
- [ ] Proportions are between 0 and 1
- [ ] High marker gene expression correlates with high cell-type proportion
- [ ] Spatial patterns make biological sense (e.g., immune cells at tumor borders)
- [ ] No spots with all zero proportions
- [ ] Distribution of proportions looks reasonable (not all 0 or 1)

---

## For Your Assignment

### Minimum Requirements (Task 2)

**Task 2a: Deconvolution**
1. Use RCTD (paper method) OR another algorithm (Cell2location, Tangram, etc.)
2. Document which method you used and why
3. Report:
   - Number of cell types identified
   - Proportion ranges for each cell type
   - Quality metrics

**Task 2b: Visualization**
1. Plot marker genes spatially (like Fig. 1c)
2. Plot deconvolution results (like Fig. 1d)
3. Plot gene set scores in situ (like Fig. 1e)
4. Compare marker expression with deconvolution results

### Recommended Workflow

1. **Start simple**: Use one sample first to test your pipeline
2. **Validate**: Check if results match known biology
3. **Scale up**: Apply to all samples once validated
4. **Compare**: If time permits, try multiple methods and compare

---

## Additional Resources

### Papers
1. **RCTD**: Cable et al., Nature Biotechnology 2022
2. **Cell2location**: Kleshchevnikov et al., Nature Biotechnology 2022
3. **Tangram**: Biancalani et al., Nature Methods 2021
4. **Review**: Dries et al., Nature Methods 2021

### Tutorials
1. RCTD: https://raw.githack.com/dmcable/spacexr/master/vignettes/spatial-transcriptomics.html
2. Cell2location: https://cell2location.readthedocs.io/
3. Squidpy visualization: https://squidpy.readthedocs.io/

### Tips for the Assignment
- **Time management**: Deconvolution can take hours - start early!
- **Computational resources**: Use a subset of data for testing
- **Language choice**: R for RCTD (exact replication), Python for alternatives
- **Documentation**: Keep notes on parameters and decisions

---

**Good luck with your deconvolution analysis!**
