# Quality Control (QC) Tutorial for Spatial Transcriptomics

## Table of Contents
1. [Overview](#overview)
2. [Understanding QC Metrics](#understanding-qc-metrics)
3. [FFPE vs Fresh Frozen Samples](#ffpe-vs-fresh-frozen-samples)
4. [How to Determine Thresholds](#how-to-determine-thresholds)
5. [Code Examples](#code-examples)
6. [Common Pitfalls](#common-pitfalls)

---

## Overview

Quality control is a critical first step in spatial transcriptomics analysis. Unlike single-cell RNA-seq where you filter individual cells, in Visium data you filter **spots** (each containing multiple cells).

**Key principle**: Balance removing low-quality data while preserving biological information and tissue architecture.

---

## Understanding QC Metrics

### 1. Total Counts (UMI counts)
- **What**: Total number of transcripts detected in a spot
- **Interpretation**:
  - Low counts → Poor RNA capture, low cell density, or degraded RNA
  - Very high counts → Potential doublets or highly transcriptionally active regions
- **Typical range (Visium)**: 500-10,000 for fresh frozen, 200-5,000 for FFPE

### 2. Number of Genes Detected (n_genes_by_counts)
- **What**: Count of unique genes with at least 1 UMI in a spot
- **Interpretation**:
  - Low gene count + high total counts → Technical artifact (same gene amplified many times)
  - Low gene count + low total counts → Poor quality spot
- **Typical range**: 200-3,000 genes per spot

### 3. Mitochondrial Percentage (pct_counts_mt)
- **What**: Percentage of reads mapping to mitochondrial genes
- **Interpretation**:
  - High MT% (>20%) → Dying/stressed cells, cytoplasmic RNA leaked out
  - **FFPE samples**: Often 0% because MT genes are degraded during fixation
- **Typical threshold**: <10-25% for fresh frozen, often N/A for FFPE

### 4. In-Tissue Status
- **What**: Whether the spot overlaps with tissue (from histology image)
- **Use**: Pre-filter to exclude spots clearly outside tissue boundaries

---

## FFPE vs Fresh Frozen Samples

### Fresh Frozen (Standard Protocol)

**Characteristics**:
- High RNA quality
- MT genes detectable
- Higher UMI counts
- More consistent quality across spots

**Typical thresholds**:
```python
sc.pp.filter_cells(adata, min_counts=500-1000)
adata = adata[adata.obs['n_genes_by_counts'] >= 200-500, :]
adata = adata[adata.obs['pct_counts_mt'] < 10-20, :]
```

### FFPE (Our Dataset)

**Characteristics**:
- RNA degradation from fixation process
- MT genes often undetectable (MT% = 0%)
- Lower overall UMI counts
- More heterogeneous quality across tissue sections
- Spatial context is critical (degradation may be biological, e.g., necrotic regions)

**Typical thresholds** (more lenient):
```python
sc.pp.filter_cells(adata, min_counts=300-500)
adata = adata[adata.obs['n_genes_by_counts'] >= 150-300, :]
# No MT% filter (unreliable)
```

**Key differences**:
1. **Lower thresholds**: FFPE requires 30-50% lower thresholds than fresh frozen
2. **No MT% filtering**: MT genes are degraded, so this metric is unreliable
3. **Spatial visualization essential**: Must check if low-quality spots represent biological regions or technical artifacts

---

## How to Determine Thresholds

### Step 1: Calculate QC Metrics

```python
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

# Calculate metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Examine distributions
print("=== QC Metrics Summary ===")
print(f"Total counts - Median: {adata.obs['total_counts'].median():.0f}, "
      f"Range: {adata.obs['total_counts'].min():.0f}-{adata.obs['total_counts'].max():.0f}")
print(f"Genes detected - Median: {adata.obs['n_genes_by_counts'].median():.0f}, "
      f"Range: {adata.obs['n_genes_by_counts'].min():.0f}-{adata.obs['n_genes_by_counts'].max():.0f}")
print(f"MT% - Median: {adata.obs['pct_counts_mt'].median():.2f}%, "
      f"Range: {adata.obs['pct_counts_mt'].min():.2f}-{adata.obs['pct_counts_mt'].max():.2f}%")
```

### Step 2: Visualize Distributions (Histograms)

**What to look for**:

```
Normal Distribution:              Bimodal Distribution:
     |                                 |
  *  |  *                           * | *      *
 *   | *   *                        * |*  *   *
*    |*    *                       *  |  *  *
_____|______                       ____|____*____
     ^                                 ^    ^
  threshold                         noise  signal
```

**Analysis guide**:
- Identify the main peak (represents majority of spots)
- Look for long tails on the left (low-quality spots)
- Set threshold to cut the tail without touching the main peak
- Typically remove bottom 5-15% of spots

```python
# Plot histograms
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Total counts
axes[0].hist(adata.obs['total_counts'], bins=50, alpha=0.7, edgecolor='black')
axes[0].axvline(adata.obs['total_counts'].median(), color='green',
                linestyle='--', label='Median')
axes[0].axvline(np.percentile(adata.obs['total_counts'], 10),
                color='red', linestyle='--', label='10th percentile')
axes[0].set_xlabel('Total counts')
axes[0].set_ylabel('Frequency')
axes[0].legend()
axes[0].set_title('Total Counts Distribution')

# Genes detected
axes[1].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7, edgecolor='black')
axes[1].axvline(adata.obs['n_genes_by_counts'].median(), color='green',
                linestyle='--', label='Median')
axes[1].axvline(np.percentile(adata.obs['n_genes_by_counts'], 10),
                color='red', linestyle='--', label='10th percentile')
axes[1].set_xlabel('Number of genes')
axes[1].legend()
axes[1].set_title('Genes Detected Distribution')

# Scatter plot: counts vs genes
axes[2].scatter(adata.obs['total_counts'], adata.obs['n_genes_by_counts'],
                alpha=0.3, s=2)
axes[2].set_xlabel('Total counts')
axes[2].set_ylabel('Number of genes')
axes[2].set_title('Counts vs Genes')

plt.tight_layout()
plt.savefig('qc_distributions.png', dpi=300)
```

**Percentile-based approach**:
```python
# Calculate percentiles
print("\n=== Percentile Analysis ===")
for p in [5, 10, 25]:
    print(f"{p}th percentile - Counts: {np.percentile(adata.obs['total_counts'], p):.0f}, "
          f"Genes: {np.percentile(adata.obs['n_genes_by_counts'], p):.0f}")
```

### Step 3: Visualize Spatial Patterns (CRITICAL!)

**This is unique to spatial data and absolutely essential.**

```python
# Spatial visualization of QC metrics
sc.pl.spatial(adata, color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'],
              spot_size=150, ncols=3)
```

**Spatial Pattern Interpretation**:

**Pattern A - Edge Artifacts (SAFE to filter)**:
```
Tissue cross-section:
┌─────────────────┐
│  High quality   │  Interior: High counts
│    ████████     │
│    ████████     │
│   Low ████      │  Edge: Low counts (technical artifact)
└─────────────────┘
```
- Low counts concentrated at tissue edges
- Result of tissue damage during sectioning
- **Action**: Filter these spots

**Pattern B - Biological Regions (DO NOT filter)**:
```
Tissue cross-section:
┌─────────────────┐
│ Tumor (High)    │  Different regions with
│   ████████      │  different biology
│ Stroma (Low)    │
│   ░░░░░░░░      │
└─────────────────┘
```
- Low counts in specific histological regions (e.g., stroma, fibrotic areas)
- Represents real biological differences
- **Action**: Keep these spots! They contain important information

**Pattern C - Necrotic/Damaged Areas (Contextual)**:
```
Tissue cross-section:
┌─────────────────┐
│  ████    ░░░    │  Necrotic region
│  ████    ░░░    │  (very low counts)
│  ████████       │
└─────────────────┘
```
- Obvious tissue damage or necrosis on H&E image
- **Action**: Decide based on research question
  - Keep if studying tumor necrosis
  - Filter if necrosis is not relevant

### Step 4: Test Thresholds Iteratively

**Create a function to test different thresholds**:

```python
def test_qc_threshold(adata, min_counts, min_genes):
    """
    Test impact of QC thresholds without modifying original data
    """
    # Create pass/fail annotation
    pass_qc = (
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['n_genes_by_counts'] >= min_genes)
    )

    # Statistics
    n_pass = pass_qc.sum()
    n_fail = (~pass_qc).sum()

    print(f"\n=== Threshold Test: {min_counts} counts, {min_genes} genes ===")
    print(f"Pass: {n_pass} spots ({n_pass/adata.n_obs*100:.1f}%)")
    print(f"Fail: {n_fail} spots ({n_fail/adata.n_obs*100:.1f}%)")

    print(f"\nFiltered spots (mean ± std):")
    print(f"  Counts: {adata.obs[~pass_qc]['total_counts'].mean():.0f} ± "
          f"{adata.obs[~pass_qc]['total_counts'].std():.0f}")
    print(f"  Genes: {adata.obs[~pass_qc]['n_genes_by_counts'].mean():.0f} ± "
          f"{adata.obs[~pass_qc]['n_genes_by_counts'].std():.0f}")

    print(f"\nRetained spots (mean ± std):")
    print(f"  Counts: {adata.obs[pass_qc]['total_counts'].mean():.0f} ± "
          f"{adata.obs[pass_qc]['total_counts'].std():.0f}")
    print(f"  Genes: {adata.obs[pass_qc]['n_genes_by_counts'].mean():.0f} ± "
          f"{adata.obs[pass_qc]['n_genes_by_counts'].std():.0f}")

    # Visualize
    adata.obs['qc_status'] = 'Pass'
    adata.obs.loc[~pass_qc, 'qc_status'] = 'Fail'

    sc.pl.spatial(adata, color='qc_status', spot_size=150,
                  title=f'QC: {min_counts} counts, {min_genes} genes')

    return pass_qc

# Test different thresholds
test_qc_threshold(adata, min_counts=200, min_genes=100)  # Lenient
test_qc_threshold(adata, min_counts=300, min_genes=150)  # Moderate
test_qc_threshold(adata, min_counts=500, min_genes=200)  # Conservative
```

### Step 5: Make Decision

**Decision criteria**:

1. **Statistical**: Threshold should be in the tail (5-15th percentile), not cutting main distribution
2. **Spatial**: Filtered spots should be edges/artifacts, not entire biological regions
3. **Percentage**: Typically filter 5-20% of spots
   - <5%: Too lenient, keeping noise
   - >30%: Too strict, losing biology
4. **Downstream**: Consider your analysis goals
   - Deconvolution: Can tolerate lower quality
   - Gene expression analysis: Need higher quality

---

## Code Examples

### Complete QC Pipeline for FFPE Visium

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def perform_qc(adata, min_counts=500, min_genes=200, output_dir='./figures'):
    """
    Complete QC pipeline for FFPE Visium data

    Parameters:
    -----------
    adata : AnnData
        Raw spatial transcriptomics data
    min_counts : int
        Minimum total counts per spot
    min_genes : int
        Minimum genes detected per spot
    output_dir : str
        Directory to save QC plots

    Returns:
    --------
    adata_filtered : AnnData
        Filtered data passing QC
    """

    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Step 1: Calculate QC metrics
    print("Calculating QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Step 2: Visualize BEFORE filtering
    print("Creating pre-QC visualizations...")

    # Histograms
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    axes[0].hist(adata.obs['total_counts'], bins=50, alpha=0.7, edgecolor='black')
    axes[0].axvline(min_counts, color='red', linestyle='--',
                    label=f'Threshold: {min_counts}')
    axes[0].axvline(adata.obs['total_counts'].median(), color='green',
                    linestyle=':', label=f'Median: {adata.obs["total_counts"].median():.0f}')
    axes[0].set_xlabel('Total counts')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Total Counts Distribution')
    axes[0].legend()

    axes[1].hist(adata.obs['n_genes_by_counts'], bins=50, alpha=0.7, edgecolor='black')
    axes[1].axvline(min_genes, color='red', linestyle='--',
                    label=f'Threshold: {min_genes}')
    axes[1].axvline(adata.obs['n_genes_by_counts'].median(), color='green',
                    linestyle=':', label=f'Median: {adata.obs["n_genes_by_counts"].median():.0f}')
    axes[1].set_xlabel('Number of genes')
    axes[1].set_title('Genes Detected Distribution')
    axes[1].legend()

    # Scatter plot with QC thresholds
    pass_qc = (
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['n_genes_by_counts'] >= min_genes)
    )
    axes[2].scatter(adata.obs.loc[pass_qc, 'total_counts'],
                    adata.obs.loc[pass_qc, 'n_genes_by_counts'],
                    alpha=0.5, s=2, label='Pass QC', color='blue')
    axes[2].scatter(adata.obs.loc[~pass_qc, 'total_counts'],
                    adata.obs.loc[~pass_qc, 'n_genes_by_counts'],
                    alpha=0.5, s=2, label='Fail QC', color='red')
    axes[2].axvline(min_counts, color='red', linestyle='--', linewidth=1)
    axes[2].axhline(min_genes, color='red', linestyle='--', linewidth=1)
    axes[2].set_xlabel('Total counts')
    axes[2].set_ylabel('Number of genes')
    axes[2].set_title('Counts vs Genes')
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(f"{output_dir}/qc_distributions_before.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Spatial plots
    adata.obs['qc_status'] = 'Pass'
    adata.obs.loc[~pass_qc, 'qc_status'] = 'Fail'

    sc.pl.spatial(adata, color=['total_counts', 'n_genes_by_counts', 'qc_status'],
                  spot_size=150, ncols=3, show=False)
    plt.savefig(f"{output_dir}/qc_spatial_before.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Step 3: Print statistics
    print(f"\n=== QC Statistics (Thresholds: {min_counts} counts, {min_genes} genes) ===")
    print(f"Total spots: {adata.n_obs}")
    print(f"Pass QC: {pass_qc.sum()} ({pass_qc.sum()/adata.n_obs*100:.1f}%)")
    print(f"Fail QC: {(~pass_qc).sum()} ({(~pass_qc).sum()/adata.n_obs*100:.1f}%)")

    print(f"\nFiltered spots - Mean ± SD:")
    print(f"  Total counts: {adata.obs[~pass_qc]['total_counts'].mean():.0f} ± "
          f"{adata.obs[~pass_qc]['total_counts'].std():.0f}")
    print(f"  Genes detected: {adata.obs[~pass_qc]['n_genes_by_counts'].mean():.0f} ± "
          f"{adata.obs[~pass_qc]['n_genes_by_counts'].std():.0f}")

    print(f"\nRetained spots - Mean ± SD:")
    print(f"  Total counts: {adata.obs[pass_qc]['total_counts'].mean():.0f} ± "
          f"{adata.obs[pass_qc]['total_counts'].std():.0f}")
    print(f"  Genes detected: {adata.obs[pass_qc]['n_genes_by_counts'].mean():.0f} ± "
          f"{adata.obs[pass_qc]['n_genes_by_counts'].std():.0f}")

    # Step 4: Filter genes and spots
    print("\nApplying filters...")
    sc.pp.filter_genes(adata, min_counts=3)  # Remove genes with <3 total counts
    sc.pp.filter_cells(adata, min_counts=min_counts)
    adata = adata[adata.obs['n_genes_by_counts'] >= min_genes, :].copy()

    print(f"After filtering: {adata.n_obs} spots, {adata.n_vars} genes")

    # Step 5: Normalization and log transformation
    print("Normalizing and log-transforming...")
    sc.pp.normalize_total(adata, target_sum=1e6)  # CPM normalization
    sc.pp.log1p(adata)  # log(x+1) transformation

    # Step 6: PCA for dimensionality reduction
    print("Running PCA...")
    sc.pp.pca(adata, n_comps=50)

    # Step 7: Visualize AFTER filtering
    print("Creating post-QC visualizations...")

    sc.pl.spatial(adata, color=['total_counts', 'n_genes_by_counts'],
                  spot_size=150, ncols=2, show=False)
    plt.savefig(f"{output_dir}/qc_spatial_after.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("QC complete!")

    return adata

# Example usage
adata_filtered = perform_qc(
    adata,
    min_counts=500,  # FFPE-appropriate threshold
    min_genes=200,   # FFPE-appropriate threshold
    output_dir='./figures'
)

# Save filtered data
adata_filtered.write('results/adata_qc_filtered.h5ad')
```

### Comparing Multiple Thresholds

```python
def compare_thresholds(adata, threshold_sets, output_dir='./figures'):
    """
    Compare multiple threshold combinations

    Parameters:
    -----------
    adata : AnnData
        Raw data
    threshold_sets : list of dict
        Each dict should have 'min_counts', 'min_genes', and 'name' keys
    """

    fig, axes = plt.subplots(2, len(threshold_sets),
                             figsize=(5*len(threshold_sets), 10))

    for i, thresh in enumerate(threshold_sets):
        min_counts = thresh['min_counts']
        min_genes = thresh['min_genes']
        name = thresh['name']

        # Calculate pass/fail
        pass_qc = (
            (adata.obs['total_counts'] >= min_counts) &
            (adata.obs['n_genes_by_counts'] >= min_genes)
        )

        adata.obs['qc_temp'] = 'Pass'
        adata.obs.loc[~pass_qc, 'qc_temp'] = 'Fail'

        # Histogram
        axes[0, i].hist(adata.obs['total_counts'], bins=50, alpha=0.7)
        axes[0, i].axvline(min_counts, color='red', linestyle='--', linewidth=2)
        axes[0, i].set_title(f'{name}\n{pass_qc.sum()}/{adata.n_obs} spots '
                             f'({pass_qc.sum()/adata.n_obs*100:.1f}%)')
        axes[0, i].set_xlabel('Total counts')

        # Spatial (can't use sc.pl.spatial with subplots easily, so describe)
        # In practice, you'd save separate spatial plots
        axes[1, i].scatter(
            adata.obs.loc[pass_qc, 'array_col'],
            adata.obs.loc[pass_qc, 'array_row'],
            s=5, alpha=0.5, label='Pass', color='blue'
        )
        axes[1, i].scatter(
            adata.obs.loc[~pass_qc, 'array_col'],
            adata.obs.loc[~pass_qc, 'array_row'],
            s=5, alpha=0.5, label='Fail', color='red'
        )
        axes[1, i].set_title(f'{name} - Spatial Distribution')
        axes[1, i].set_xlabel('Array column')
        axes[1, i].set_ylabel('Array row')
        axes[1, i].legend()
        axes[1, i].invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{output_dir}/threshold_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()

# Example usage
threshold_sets = [
    {'min_counts': 200, 'min_genes': 100, 'name': 'Very Lenient'},
    {'min_counts': 300, 'min_genes': 150, 'name': 'Lenient'},
    {'min_counts': 500, 'min_genes': 200, 'name': 'Moderate'},
    {'min_counts': 700, 'min_genes': 300, 'name': 'Strict'}
]

compare_thresholds(adata, threshold_sets)
```

---

## Common Pitfalls

### 1. Using Single-Cell QC Thresholds on Spatial Data
**Problem**: Single-cell thresholds are too strict for spatial spots (which contain multiple cells)

**Wrong**:
```python
# DON'T use scRNA-seq thresholds!
sc.pp.filter_cells(adata, min_counts=1000)  # Too high for FFPE
adata = adata[adata.obs['n_genes_by_counts'] >= 500, :]  # Too high
```

**Right**:
```python
# Use spatial-appropriate thresholds
sc.pp.filter_cells(adata, min_counts=300-500)  # FFPE
adata = adata[adata.obs['n_genes_by_counts'] >= 150-300, :]
```

### 2. Ignoring Spatial Patterns
**Problem**: Filtering based on distributions alone without checking spatial patterns

**Wrong**:
```python
# Only looking at histograms, not spatial patterns
sc.pp.filter_cells(adata, min_counts=500)
# Might remove entire biological regions!
```

**Right**:
```python
# Always visualize spatially BEFORE deciding
sc.pl.spatial(adata, color='total_counts', spot_size=150)
# Then decide if low-count regions are technical or biological
```

### 3. Over-filtering FFPE Data
**Problem**: Applying fresh frozen thresholds to FFPE samples

**Impact**: Lose 50%+ of spots, destroy tissue architecture

**Solution**: Use FFPE-specific, more lenient thresholds

### 4. Using MT% for FFPE Samples
**Problem**: MT genes are degraded in FFPE, so MT% is always ~0%

**Wrong**:
```python
# This filter does nothing in FFPE!
adata = adata[adata.obs['pct_counts_mt'] < 20, :]
```

**Right**:
```python
# For FFPE, focus on counts and gene detection
# Skip MT% filtering
```

### 5. Not Documenting Threshold Choices
**Problem**: Using arbitrary thresholds without justification

**Wrong**: Just applying filters with no explanation

**Right**:
- Show distributions with thresholds marked
- Show spatial patterns of filtered spots
- Report percentages filtered
- Justify choices based on FFPE characteristics and biological context

### 6. Filtering Too Early
**Problem**: Filtering before examining the data

**Wrong**:
```python
# Filter immediately without looking
adata = sc.read_h5ad('data.h5ad')
sc.pp.filter_cells(adata, min_counts=500)
```

**Right**:
```python
# Always examine first
adata = sc.read_h5ad('data.h5ad')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Visualize
sc.pl.spatial(adata, color=['total_counts', 'n_genes_by_counts'])

# Then decide on thresholds
```

### 7. Not Filtering Genes
**Problem**: Only filtering spots, not low-count genes

**Wrong**:
```python
sc.pp.filter_cells(adata, min_counts=500)
# Forgot to filter genes!
```

**Right**:
```python
# Filter BOTH genes and spots
sc.pp.filter_genes(adata, min_counts=3)  # Remove genes with very few reads
sc.pp.filter_cells(adata, min_counts=500)  # Remove low-quality spots
```

---

## Summary Checklist

For your assignment (Task 1), make sure you:

- [ ] Calculate QC metrics using `sc.pp.calculate_qc_metrics()`
- [ ] Create histograms showing distributions with threshold lines
- [ ] Create spatial plots showing QC metrics overlaid on tissue
- [ ] Test multiple threshold combinations
- [ ] Report percentage of spots filtered
- [ ] Show spatial distribution of pass/fail spots
- [ ] Justify thresholds based on:
  - Distribution analysis (percentiles)
  - Spatial patterns (edge vs interior)
  - FFPE-specific considerations
  - Biological context (tumor heterogeneity)
- [ ] Document differences between FFPE and fresh frozen QC
- [ ] Show before/after filtering comparisons

---

## Recommended Thresholds for This Assignment

Based on the PDAC FFPE dataset characteristics:

**Conservative (Recommended)**:
```python
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=500)
adata = adata[adata.obs['n_genes_by_counts'] >= 200, :]
```

**Moderate**:
```python
sc.pp.filter_genes(adata, min_counts=3)
sc.pp.filter_cells(adata, min_counts=300)
adata = adata[adata.obs['n_genes_by_counts'] >= 150, :]
```

Choose based on your spatial visualization results!

---

## References

1. Scanpy tutorials: https://scanpy-tutorials.readthedocs.io/
2. Seurat spatial vignettes: https://satijalab.org/seurat/articles/spatial_vignette.html
3. Original paper: Nature Genetics (2024) - PDAC spatial atlas
4. FFPE spatial considerations: 10x Genomics technical notes

---

**Good luck with your analysis!**
