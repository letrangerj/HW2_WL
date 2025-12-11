# Tutorial: Approaching the Spatial Transcriptomics Assignment

## Overview
This assignment analyzes PDAC spatial transcriptomics data using Python (Scanpy). The workflow follows the paper's methodology but can be adapted from your provided scripts.

---

## Task 0: Environment Setup (20 points)

**What you need:**
- Miniconda/Conda environment
- Key packages: `scanpy`, `pandas`, `numpy`, `matplotlib`, `seaborn`, `PIL`
- For deconvolution: Consider `cell2location` or implement RCTD-like approach
- For spatial analysis: `squidpy` (spatial statistics)

**Approach from your scripts:**
```python
# Your hd.py shows the basic imports needed
import scanpy as sc
import pandas as pd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import seaborn as sns
```

---

## Task 1: Preprocessing and Quality Control (30 points)

**Key steps from hd.py:**

### 1.1 Load Data

**Understanding the Data Structure:**

The dataset is organized in `data/GSE272362/HM_files/` with **30 samples** (HM3332 to HM3361). Each sample directory contains:

```
HM3332/
├── expression.h5ad                           # Gene expression matrix (spots × genes)
└── spatial/
    ├── GSM8399450_*_tissue_positions.csv     # Spatial coordinates of spots
    ├── GSM8399450_*_scalefactors_json.json   # Scaling factors for image alignment
    └── GSM8399450_*_tissue_hires_image.png   # H&E histology image
```

**Key Components Explained:**

1. **expression.h5ad**: Pre-formatted AnnData object containing:
   - `adata.X`: Raw gene expression counts (spots × genes matrix)
   - `adata.obs`: Spot metadata (will be enriched with spatial info)
   - `adata.var`: Gene metadata

2. **tissue_positions.csv**: Contains spatial coordinates for each barcode:
   - `barcode`: Unique spot identifier (e.g., "ACGCCTGACACGCGCT-1")
   - `in_tissue`: Binary flag (1=on tissue, 0=background)
   - `array_row`, `array_col`: Array grid coordinates
   - `pxl_row_in_fullres`, `pxl_col_in_fullres`: Pixel coordinates in full-resolution image

3. **scalefactors_json.json**: Scaling factors for image-to-spot alignment:
   - `tissue_hires_scalef`: Scale factor for high-resolution image (e.g., 0.08126)
   - `spot_diameter_fullres`: Spot diameter in pixels (~176 pixels for Visium)

4. **tissue_hires_image.png**: H&E stained histology image (usually ~2000×2000 pixels)

**Data Loading Function:**

Here's a robust function to load a single sample:

```python
import scanpy as sc
import pandas as pd
import numpy as np
import json
from PIL import Image
from pathlib import Path

def load_single_sample(sample_id, base_dir="data/GSE272362/HM_files"):
    """
    Load a single spatial transcriptomics sample with all spatial information.

    Parameters:
    -----------
    sample_id : str
        Sample identifier (e.g., "HM3332")
    base_dir : str
        Base directory containing HM_files

    Returns:
    --------
    adata : AnnData
        AnnData object with spatial information in:
        - adata.obs: includes in_tissue, array coordinates, pixel coordinates
        - adata.obsm['spatial']: spatial coordinates for plotting
        - adata.uns['spatial'][sample_id]: images and scalefactors
    """
    data_dir = Path(base_dir) / sample_id

    # 1. Load h5ad file with gene expression data
    print(f"Loading {sample_id}...")
    adata = sc.read_h5ad(data_dir / "expression.h5ad")
    adata.var_names_make_unique()  # Ensure unique gene names

    # 2. Auto-detect spatial files (handles different GSM IDs per sample)
    spatial_dir = data_dir / "spatial"
    tissue_positions_path = list(spatial_dir.glob("*_tissue_positions.csv"))[0]
    scalefactors_path = list(spatial_dir.glob("*_scalefactors_json.json"))[0]
    image_path = list(spatial_dir.glob("*_tissue_hires_image.png"))[0]

    # 3. Load and merge spatial coordinates
    positions = pd.read_csv(tissue_positions_path)
    positions.set_index("barcode", inplace=True)
    # Merge spatial info into adata.obs (left join to keep all spots)
    adata.obs = adata.obs.merge(positions, left_index=True, right_index=True, how="left")

    # 4. Load scalefactors for image alignment
    with open(scalefactors_path, "r") as f:
        scalefactors = json.load(f)

    # 5. Load H&E histology image
    image = np.array(Image.open(image_path))

    # 6. Store spatial information in standard Scanpy format
    adata.uns["spatial"] = {
        sample_id: {
            "images": {"hires": image},
            "scalefactors": scalefactors
        }
    }

    # 7. Add spatial coordinates to obsm (required for sc.pl.spatial)
    # Use pixel coordinates for accurate spatial plotting
    adata.obsm["spatial"] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].to_numpy()

    # 8. Add sample ID for tracking when combining multiple samples
    adata.obs["sample_id"] = sample_id

    print(f"  Loaded: {adata.n_obs} spots, {adata.n_vars} genes")
    print(f"  In-tissue spots: {adata.obs['in_tissue'].sum()}")
    print(f"  Image shape: {image.shape}")

    return adata
```

**Usage for Single Sample:**

```python
# Load one sample
adata = load_single_sample("HM3332")

# Verify the structure
print(adata)
print(adata.obs.columns)  # Should include spatial coordinates
print(adata.obsm.keys())  # Should include 'spatial'
print(adata.uns['spatial'].keys())  # Should include 'HM3332'

# Optional: Filter to in-tissue spots only
adata = adata[adata.obs["in_tissue"] == 1].copy()

# Quick visualization
sc.pl.spatial(adata, img_key="hires", color="in_tissue", spot_size=1.5)
```

**Usage for Multiple Samples (Task 3):**

```python
# Load multiple samples
sample_ids = ["HM3332", "HM3333", "HM3334"]
adatas = [load_single_sample(sid) for sid in sample_ids]

# Option A: Concatenate for batch analysis
adata_combined = sc.concat(adatas, label="batch", keys=sample_ids)

# Option B: Process separately and compare
for adata in adatas:
    # Perform analysis on each sample
    pass
```

**Why Use `glob()` Pattern Matching?**

Each sample has a different GSM ID in the filename (e.g., GSM8399450, GSM8399451, etc.). Using `glob()` auto-detects the correct files without hardcoding IDs, making the code more robust and reusable.

**Common Issues and Solutions:**

1. **Missing spatial coordinates**: Ensure tissue_positions.csv barcodes match adata.obs_names
2. **Image not displaying**: Check that scalefactors and image are in adata.uns['spatial'][sample_id]
3. **Index mismatch**: Use `how="left"` in merge to keep all spots from adata

Now you're ready to proceed with QC!

### 1.2 QC Metrics (lines 82-100 in hd.py)
```python
# Calculate QC metrics
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# Visualize distributions
- total_counts histogram
- n_genes_by_counts histogram  
- pct_counts_mt
```

### 1.3 Filtering
**Key differences from scRNA-seq:**
- Spots contain multiple cells (not single cells)
- Higher count thresholds (e.g., >200 genes, not >500)
- Tissue-specific background may be higher
- Spatial patterns matter (don't just filter by counts)

**From hd.py approach:**
```python
sc.pp.filter_cells(adata, min_counts=100)  # Adjust threshold
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()
```

**What to visualize:**
- QC metrics spatially overlaid on H&E (using `sc.pl.spatial()`)
- Show which spots get filtered and verify they're low-quality regions

---

## Task 2: Cell Type Annotation via Deconvolution (30 points)

**Paper uses RCTD** - you'll need to either:
1. Use the reference scRNA-seq data provided (`.rds` file)
2. Implement deconvolution (cell2location, Stereoscope, or similar)

**Workflow:**

### 2.1 Prepare Reference
- Load scRNA-seq reference (convert .rds to h5ad if needed)
- Ensure cell type annotations are clear
- Match gene names between reference and spatial data

### 2.2 Run Deconvolution
Since RCTD is R-based, consider:
- `cell2location` (Python, similar approach)
- `tangram` (Python, maps cells to space)
- Or call R from Python for RCTD

### 2.3 Visualize Results
**Key visualizations (from paper Fig 1c,d,e):**
```python
# Cell type proportions per spot
sc.pl.spatial(adata, color=['Fibroblast_prop', 'Tumor_cell_prop', 
                             'B_cell_prop', 'T_cell_prop'])

# Marker gene expression
sc.pl.spatial(adata, color=['COL1A1',  # Fibroblasts
                             'KRT8',    # Epithelial/Tumor
                             'EPCAM',   # Epithelial
                             'CD79A'])  # B cells
```

**Gene set expression:**
- Select pathways/signatures (ECM, hypoxia, proliferation)
- Calculate scores: `sc.tl.score_genes()`
- Visualize spatially

---

## Task 3: Cell-Cell Interaction Analysis (20 points)

**Paper uses MISTy** - alternative Python approaches:

### 3.1 Spatial Neighbor Analysis (from spatial_cluster.py concept)

```python
import squidpy as sq

# Build spatial graph
sq.gr.spatial_neighbors(adata, coord_type="generic", 
                        spatial_key="spatial")

# Co-occurrence analysis
sq.gr.co_occurrence(adata, cluster_key="cell_type_dominant")

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key="cell_type_dominant")

# Visualize interaction matrix
sq.pl.nhood_enrichment(adata, cluster_key="cell_type_dominant")
```

### 3.2 Reproduce Paper's Fig 2a
**Target patterns:**
- Tumor-CAF interactions
- Immune cell clustering
- TAM-T cell proximity
- Spatial ecotypes (CC1-CC10)

### 3.3 Histology Correlation
- Compare interaction patterns to H&E morphology
- Identify tumor borders vs. interior
- Stromal-rich vs. tumor-rich regions

**From spatial_cluster.py approach:**
```python
# Spatial clustering based on location
sc.pp.neighbors(adata, n_neighbors=8, use_rep='spatial', 
                metric='euclidean')
sc.tl.leiden(adata, resolution=0.3)
```

---

## General Workflow Structure

**Recommended script organization:**

```python
# 1. Setup and Data Loading
- Import libraries
- Load h5ad and spatial info
- Create adata object

# 2. QC and Preprocessing  
- Calculate metrics
- Filter spots/genes
- Normalize and log-transform
- Identify highly variable genes

# 3. Deconvolution
- Load reference
- Run deconvolution algorithm
- Add results to adata.obs

# 4. Spatial Analysis
- Calculate cell type proportions
- Spatial clustering
- Cell-cell interactions
- Pathway enrichment

# 5. Visualization
- QC plots
- Spatial cell type maps
- Marker expression
- Interaction networks
```

---

## Key Considerations

### Multiple Sections
- The assignment notes you can use one section OR multiple
- For Task 3, multiple sections give better statistics
- Loop through samples if analyzing multiple

### Differences from Paper
- You don't need exact replication
- Parameter sensitivity is expected
- Focus on major trends (desmoplastic vs. proliferative ecotypes)

### Figure Quality
- Use `dpi=300` for publication quality
- `bbox_inches='tight'` to avoid cropping
- Save as PNG or PDF

### Computing Resources
- Deconvolution is memory-intensive
- May need HPC/cluster for multiple samples
- Consider subsetting data for testing

---

## Report Structure Suggestion

1. **Introduction** (brief background on PDAC spatial biology)
2. **Methods** (QC criteria, tools used, parameters)
3. **Results:**
   - QC and filtering rationale
   - Cell type distributions
   - Spatial patterns and interactions
4. **Discussion** (compare to paper, biological insights)
5. **Figures** (well-labeled, high-quality)

---

## Useful References from Your Scripts

**From hd.py:**
- Lines 82-100: QC metrics calculation
- Lines 103-110: Filtering approach
- Lines 113-123: Standard normalization/PCA
- Lines 137-143: Clustering
- Lines 172-175: Spatial visualization

**From spatial_cluster.py:**
- Complete CLI structure if you want arguments
- Spatial clustering function (lines 94-112)
- Modular function design

---

## Tips

1. **Start small:** Test on one section first
2. **Check intermediate outputs:** Save QC plots before filtering
3. **Validate deconvolution:** Compare with known markers
4. **Iterate on thresholds:** QC cutoffs affect downstream results
5. **Document choices:** Why you filtered at certain thresholds

Good luck! The key is understanding the biological context from the paper while adapting the computational workflow from your scripts.
