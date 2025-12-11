# Import core packages
import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
import matplotlib.pyplot as plt
from PIL import Image
import geopandas as gpd
from matplotlib.patches import Polygon
import seaborn as sns

# Disable annoying warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# set working directory
import os
os.chdir('/lustre1/share/hd/')
output_dir = Path("/home/bjx143_pkuhpc/spatial_WL/output/")
output_dir.mkdir(exist_ok=True, parents=True)

# save plot
def save_plot(fig, name, dpi=150):
    fig.savefig(output_dir / f"{name}.png", bbox_inches="tight", dpi=dpi)
    plt.close(fig)

# load data from spaceranger output
adata = sc.read_10x_h5('spaceranger_output/outs/binned_outputs/square_008um/filtered_feature_bc_matrix.h5')
adata.var_names_make_unique()
adata = sc.read_h5ad('hd.adata')

# load coordinates for each 8-um bin
coords = pd.read_parquet('spaceranger_output/outs/binned_outputs/square_008um/spatial/tissue_positions.parquet')
coords.index = coords['barcode']
coords = coords.drop(columns='barcode')
coords.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]

# add the spatial coordinates to the adata object
adata.obs = pd.merge(adata.obs, coords, how="left", left_index=True, right_index=True)
adata.obsm['spatial'] = adata.obs[["pxl_row_in_fullres", "pxl_col_in_fullres"]].values

# add H&E image and scalefactors to the adata object
img = np.array(Image.open("spaceranger_output/outs/binned_outputs/square_008um/spatial/tissue_hires_image.png"))
adata.uns['spatial'] = {
    "HD": {
        "images": {
            "hires": img
        }
    }
}
adata.uns['spatial']['HD']['scalefactors'] = json.loads(Path('spaceranger_output/outs/binned_outputs/square_008um/spatial/scalefactors_json.json').read_bytes())

print("\n=== AnnData object ===")
print(adata)

print("\n=== Gene names ===")
print(adata.var_names)

print("\n=== Properties of each gene ===")
print(adata.var.head())

print("\n=== Cell names ===")
print(adata.obs_names)

print("\n=== Properties of each cell ===")
print(adata.obs.head())

print("\n=== Coordinates of each cell ===")
print(adata.obsm['spatial'])

print("\n=== AnnData uns ===")
print(adata.uns)

# load cell segmentation results
gdf = gpd.read_file(f"spaceranger_output/outs/segmented_outputs/cell_segmentations.geojson")
gdf["centroid_x"] = gdf.centroid.x
gdf["centroid_y"] = gdf.centroid.y
# each row represents a cell, and the geometry column contains the polygon of the cell
print("\n=== First cell geometry ===")
print(gdf.loc[0, 'geometry'])

gdf_subset = gdf[
    (gdf['centroid_x'] >= 10000) & (gdf['centroid_x'] <= 10500) &
    (gdf['centroid_y'] >= 10000) & (gdf['centroid_y'] <= 10500)
]
fig, ax = plt.subplots(figsize=(4, 4))
for geom in gdf_subset.geometry:
    coords = list(geom.exterior.coords)
    patch = Polygon(coords, closed=True, facecolor='skyblue', edgecolor='blue', alpha=0.5)
    ax.add_patch(patch)
ax.set_aspect('equal')
ax.autoscale_view()
ax.set_title('Segmented cells')
fig.savefig(output_dir / "segmented_cells_subset.png", bbox_inches="tight", dpi=150)
plt.close(fig)  

# standards QC metrics with pp.calculate_qc_metrics and percentage of mitochondrial read counts per sample.
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

# QC metrics histograms
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
sns.histplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
sns.histplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)
save_plot(fig, "qc_metrics_histograms")

# perform some basic filtering of spots based on total counts and expressed genes
sc.pp.filter_cells(adata, min_counts=100)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs["pct_counts_mt"] < 10].copy()

# proceed to normalize Visium counts data with the built-in normalize_total method from Scanpy, and detect highly-variable genes
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)

# reduce the dimensionality of the data by running principal component analysis (PCA), which reveals the main axes of variation and denoises the data.
sc.pp.pca(adata)

# compute the neighborhood graph of cells using the PCA representation of the data matrix
sc.pp.neighbors(adata)

# embedding the graph in two dimensions using UMAP
sc.tl.umap(adata)

# leiden clustering
sc.tl.leiden(
    adata, key_added="clusters", n_iterations=2, resolution=0.2
)
# plot the clusters
plt.rcParams["figure.figsize"] = (4, 4)
sc.pl.umap(adata, color=["clusters"], wspace=0.4, show=False)
plt.savefig(output_dir / "umap_clusters.png", bbox_inches="tight", dpi=150)
plt.close()

# compute a ranking for the highly differential genes in each cluster.
sc.tl.rank_genes_groups(
    adata,
    groupby='clusters',  
    method='wilcoxon',   
    pts=True             
)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
plt.savefig(output_dir / "rank_genes_groups.png", bbox_inches="tight", dpi=150)
plt.close()

# assign cell types for each cluster
adata.obs['cell_type'] = adata.obs['clusters'].map({
    '0': 'Fibroblast',
    '1': 'Epithelial',
    '2': 'Epithelial'
})

# check the expression of COL1A2 and EPCAM across cells
sc.pl.umap(adata, color=["COL1A2","EPCAM","cell_type"], wspace=0.4)
plt.savefig(output_dir / "umap_col1a2_epcam.png", bbox_inches="tight", dpi=150)
plt.close()

# We will overlay the spots on top of the Hematoxylin and eosin stain (H&E) image to see how the cell types are spatially organized, using the function sc.pl.spatial
sc.pl.spatial(adata, img_key="hires", color=["cell_type", "EPCAM", "COL1A2"], wspace=0.4, show=False)
plt.savefig(output_dir / "spatial_celltype_epcam_col1a2.png", bbox_inches="tight", dpi=150)
plt.close()

# clustering resolution (granularity) in single-cell analysis controls how finely or coarsely the cells are grouped into clusters.
sc.tl.leiden(
    adata, key_added="clusters", n_iterations=2, resolution = 0.6
)
plt.rcParams["figure.figsize"] = (4, 4)
sc.tl.leiden(adata, key_added="clusters_highres", n_iterations=2, resolution=0.6)
sc.pl.umap(adata, color=["clusters_highres"], show=False)
plt.savefig(output_dir / "umap_clusters_highres.png", bbox_inches="tight", dpi=150)
plt.close()
sc.pl.spatial(adata, img_key="hires", color=["clusters_highres"], show=False)
plt.savefig(output_dir / "spatial_clusters_highres.png", bbox_inches="tight", dpi=150)
plt.close()



