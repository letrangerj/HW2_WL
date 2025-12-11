# %%
import json
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from PIL import Image
from scipy.io import mmwrite

warnings.filterwarnings("ignore", category=FutureWarning, module="scanpy")
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

print("Loading Data...")


## Load Data
def load_sample(sample_id, base_dir="data/GSE272362/HM_files/"):
    data_dir = Path(base_dir) / sample_id

    adata = sc.read_h5ad(data_dir / "expression.h5ad")
    adata.var_names_make_unique()  # unique gene names

    ## Images/coordinates
    spatial_dir = data_dir / "spatial"
    tissue_positions_path = list(spatial_dir.glob("*_tissue_positions.csv"))[0]
    scalefactors_path = list(spatial_dir.glob("*_scalefactors_json.json"))[0]
    image_path = list(spatial_dir.glob("*_tissue_hires_image.png"))[0]

    positions = pd.read_csv(tissue_positions_path)
    positions.set_index("barcode", inplace=True)
    # Merge spatial info into adata.obs (left join to keep all spots)
    adata.obs = adata.obs.merge(
        positions, left_index=True, right_index=True, how="left"
    )

    with open(scalefactors_path, "r") as f:
        scalefactors = json.load(f)

    image = np.array(Image.open(image_path))

    # Store spatial information in standard Scanpy format
    adata.uns["spatial"] = {
        sample_id: {"images": {"hires": image}, "scalefactors": scalefactors}
    }

    # Use pixel coordinates for accurate spatial plotting
    adata.obsm["spatial"] = adata.obs[
        ["pxl_row_in_fullres", "pxl_col_in_fullres"]
    ].to_numpy()

    adata.obs["sample_id"] = sample_id

    return adata


sample_ids = [f"HM33{i}" for i in range(32, 62)]
adatas = [load_sample(id) for id in sample_ids]

adata = sc.concat(adatas, label="batch", keys=sample_ids)
adata.obs_names_make_unique()  # 混合时会使得obs名字重复
adata.var_names_make_unique()

print("Load data complete")
print(f"Loaded: {adata.n_obs} spots, {adata.n_vars} genes")

# %%
# QC

## 输出位置
figures_dir = Path("./figures")
h5ad_dir = Path("./results")

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


## Plot processing
def plot_qc(adata, group):
    plt.rcParams["figure.figsize"] = (5, 5)
    sc.pl.spatial(
        adata,
        img_key="hires",
        spot_size=175,  # spot_diameter_fullres from scalefactors
        color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
        ncols=3,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(f"./figures/{group}_QC.png")

    # Add histograms for QC metrics distribution
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    axes[0].hist(
        adata.obs["total_counts"],
        bins=50,
        alpha=0.7,
        color="steelblue",
        edgecolor="black",
    )
    axes[0].axvline(3, color="red", linestyle="--", label="Filter threshold")
    axes[0].set_title("Total counts per spot")
    axes[0].set_xlabel("Total counts")
    axes[0].set_ylabel("Frequency")
    axes[0].legend()

    axes[1].hist(
        adata.obs["n_genes_by_counts"],
        bins=50,
        alpha=0.7,
        color="seagreen",
        edgecolor="black",
    )
    axes[1].axvline(3, color="red", linestyle="--", label="Filter threshold")
    axes[1].set_title("Genes per spot")
    axes[1].set_xlabel("Number of genes")
    axes[1].legend()

    axes[2].hist(
        adata.obs["pct_counts_mt"],
        bins=50,
        alpha=0.7,
        color="darkorange",
        edgecolor="black",
    )
    axes[2].set_title("Mitochondrial percentage")
    axes[2].set_xlabel("% MT counts")

    plt.tight_layout()
    plt.savefig(f"./figures/{group}_qc_histograms.png", dpi=300, bbox_inches="tight")
    plt.close()
    return


## 预处理函数
def process_qc(adata):
    """QC处理"""
    print("Processing QC metrics...")
    # Filter genes and cells
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=500)
    adata = adata[
        adata.obs["n_genes_by_counts"] > 300, :
    ]  # Filter out extreme outliers

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    # PCA
    sc.pp.pca(adata)

    return adata


## Plot before/after processing
plot_qc(adata, "1b_before")

# RCTD 需要原始的表达数据（整数），储存适合其的输入文件，参考ChatGPT
print("\nSaving RCTD-compatible raw counts (before normalization)...")
rctd_dir = h5ad_dir / "rctd_input"
rctd_dir.mkdir(parents=True, exist_ok=True)

adata_raw = adata.copy()

# Save raw count matrix and metadata for RCTD
import scipy.sparse as sp

raw_counts = sp.csr_matrix(adata_raw.X) if not sp.issparse(adata_raw.X) else adata_raw.X
sp.save_npz(rctd_dir / "spatial_counts_raw.npz", raw_counts)

raw_coo = raw_counts.tocoo()
mmwrite(str(rctd_dir / "spatial_counts_raw.mtx"), raw_coo)

# Save gene names (rownames for RCTD - genes as rows)
pd.DataFrame({"gene": adata_raw.var_names}).to_csv(
    rctd_dir / "genes.csv", index=False, header=False
)

# Save spot barcodes with metadata
adata_raw.obs.to_csv(rctd_dir / "spot_metadata.csv")

# Save spatial coordinates
spatial_coords_rctd = pd.DataFrame(
    adata_raw.obsm["spatial"],
    columns=["pxl_row_in_fullres", "pxl_col_in_fullres"],
    index=adata_raw.obs_names,
)
spatial_coords_rctd.to_csv(rctd_dir / "spatial_coords_rctd.csv")

print(f"RCTD raw count files saved to {rctd_dir}")

# 现在处理并绘制QC后图像
adata = process_qc(adata)
plotgqc(adata, "1b_after")

# Save processed data
adata.write(h5ad_dir / "processed_data.h5ad")
