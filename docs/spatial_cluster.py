#!/usr/bin/env python3
"""
Spatial Transcriptomics Analysis Pipeline
Usage: python spatial_analysis.py --path /path/to/visium/data
"""

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 非交互式环境使用
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Disable annoying warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(description='Spatial transcriptomics analysis pipeline')
    parser.add_argument('--path', type=str, required=True,
                       help='Path to Visium dataset directory')
    parser.add_argument('--output', type=str, default='./results',
                       help='Output directory for results (default: ./results)')
    parser.add_argument('--sample_id', type=str, default=None,
                       help='Sample ID for downloading dataset (optional)')
    parser.add_argument('--resolution', type=float, default=0.3,
                       help='Resolution for Leiden clustering (default: 0.3)')
    parser.add_argument('--n_clusters', type=int, default=10,
                       help='Number of clusters for GMM (default: 10)')
    parser.add_argument('--n_neighbors', type=int, default=8,
                       help='Number of neighbors for spatial graph (default: 8)')
    
    return parser.parse_args()

def create_output_directory(output_dir):
    """创建输出目录"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    return output_dir

def load_data(path, sample_id=None):
    """加载数据"""
    print(f"Loading data from: {path}")
    
    # 检查路径是否存在
    if not os.path.exists(path):
        print(f"Error: Path does not exist: {path}")
        
        # 如果提供了sample_id，尝试下载数据
        if sample_id:
            print(f"Attempting to download sample: {sample_id}")
            try:
                adata = sc.datasets.visium_sge(sample_id=sample_id)
                print(f"Successfully downloaded sample: {sample_id}")
                return adata
            except Exception as e:
                print(f"Error downloading sample: {e}")
                sys.exit(1)
        else:
            sys.exit(1)
    
    try:
        # 尝试读取Visium数据
        adata = sc.read_visium(path)
        print(f"Successfully loaded data with {adata.n_obs} spots and {adata.n_vars} genes")
        return adata
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

def process_data(adata, output_dir):
    """数据处理流程"""
    print("Processing data...")
    
    # 基本处理
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)
    
    # 标准化
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    
    # PCA
    sc.pp.pca(adata)
    
    # 保存原始数据
    # adata.write(os.path.join(output_dir, 'adata_processed.h5ad'))
    # print(f"Processed data saved to: {os.path.join(output_dir, 'adata_processed.h5ad')}")
    
    return adata

def create_qc_plots(adata, output_dir):
    """创建QC图"""
    print("Creating QC plots...")
    
    plt.rcParams["figure.figsize"] = (3, 3)
    fig = sc.pl.pca(adata, color=["total_counts", "n_genes_by_counts", "pct_counts_mt", "array_row"],
                    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)], ncols=2, size=2,
                    show=False, return_fig=True)
    fig.savefig(os.path.join(output_dir, 'pca_qc_plots.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    plt.rcParams["figure.figsize"] = (4, 4)
    fig = sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"],
                       show=False, return_fig=True)
    fig.savefig(os.path.join(output_dir, 'spatial_qc_plots.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

def spatial_clustering(adata, n_neighbors=8, resolution=0.3, output_dir='./'):
    """空间聚类分析"""
    print("Performing spatial clustering...")
    
    # 基于空间位置的聚类
    sc.pp.neighbors(adata, n_pcs=2, n_neighbors=n_neighbors, method='gauss', 
                    use_rep='spatial', metric='euclidean', key_added='loc')
    
    sc.tl.leiden(adata, neighbors_key='loc', resolution=resolution)
    
    # 保存leiden聚类结果图
    plt.rcParams["figure.figsize"] = (4, 4)
    fig = sc.pl.spatial(adata, color=["leiden"], show=False, return_fig=True)
    fig.savefig(os.path.join(output_dir, 'spatial_leiden_clusters.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # 保存聚类结果到CSV
    # adata.obs[['leiden']].to_csv(os.path.join(output_dir, 'leiden_clusters.csv'))
    
    return adata

def cellcharter_analysis(adata, n_neighbors=8, n_clusters=10, output_dir='./'):
    """CellCharter分析"""
    print("Performing CellCharter analysis...")
    
    try:
        import cellcharter as cc
        
        # 重新计算邻居（确保使用相同的参数）
        sc.pp.neighbors(adata, n_pcs=2, n_neighbors=n_neighbors, method='gauss', 
                       use_rep='spatial', metric='euclidean', key_added='loc')
        
        # CellCharter聚合
        cc.gr.aggregate_neighbors(adata, n_layers=2, use_rep='X_pca', 
                                 out_key='X_cellcharter', connectivity_key='loc')
        
        # 自动K聚类
        autok = cc.tl.ClusterAutoK(
            n_clusters=(2,12), 
            max_runs=3,
            convergence_tol=0.001
        )
        autok.fit(adata, use_rep='X_cellcharter')
        save_path = os.path.join(output_dir, "cellcharter_autok_stability.png")
        cc.pl.autok_stability(autok, save=save_path)
        adata.obs['cluster_cellcharter'] = autok.predict(adata, use_rep='X_cellcharter')
        
        # 保存自动聚类结果
        fig = sc.pl.spatial(adata, color=['cluster_cellcharter'], img=None,
                          size=1.5, show=False, return_fig=True)
        fig.savefig(os.path.join(output_dir, 'cellcharter_autok_clusters.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close(fig)
    
        # GMM聚类
        gmm = cc.tl.Cluster(
            n_clusters=n_clusters, 
            random_state=12345
        )
        gmm.fit(adata, use_rep='X_cellcharter')
        adata.obs[f'cluster_cellcharter_{n_clusters}'] = gmm.predict(adata, use_rep='X_cellcharter')
        
        # 保存GMM聚类结果
        fig = sc.pl.spatial(adata, color=[f'cluster_cellcharter_{n_clusters}'], 
                          img_key="hires", size=1.5, show=False, return_fig=True)
        fig.savefig(os.path.join(output_dir, f'cellcharter_gmm_{n_clusters}_clusters.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close(fig)
        
        # 保存CellCharter结果到CSV
        cellcharter_cols = [col for col in adata.obs.columns if 'cluster_cellcharter' in col]
        adata.obs[cellcharter_cols].to_csv(os.path.join(output_dir, 'cellcharter_clusters.csv'))
        
    except ImportError:
        print("Warning: cellcharter package not found. Skipping CellCharter analysis.")
        print("To install: pip install cellcharter")
    except Exception as e:
        print(f"Error in CellCharter analysis: {e}")

def save_results(adata, output_dir):
    """保存所有结果"""
    print("Saving results...")
    
    # 保存完整的AnnData对象
    # adata.write(os.path.join(output_dir, 'adata_final.h5ad'))
    
    # 保存观察数据（包含所有聚类结果）
    # adata.obs.to_csv(os.path.join(output_dir, 'all_clusters.csv'))
    
    # 保存变量数据
    # adata.var.to_csv(os.path.join(output_dir, 'genes_info.csv'))
    
    # 保存聚类统计
    cluster_stats = {}
    for col in adata.obs.columns:
        if 'cluster' in col.lower() or 'leiden' in col.lower():
            cluster_stats[col] = adata.obs[col].value_counts().to_dict()
    
    import json
    with open(os.path.join(output_dir, 'cluster_statistics.json'), 'w') as f:
        json.dump(cluster_stats, f, indent=2)
    
    print(f"All results saved to: {output_dir}")

def main():
    """主函数"""
    # 解析参数
    args = parse_arguments()
    
    # 创建输出目录
    output_dir = create_output_directory(args.output)
    
    # 加载数据
    adata = load_data(args.path, args.sample_id)
    
    # 处理数据
    adata = process_data(adata, output_dir)
    
    # 创建QC图
    create_qc_plots(adata, output_dir)
    
    # 空间聚类
    # adata = spatial_clustering(adata, args.n_neighbors, args.resolution, output_dir)
    
    # CellCharter分析
    cellcharter_analysis(adata, args.n_neighbors, args.n_clusters, output_dir)
    
    # 保存结果
    save_results(adata, output_dir)
    
    print("Analysis completed successfully!")
    print(f"Results saved in: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()

