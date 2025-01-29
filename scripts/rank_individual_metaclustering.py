import os, sys
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import scanpy as sc

sample = sys.argv[1]
output_pdf_file = f"meta_ranking_{sample}.pdf"
adata = sc.read_h5ad(f"all_meta_clustering_{sample}.h5ad") # Replace with your file name

# 1. Check if var is empty
print(adata)

"""Explores, summarizes, and creates a PDF report for an AnnData object."""
with PdfPages(output_pdf_file) as pdf:
    # 1. Basic Information
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.text(0.05, 0.95, str(adata), fontsize=10, va='top')  # Print AnnData structure
    ax.axis('off')
    pdf.savefig(fig)
    plt.close(fig)

    # 2. Observation (obs) Summary
    fig, ax = plt.subplots(figsize=(10, 8))
    obs_summary = adata.obs.describe(include='all')
    ax.table(cellText=obs_summary.values, colLabels=obs_summary.columns, loc='center', cellLoc='center')
    ax.axis('off')
    ax.set_title("Observation (obs) Summary")
    pdf.savefig(fig)
    plt.close(fig)

    # 3. Variable (var) Summary
    fig, ax = plt.subplots(figsize=(10, 8))
    var_summary = adata.var.describe(include='all')
    ax.table(cellText=var_summary.values, colLabels=var_summary.columns, loc='center', cellLoc='center')
    ax.axis('off')
    ax.set_title("Variable (var) Summary")
    pdf.savefig(fig)
    plt.close(fig)

    # 4. Distribution of Key Metadata (e.g., cell type, batch)
    for col in adata.obs.columns:
        if adata.obs[col].dtype == 'category' or adata.obs[col].dtype == 'object': #check if it is categorical
            fig, ax = plt.subplots(figsize=(8, 6))
            sns.countplot(x=col, data=adata.obs, ax=ax)
            ax.set_title(f"Distribution of {col}")
            ax.tick_params(axis='x', rotation=45) #rotate labels if needed
            pdf.savefig(fig)
            plt.close(fig)

    # 5. Expression Distribution (Histograms or Boxplots)
    fig, axes = plt.subplots(2, 2, figsize=(10, 8)) # Example: 4 genes
    genes_to_plot = adata.var_names[:4]  # Adjust as needed
    for i, gene in enumerate(genes_to_plot):
        row = i // 2
        col = i % 2
        sns.histplot(adata[:, gene].X.flatten(), ax=axes[row, col]) # Plotting the expression values
        axes[row, col].set_title(f"Expression of {gene}")

    pdf.savefig(fig)
    plt.close(fig)

    # 6. Spatial Visualization (if spatial data is present)
    if "spatial" in adata.obsm.keys() and "Slide" in adata.uns.keys():  # Check for spatial data
        fig, ax = plt.subplots(figsize=(8, 6))
        sc.pl.spatial(adata, img_key="Slide", color="cell_type", ax=ax, show=False) # Example: color by cell type
        ax.set_title("Spatial Distribution")
        pdf.savefig(fig)
        plt.close(fig)

    # 7. UMAP or t-SNE Visualization
    sc.pp.neighbors(adata) #compute neighbors for umap
    sc.tl.umap(adata) #computes the umap
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color="cell_type", ax=ax, show=False) # Example: color by cell type
    ax.set_title("UMAP Visualization")
    pdf.savefig(fig)
    plt.close(fig)

    # 8. Highly Variable Genes
    sc.pp.highly_variable_genes(adata)
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.highly_variable_genes(adata, ax=ax, show=False)
    ax.set_title("Highly Variable Genes")
    pdf.savefig(fig)
    plt.close(fig)

    # 9. Add more visualizations or analyses as needed...




print(f"Report saved to {output_pdf_file}")
