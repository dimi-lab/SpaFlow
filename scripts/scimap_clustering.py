# Load necessary libraries
import sys
import os
import anndata as ad
import pandas as pd
import scanpy as sc

# Import Scimap
import scimap as sm

sc.settings.figdir = "./"

# Read in selected markers
marker_configs = pd.read_csv("marker_configs.csv")
marker_subset = marker_configs['marker'].tolist()

# Read in configs
configs = pd.read_csv("configs.csv")
clustering_res = configs.loc[configs['object'] == "scimap_resolution", 'value'].item()
clustering_res = float(clustering_res)

# Read in the cleaned marker file
roi_path = [x for x in os.listdir() if x.startswith("all_markers_clean")][0]
roi_df = pd.read_csv(roi_path)

# Get ROI name for file paths later
roi = roi_path.replace("all_markers_clean_", "")
roi = roi.replace(".csv", "")

# Get the feature columns only, exclude DAPI
noBad = roi_df[~roi_df['qc'].str.contains("Artifact", na=False)] # remove artifacts from roi_df
singleJustVars = noBad.filter(regex='(Cell: Median)',axis=1) # Get only markers
singleJustVars.columns = singleJustVars.columns.str.replace(": Cell: Median", "")
singleJustVars = singleJustVars.filter(regex='^((?!DAPI).)*$',axis=1)
#singleJustVars = singleJustVars[marker_subset]
singleJustVars


# Create anndata object
adata = ad.AnnData(singleJustVars) # create AnnData object
adata.var_names = singleJustVars.columns.to_list() # Set variable names
adata.obsm={ "spatial": noBad[['Centroid X um','Centroid Y um']].to_numpy() # Add coordinates
           }  
adata.obs["imageid"] = pd.Categorical( noBad["roi"] ) # Add ROI number
adata.obs["X_centroid"] = noBad[['Centroid X um']].to_numpy()
adata.obs["Y_centroid"] = noBad[['Centroid Y um']].to_numpy()

# Run UMAP - want to specify n_neighbors and n_pcs in configs?
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata) # Build a UMAP to visualize the neighbourhood graph

# Run leiden clustering with selected genes
adata = sm.tl.cluster(adata, method = 'leiden', subset_genes = marker_subset, resolution = clustering_res, use_raw=False)

# Plot the UMAP and save figure
sc.pl.umap(adata, color=['leiden'], cmap= 'vlag', use_raw=False, s=30, save = f'_{roi}.png', show=False)

# Plot heatmap of clusters vs markers
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='leiden', dendrogram=True, #### Add another with only subset genes
use_raw=False, cmap="vlag", standard_scale='var', save = f'{roi}.png', show = False)

sm.pl.spatial_scatterPlot(adata, colorBy='leiden', s=2, outputDir = './', outputFileName = f'spatialplot_{roi}.png')

# Write out coordinates and clusters
adata.obs.to_csv(f"scimap_clusters_{roi}.csv", index=False)












