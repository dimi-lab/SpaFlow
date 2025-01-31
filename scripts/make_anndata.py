# Load necessary libraries
import sys, os
import pandas as pd
import anndata as ad
import numpy as np
import glob

## Should be able to store very large sizes:
## https://github.com/scverse/scanpy_usage/tree/master/170522_visualizing_one_million_cells

sample = sys.argv[1]

def clean_headers(df):
    header = [e.replace(':', '') for e in df.columns.values.tolist()]
    header = [e.replace('/', '') for e in header]
    header = [e.replace('^', '') for e in header]
    header = [e.replace('.', '') for e in header]
    header = [e.replace('µ', 'u') for e in header]
    header = [e.replace(' ', '_') for e in header]
    header = [e.replace('-02_', '_') for e in header]
    df.columns = header
    return df
            
# 1. Load all_markers_clean
markers_file = glob.glob(f"all_markers_clean_{sample}.csv")[0]  # Ensure only one file is matched
originaldf = pd.read_csv(markers_file, sep=',', low_memory=False) # changed to comma separated
all_markers_clean_means = originaldf.filter(regex=": Mean|: Min") # changed filter to like

originaldf = clean_headers(originaldf)
all_markers_clean_means = clean_headers(all_markers_clean_means)


# 2. Create AnnData object
adata = ad.AnnData(all_markers_clean_means, dtype=all_markers_clean_means.values.dtype) # added dtype
adata.var_names = all_markers_clean_means.columns.to_list()


# 3. Add spatial coordinates and other obsm data
x_col = 'Centroid_X_um' if 'Centroid_X_um' in originaldf.columns else 'x'
y_col = 'Centroid_Y_um' if 'Centroid_Y_um' in originaldf.columns else 'y'
adata.obsm["spatial"] = originaldf[[x_col, y_col]].to_numpy()

area_cols = ['Nucleus_Area_um2', 'Cell_Area_um2', 'Cell_Length_um', 'Cell_Circularity', 'Cell_Solidity', 'Cell_Max_diameter_um','Cell_Min_diameter_um']
for col in area_cols:
    if col in originaldf.columns:
        sNom = col.replace('_um2','').replace('_um','')
        adata.obsm[sNom] = originaldf[[col]].to_numpy()
        continue #stop on first hit

# 4. Add metadata from other CSVs
metacluster_files = glob.glob(f"*metaclusters_{sample}*.csv")
for file in metacluster_files:

    df_meta = pd.read_csv(file, sep=',', low_memory=False) # changed to comma separated
    df_meta = clean_headers(df_meta)

    x_col_meta = 'X_centroid' if 'X_centroid' in df_meta.columns else 'x'
    y_col_meta = 'Y_centroid' if 'Y_centroid' in df_meta.columns else 'y'

    #check if centroids match
    merged_df = pd.merge(originaldf[[x_col, y_col]], df_meta[[x_col_meta, y_col_meta]], left_on=[x_col,y_col], right_on=[x_col_meta, y_col_meta], how='left')
    if merged_df.shape[0] != originaldf.shape[0]:
        raise ValueError(f"Centroids in {file} do not match all_markers_clean centroids")

    for col in df_meta.columns:
        if "cluster" in col.lower(): #case insensitive matching
            adata.obs[col] = df_meta[col].values

# 5. Add imageid and spatial coordinates to obs
adata.uns["Slide"] = sample # changed to sample
adata.obs["imageid"] = pd.Categorical(originaldf["roi"])
adata.obs["X"] = adata.obsm['spatial'][:, 0]
adata.obs["Y"] = adata.obsm['spatial'][:, 1]

print(adata)

adata.write_h5ad(f"all_meta_clustering_{sample}.h5ad")
print(f"AnnData object for {sample} saved.")
