# Load necessary libraries
import sys
import os
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
sc.settings.verbosity = 1
sc.settings.figdir = "./"

# Import Scimap
import scimap as sm

from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans

# pip install kneed
from kneed import KneeLocator

min_k = 3
max_k = 38
desired_k_method = "Calinski-Harabasz Score"
# desired_k_method = "Elbow Method"

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
datasetMarkerNames = singleJustVars.columns.str.replace(": Cell: Median", "")
singleJustVars.columns = datasetMarkerNames
singleJustVars = singleJustVars.filter(regex='^((?!DAPI).)*$',axis=1)
# singleJustVars = singleJustVars[marker_subset]
# singleJustVars

diffMarkers = set(marker_subset).difference(datasetMarkerNames)
if len(diffMarkers) > 0:
    sys.exit("Config Contains a Marker NOT in the dataset: "+', '.join(diffMarkers) )

# Create anndata object
adata = ad.AnnData(singleJustVars) # create AnnData object
adata.var_names = singleJustVars.columns.to_list() # Set variable names
adata.obsm={ "spatial": noBad[['Centroid X um','Centroid Y um']].to_numpy() }# Add coordinates
adata.obs["imageid"] = pd.Categorical( noBad["roi"] ) # Add ROI number
adata.obs["X_centroid"] = noBad[['Centroid X um']].to_numpy()
adata.obs["Y_centroid"] = noBad[['Centroid Y um']].to_numpy()


 ##### Want to get Optimal # of Clusters #####
 #############################################
K = range(min_k, max_k)
# Elbow method to determine optimal k
inertia = []
scoresCH = {}

selectData = singleJustVars[marker_subset]

for k in K:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(selectData)
    inertia.append(kmeans.inertia_)
    score = calinski_harabasz_score(selectData, labels)
    scoresCH[k] = score

# Use KneeLocator to find the "elbow"
kn = KneeLocator(K, inertia, curve='convex', direction='decreasing')
optimal_k_elbow  = kn.knee

# Determine optimal k from Calinski-Harabasz scores (max score)
ch_scores = list(scoresCH.values())
optimal_k_ch = K[np.argmax(ch_scores)]

# Plotting side-by-side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Elbow plot
axes[0].plot(K, inertia, marker='o')
axes[0].axvline(x=optimal_k_elbow, color='r', linestyle='--', label=f'Optimal k = {optimal_k_elbow}')
axes[0].set_title('Elbow Method (Inertia)')
axes[0].set_xlabel('Number of Clusters (k)')
axes[0].set_ylabel('Inertia')
axes[0].legend()
axes[0].grid(True)

# Calinski-Harabasz plot
axes[1].plot(K, ch_scores, marker='o')
axes[1].axvline(x=optimal_k_ch, color='g', linestyle='--', label=f'Optimal k = {optimal_k_ch}')
axes[1].set_title('Calinski-Harabasz Scores')
axes[1].set_xlabel('Number of Clusters (k)')
axes[1].set_ylabel('CH Score')
axes[1].legend()
axes[1].grid(True)

plt.suptitle('Cluster Analysis: Elbow and Calinski-Harabasz Score Methods')
plt.tight_layout()
#plt.show()
plt.savefig("optimal_cluster_{roi}.png", dpi=300, bbox_inches='tight')

optimalKoverall = 3
if desired_k_method == "Calinski-Harabasz Score":
    print(f'Optimal k (Calinski-Harabasz Score): {optimal_k_ch}')
    optimalKoverall = optimal_k_ch    
else:
    print(f'Optimal k (Elbow Method):'.format(optimal_k_elbow))
    optimalKoverall = optimal_k_elbow


kmeans = KMeans(n_clusters=optimalKoverall, random_state=42, n_init=10)
adata.obs['kmeans'] = kmeans.fit_predict(selectData)
adata.obs['kClusters'] = pd.Categorical( adata.obs['kmeans'])
sc.pl.spatial(adata, color='kClusters', spot_size=20, title='K-means Clustering', save = f'_kmeans_{roi}.png', show = False)


# Run leiden clustering with selected genes
adata = sm.tl.cluster(adata, method = 'leiden', subset_genes = marker_subset, resolution = clustering_res, use_raw=False)

# Extract cluster labels
kmeans_labels = adata.obs['kmeans']
leiden_labels = adata.obs['leiden']

# Calculate Adjusted Rand Index (ARI) and Normalized Mutual Information (NMI)
ari_score = adjusted_rand_score(kmeans_labels, leiden_labels)
nmi_score = normalized_mutual_info_score(kmeans_labels, leiden_labels)

results = {
    'ARI': ari_score,
    'NMI': nmi_score
}

with open('cluster_scores_{roi}.json', 'w') as fp:
  json.dump(results, fp, indent=4, ensure_ascii=False)
  fp.close()


# Run UMAP - want to specify n_neighbors and n_pcs in configs?
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10) # Computing the neighborhood graph
sc.tl.umap(adata) # Build a UMAP to visualize the neighbourhood graph

# Plot the UMAP and save figure
sc.pl.umap(adata, color=['leiden'], cmap= 'vlag', use_raw=False, s=30, save = f'_{roi}.png', show=False)

# Plot heatmap of clusters vs markers
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='leiden', dendrogram=True, #### Add another with only subset genes
use_raw=False, cmap="vlag", standard_scale='var', save = f'_leiden_{roi}.png', show = False)

sm.pl.spatial_scatterPlot(adata, colorBy='leiden', s=2, outputDir = './', outputFileName = f'spatialplot_{roi}.png', show = False)

# Plot heatmap of clusters vs markers
sc.pl.matrixplot(adata, var_names= adata.var.index, groupby='kmeans', dendrogram=True, #### Add another with only subset genes
use_raw=False, cmap="icefire", standard_scale='var', save = f'_kmeans_{roi}.png', show = False)

# Write out coordinates and clusters
adata.obs.to_csv(f"scimap_clusters_{roi}.csv", index=False)












