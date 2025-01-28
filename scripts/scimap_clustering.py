# Optimized and refactored code with warnings addressed

# Load necessary libraries
import sys
import os
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
import scimap as sm
from sklearn.metrics import silhouette_score, calinski_harabasz_score
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from sklearn.cluster import KMeans
from kneed import KneeLocator

# Settings
sc.settings.verbosity = 1
sc.settings.figdir = "./"

# Constants
MIN_K = 3
MAX_K = 26
desired_k_method = "Calinski-Harabasz Score"

# Load configurations
marker_configs = pd.read_csv("marker_configs.csv")
marker_subset = marker_configs['marker'].tolist()

configs = pd.read_csv("configs.csv")
clustering_res = float(configs.loc[configs['object'] == "scimap_resolution", 'value'].item())

# Load ROI data
roi_path = [x for x in os.listdir() if x.startswith("all_markers_clean")][0]
roi_df = pd.read_csv(roi_path)
roi = os.path.splitext(roi_path.replace("all_markers_clean_", ""))[0]

# Preprocess data
no_bad = roi_df[~roi_df['qc'].str.contains("Artifact", na=False)]
single_vars = no_bad.filter(regex='(Cell: Median)', axis=1)
single_vars.columns = single_vars.columns.str.replace(": Cell: Median", "")
single_vars = single_vars.filter(regex='^(?!DAPI).')

# Check marker compatibility
diff_markers = set(marker_subset).difference(single_vars.columns)
if diff_markers:
    sys.exit("Config Contains a Marker NOT in the dataset: " + ', '.join(diff_markers))

# Create AnnData object
adata = ad.AnnData(single_vars[marker_subset])
adata.var_names = marker_subset
adata.obsm = {"spatial": no_bad[['Centroid X um', 'Centroid Y um']].to_numpy()}
adata.obs["imageid"] = pd.Categorical(no_bad["roi"])
adata.obs["X_centroid"] = no_bad[['Centroid X um']].to_numpy()
adata.obs["Y_centroid"] = no_bad[['Centroid Y um']].to_numpy()

# Optimal cluster determination
K = range(MIN_K, MAX_K)
inertia = []
scores_ch = {}

select_data = adata.X

for k in K:
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(select_data)
    inertia.append(kmeans.inertia_)
    scores_ch[k] = calinski_harabasz_score(select_data, labels)

kn = KneeLocator(K, inertia, curve='convex', direction='decreasing')
optimal_k_elbow = kn.knee
optimal_k_ch = max(scores_ch, key=scores_ch.get)

# Plot cluster evaluation
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].plot(K, inertia, marker='o')
axes[0].axvline(optimal_k_elbow, color='r', linestyle='--', label=f'Optimal k = {optimal_k_elbow}')
axes[0].set_title('Elbow Method (Inertia)')
axes[0].set_xlabel('Number of Clusters (k)')
axes[0].set_ylabel('Inertia')
axes[0].legend()
axes[0].grid(True)

axes[1].plot(K, scores_ch.values(), marker='o')
axes[1].axvline(optimal_k_ch, color='g', linestyle='--', label=f'Optimal k = {optimal_k_ch}')
axes[1].set_title('Calinski-Harabasz Scores')
axes[1].set_xlabel('Number of Clusters (k)')
axes[1].set_ylabel('CH Score')
axes[1].legend()
axes[1].grid(True)

plt.tight_layout()
plt.savefig(f"optimal_cluster_{roi}.png", dpi=300, bbox_inches='tight')

optimal_k = optimal_k_ch if desired_k_method == "Calinski-Harabasz Score" else optimal_k_elbow

# K-means clustering
kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init=10)
adata.obs['kmeans'] = kmeans.fit_predict(select_data)
adata.obs['kClusters'] = pd.Categorical(adata.obs['kmeans'])
sc.pl.spatial(adata, color='kClusters', spot_size=20, title='K-means Clustering', save=f'_kmeans_{roi}.png', show=False)

# Leiden clustering
adata = sm.tl.cluster(adata, method='leiden', subset_genes=marker_subset, resolution=clustering_res, use_raw=False)
sc.pl.spatial(adata, color='leiden', spot_size=20, title='Leiden Clustering', save=f'_leiden_{roi}.png', show=False)

# Calculate metrics
ari_score = adjusted_rand_score(adata.obs['kmeans'], adata.obs['leiden'])
nmi_score = normalized_mutual_info_score(adata.obs['kmeans'], adata.obs['leiden'])

results = {'ARI': ari_score, 'NMI': nmi_score, 'Elbow_N':optimal_k_elbow, 'CH_N':optimal_k_ch}
# Ensure all values in results are converted to native Python types
results = {key: int(value) if isinstance(value, np.integer) else value for key, value in results.items()}
with open(f'cluster_scores_{roi}.json', 'w') as fp:
    json.dump(results, fp, indent=4, ensure_ascii=False)

# UMAP
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=10)
sc.tl.umap(adata)
sc.pl.umap(adata, color=['kClusters'], cmap='vlag', use_raw=False, s=30, save=f'_kmeans_{roi}.png', show=False)
sc.pl.umap(adata, color=['leiden'], cmap='vlag', use_raw=False, s=30, save=f'_leiden_{roi}.png', show=False)

# Heatmaps
sc.pl.matrixplot(adata, var_names=adata.var.index, groupby='leiden', dendrogram=True,
                 use_raw=False, cmap="vlag", standard_scale='var', save=f'leiden_{roi}.png', show=False)
sc.pl.matrixplot(adata, var_names=adata.var.index, groupby='kClusters', dendrogram=True,
                 use_raw=False, cmap="icefire", standard_scale='var', save=f'kmeans_{roi}.png', show=False)

# Write out coordinates and clusters
# adata.obs.to_csv(f"scimap_clusters_{roi}.csv", index=False) # This doesn't include the dropped Artifact rows
aDf = pd.DataFrame(adata.obs)

allArtifacts = roi_df.loc[roi_df['qc'].str.contains("Artifact", na=True), ["roi", 'Centroid X um', 'Centroid Y um']]

allArtifacts.columns = ["imageid", "X_centroid", "Y_centroid"]
allArtifacts['kmeans'] = 'A'
allArtifacts['kClusters'] = 0
allArtifacts['leiden'] = 'A'


concatenated_df = pd.concat([aDf, allArtifacts], axis=0, ignore_index=True)
concatenated_df.to_csv(f"scimap_clusters_{roi}.csv", index=False)


