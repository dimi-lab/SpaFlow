import os, sys
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score
import scanpy as sc

# Settings
sc.settings.verbosity = 1
sc.settings.figdir = "./"

sample = sys.argv[1]
output_pdf_file = f"meta_ranking_{sample}.pdf"
adata = sc.read_h5ad(f"all_meta_clustering_{sample}.h5ad") # Replace with your file name
feature_to_cluster_correlate_threshold = 4

# What does a high CH index score mean?
# A high CH index score indicates that the clusters are well separated and compact 
# A high CH index score indicates that the clustering solution has low within variability and high between cluster variability 

# Silhouette Interpretation
# High score = A score close to 1 indicates that the data points are well-matched to their cluster and poorly matched to other clusters. 
# Zero score = A score close to zero indicates that clusters overlap or that data points are close to multiple clusters. 
# Negative score = A score close to -1 indicates that data points are assigned to the wrong cluster. 

# Davies-Bouldin score:
# Lower is better: Unlike many metrics, a lower Davies-Bouldin score is desirable, signifying better cluster separation. 

select_data = adata.X
mean_cols = [i for i, col in enumerate(adata.var_names) if "Mean" in col]
select_data = select_data[:, mean_cols]
obs_table = adata.obs

# 1. Split cluster column and drop rows with "A"
if 'cluster' in obs_table.columns: #check if column exists
    obs_table['cluster_type'] = obs_table['cluster'].str.rsplit('_', 1, expand=True)[1] # added expand=True
    rows_to_drop = obs_table['cluster_type'] == 'A'
    obs_table = obs_table[~rows_to_drop]
    select_data = select_data[~rows_to_drop]  # Filter select_data accordingly
    obs_table = obs_table.drop(columns=['cluster_type', 'cluster']) # drop the column

# 2. Evaluate clusterings
scores = {}
clust_cols = [col for col in obs_table.columns if 'clust' in col]

for cl in clust_cols:
    try:
        labels = obs_table[cl].astype(str) # added astype(str) to handle mixed datatypes
        n_clusters = len(np.unique(labels)) # count unique labels
        if n_clusters > 1 and n_clusters < len(labels): # check if there is more than 1 cluster and less than the total rows
            scores[cl] = {
                "Calinski-Harabasz": calinski_harabasz_score(select_data, labels),
                "Silhouette": silhouette_score(select_data, labels),
                "Davies-Bouldin": davies_bouldin_score(select_data, labels)
            }
        else:
            print(f"Clustering {cl} has only one cluster or all rows belong to one cluster, skipping...")
    except Exception as e:
        print(f"Error evaluating clustering {cl}: {e}")
        scores[cl] = {"Calinski-Harabasz": np.nan, "Silhouette": np.nan, "Davies-Bouldin": np.nan} # set to nan to avoid crashing
        

# 3. Rank clusterings (example: based on Calinski-Harabasz)
ranked_clusterings = sorted(scores.items(), key=lambda item: item[1]["Calinski-Harabasz"] if not np.isnan(item[1]["Calinski-Harabasz"]) else -np.inf, reverse=True) 

#print("Ranked Clusterings (Calinski-Harabasz):")
#for cl, metrics in ranked_clusterings:
#    print(f"{cl}: {metrics}")



#4. """Normalizes scores to a 0-1 range for each metric."""
normalized_scores = {}
for metric in scores[list(scores.keys())[0]].keys():  # Iterate over metrics
    metric_values = [scores[cl][metric] for cl in scores]
    min_val = np.nanmin(metric_values) # handle nan values
    max_val = np.nanmax(metric_values) # handle nan values
    if max_val - min_val == 0: # avoid division by zero
        normalized_scores[metric] = {cl: 0 for cl in scores}
    else:
        normalized_scores[metric] = {
            cl: (scores[cl][metric] - min_val) / (max_val - min_val) if not np.isnan(scores[cl][metric]) else 0 for cl in scores
        }

## Davies-Bouldin score needs to be inverted to match other scores positive directions
inverted_scores = normalized_scores['Davies-Bouldin'].copy()
for cluster, score in inverted_scores.items():
    if not np.isnan(score): # handle nan values
        inverted_scores[cluster] = 1 - score
    else:
        inverted_scores[cluster] = score # keep nan values as nan
normalized_scores['Davies-Bouldin'] = inverted_scores


#5.  """Calculates the average of normalized scores for each clustering."""
average_scores = {}
for cl in normalized_scores[list(normalized_scores.keys())[0]].keys(): # iterate over the clusters
    valid_scores = [normalized_scores[metric][cl] for metric in normalized_scores if not np.isnan(normalized_scores[metric][cl])] # remove nan values from the list
    if valid_scores: # check if the list is not empty
        average_scores[cl] = np.mean(valid_scores)
    else:
        average_scores[cl] = 0 # if all scores are nan, set average to 0



#6. """Generates a plot visualizing the ranking."""
cluster_names = [cl for cl, metrics in ranked_clusterings]

# Original Scores (Log-transformed Calinski-Harabasz)
ch_scores = [scores[cl]["Calinski-Harabasz"] for cl in cluster_names]
# Handle potential negative or zero values before log transformation
ch_scores = np.array(ch_scores) # convert to numpy array
ch_scores[ch_scores <= 0] = np.nan # replace 0 or negative values with nan
ch_scores = np.log10(ch_scores)  # Log transform (using log1p to avoid log(0))
    
sil_scores = [scores[cl]["Silhouette"] for cl in cluster_names]

db_scores = [scores[cl]["Davies-Bouldin"] for cl in cluster_names]
# Handle potential negative or zero values before log transformation
db_scores = np.array(db_scores) # convert to numpy array
db_scores[db_scores <= 0] = np.nan # replace 0 or negative values with nan
db_scores = np.log10(db_scores)  # Log transform (using log1p to avoid log(0))

# Normalized Scores
norm_ch_scores = [normalized_scores["Calinski-Harabasz"][cl] for cl in cluster_names]
norm_sil_scores = [normalized_scores["Silhouette"][cl] for cl in cluster_names]
norm_db_scores = [normalized_scores["Davies-Bouldin"][cl] for cl in cluster_names]

# Average Scores
avg_scores = [average_scores[cl] for cl in cluster_names]

x = np.arange(len(cluster_names))  # the label locations
width = 0.15  # the width of the bars
fig, ax1 = plt.subplots(figsize=(15, 8))
ax2 = ax1.twinx()
rects1 = ax1.bar(x - width, ch_scores, width, label='Calinski-Harabasz', color = 'skyblue')
rects2 = ax1.bar(x, sil_scores, width, label='Silhouette', color = 'lightcoral')
rects3 = ax1.bar(x + width, db_scores, width, label='Davies-Bouldin', color = 'lightgreen')
rects4 = ax2.plot(x, avg_scores, label='Average Score', color = 'black', marker = 'o')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax1.set_ylabel('Original Scores (log10)')
ax2.set_ylabel('Normalized and Average Score')
ax1.set_title(f'Clustering Performance - {sample}')
ax1.set_xticks(x)
ax1.set_xticklabels(cluster_names, rotation=45, ha = 'right')
ax1.legend(loc = 'upper center')
ax2.legend(loc = 'upper right')

fig.tight_layout()
plt.savefig(f"clustering_ranking_{sample}.png", dpi=300, bbox_inches='tight')
# plt.show()

# 7. Output metrics:
data = []
for cluster_name, metrics in ranked_clusterings:
    row = metrics.copy()  # Create a copy to avoid modifying the original
    row['Cluster'] = cluster_name
    if cluster_name in average_scores: #check if cluster name exists in average scores
        row['Average_Score'] = average_scores[cluster_name]
    else:
        row['Average_Score'] = np.nan # if not, add a nan value
    data.append(row)

ranked_df = pd.DataFrame(data)
# Reorder columns (optional, but makes it look nicer)
cols = ['Cluster', 'Average_Score'] + [col for col in ranked_df.columns if col not in ['Cluster', 'Average_Score']]
ranked_df = ranked_df[cols]
ranked_df.to_csv(f"ranked_clusterings_{sample}.csv", index=False)


#### Section 2. Keeping top auto cluster...plot visuals
highest_cluster = max(average_scores, key=average_scores.get)
print(f"Highest scoring cluster: {highest_cluster}")

if highest_cluster not in adata.obs.columns:
    print(f"Warning: Cluster column '{highest_cluster}' not found in adata.obs. Skipping spatial plot.")

fig, ax = plt.subplots(figsize=(9, 8))
try:
    sc.pl.spatial(adata, color=highest_cluster, ax=ax, show=False, spot_size=10)  # Color by highest scoring cluster
    ax.set_title(f"Colored by {highest_cluster} - {sample}")
    plt.savefig(f"spatial_highest_score_pointplot_{sample}.png")
    # plt.show()
    plt.close(fig) # close the figure to prevent it from being displayed
    print(f"Spatial plot saved to spatial_highest_score_{sample}.png")
except Exception as e:
    print(f"Error creating spatial plot: {e}")
    

"""Calculates the top N positively correlated features to each cluster."""
top_labels = obs_table[highest_cluster]
unique_labels = np.unique(top_labels)
obs_table.index = obs_table.index.astype(int)  # Convert index to integers
top_features = {}
heatmap_data = []
for label in unique_labels:
    # Correctly get indices from obs_table
    ii = obs_table[obs_table[highest_cluster] == label].index.to_numpy()
    cluster_data = select_data[ii]
    if cluster_data.shape[0] < 1:
        continue  # Skip if no data for this cluster

    cluster_mean = np.mean(cluster_data, axis=0)    
    correlations = np.corrcoef(cluster_data, rowvar=False)[0, 1:] # Correct calculation
    correlations = np.nan_to_num(correlations)  # Handle NaNs
    top_indices = np.argsort(correlations)[::-1][:feature_to_cluster_correlate_threshold]
    top_features[label] = top_indices
    heatmap_data.append(cluster_mean)

print(top_features)

"""Generates a heatmap of the top features vs. the clusters."""
all_top_features = np.unique(np.concatenate(list(top_features.values())).flatten())
allVarMeanNames = adata.var_names[mean_cols]

heatmap_df = pd.DataFrame(heatmap_data, index=unique_labels, columns=allVarMeanNames) # create dataframe for the heatmap
namedTopFeatures = allVarMeanNames[ all_top_features ]
heatmap_df = heatmap_df[namedTopFeatures]

plt.figure(figsize=(10, 12))
sns.heatmap(heatmap_df, cmap="vlag", annot=True) # plot heatmap
plt.title(f"Heatmap of Top Correlated Features")
plt.xlabel("Top Correlated Features")
plt.xticks(rotation=25, ha='right')  # Rotate x-axis labels
plt.ylabel(highest_cluster)
plt.savefig(f"heatmap_top_correlations_{sample}.png")
#plt.show()


#### Section #3. Dig into overly expressed Min values?   ####


