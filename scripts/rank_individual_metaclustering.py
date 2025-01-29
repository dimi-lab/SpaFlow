import os, sys
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.metrics import calinski_harabasz_score, silhouette_score, davies_bouldin_score
import scanpy as sc

sample = sys.argv[1]
output_pdf_file = f"meta_ranking_{sample}.pdf"
adata = sc.read_h5ad(f"all_meta_clustering_{sample}.h5ad") # Replace with your file name
# print(adata)

select_data = adata.X
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

print("Ranked Clusterings (Calinski-Harabasz):")
for cl, metrics in ranked_clusterings:
    print(f"{cl}: {metrics}")



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
ax1.legend(loc = 'upper left')
ax2.legend(loc = 'upper right')

fig.tight_layout()
plt.savefig(f"clustering_ranking_{sample}.png")
plt.show()



    

