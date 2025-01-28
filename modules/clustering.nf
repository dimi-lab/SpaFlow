process RUNSEURAT {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/seurat",
        pattern: "*.html",
        mode: "copy"
  )
    
  publishDir(
        path: "${params.output_dir}/output_tables/seurat",
        pattern: "seurat_clusters*.csv",
        mode: "copy"
  )
  
  if (params.export_intermediates) {
    publishDir(
          path: "${params.output_dir}/intermediates/seurat",
          pattern: "CLR_seurat_centroids*.csv",
          mode: "copy"
    )
  }
  
  input:
  path clustering_script
  path all_markers
  path configs
  path marker_configs
  
  output:
  path "seurat_report_${roi}.html"
  path "seurat_clusters_${roi}.csv", emit: seurat_clusters_noid
  path "CLR_seurat_centroids_${roi}.csv"
  tuple val(roi), path("seurat_clusters_${roi}.csv") , emit: seurat_clusters
    
  script:
  roi = all_markers.baseName.replace("all_markers_clean_", "")
  
  """
  Rscript -e "rmarkdown::render('${clustering_script}', output_file='seurat_report_${roi}.html')" $all_markers $configs $marker_configs
  """
}

process RUNCELESTA {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/celesta",
        pattern: "*.html",
        mode: "copy"
    )
    
  publishDir(
        path: "${params.output_dir}/output_tables/celesta",
        pattern: "*.csv",
        mode: "copy"
    )
  
  input:
  path celesta_script
  path all_markers
  path configs
  path celesta_prior_matrix
  
  output:
  path "celesta_report_${roi}.html"
  tuple val(roi), path("celesta_classes_${roi}.csv"), emit: celesta_classes
    
  script:
  roi = all_markers.baseName.replace("all_markers_clean_", "")
  
  """
  Rscript -e "rmarkdown::render('${celesta_script}', output_file='celesta_report_${roi}.html')" $all_markers $configs $celesta_prior_matrix
  """
}

process RUNSCIMAP {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_tables/scimap",
        pattern: "*.csv",
        mode: "copy"
  )
  
  input:
  path scimapscript
  path all_markers
  path configs
  path marker_configs
  
  output:
  path "scimap_clusters_${roi}.csv"
  tuple val(roi), path("scimap_clusters_${roi}.csv") , emit: scimap_clusters
  path "scimap_clusters_${roi}.csv", emit: scimap_clusters_noid
  path "matrixplot_kmeans_${roi}.png", emit: matrixplot_K
  path "matrixplot_leiden_${roi}.png", emit: matrixplot_L
  path "optimal_cluster_${roi}.png", emit: optimalplot
  path "show_kmeans_${roi}.png", emit: spatialplot_K
  path "show_leiden_${roi}.png", emit: spatialplot_L
  path "umap_kmeans_${roi}.png", emit: umap_K
  path "umap_leiden_${roi}.png", emit: umap_L
  path "cluster_scores_${roi}.json", emit: clustermeterics
    
  script:
  roi = all_markers.baseName.replace("all_markers_clean_", "")
  
  """
  python scimap_clustering.py
  """
}

process SCIMAPREPORT {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/scimap",
        pattern: "*.html",
        mode: "copy"
  )
  
  input:
  path scimap_report_script
  path matrixplot
  path spatialplot
  path umap
  
  output:
  path "*.html"
  
  script:
  roi = umap.baseName.replace("umap_", "")
  
  """
  Rscript -e "rmarkdown::render('${scimap_report_script}', 
                                  output_file='scimap_report_${roi}.html')"
  """  
}
