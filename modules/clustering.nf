process RUNSEURAT {
  maxForks 10
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
  maxForks 10
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
  Rscript -e "rmarkdown::render('${celesta_script}', output_file='celesta_report_${roi}.html')" $all_markers $configs
  """
}

process RUNSCIMAP {
  maxForks 10
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
  path "matrixplot_${roi}.png", emit: matrixplot
  path "spatialplot_${roi}.png", emit: spatialplot
  path "umap_${roi}.png", emit: umap
    
  script:
  roi = all_markers.baseName.replace("all_markers_clean_", "")
  
  """
  python scimap_clustering.py
  """
}

process SCIMAPREPORT {
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