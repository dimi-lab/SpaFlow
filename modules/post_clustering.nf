process RUNMETACLUSTERS {
  cpus 8
  memory '24 GB'


  publishDir(
        path: "${params.output_dir}/output_reports/metacluster",
        pattern: "*.html",
        mode: "copy"
  )
    
  publishDir(
        path: "${params.output_dir}/output_tables/metacluster",
        pattern: "seurat_metaclusters*.csv",
        mode: "copy"
  )
  
  if (params.export_intermediates) {
    publishDir(
          path: "${params.output_dir}/intermediates/metacluster",
          pattern: "arcsin_zscore_seurat_centroids.csv",
          mode: "copy"
    )
  }
  
  input:
  path seurat_metacluster_script
  path configs
  path marker_configs
  path allmarkers_collected
  path seurat_output
  
  output:
  path "arcsin_zscore_seurat_centroids.csv"
  path "seurat_metacluster_report.html"
  path "seurat_metaclusters*.csv"
  
  script:
  """
  Rscript -e "rmarkdown::render('${seurat_metacluster_script}', 
                                output_file='seurat_metacluster_report.html')"
  """
  
}

process RUNSOMCLUSTERS {
  cpus 8
  memory '24 GB'


  publishDir(
        path: "${params.output_dir}/output_reports/som_metacluster",
        pattern: "*.html",
        mode: "copy"
  )

  publishDir(
        path: "${params.output_dir}/output_tables/som_metacluster",
        pattern: "som_metaclusters*.csv",
        mode: "copy"
  )

  input:
  path seurat_som_script
  path configs
  path marker_configs
  path allmarkers_collected
  path allcentroids_collected
  path seurat_clusters_noid_collected
  val num_clusters
  

  output:
  path "som_metaclusters*.html"
  path "som_metaclusters*.csv", optional: true

  script:
  """
  Rscript -e "rmarkdown::render('${seurat_som_script}', params = list(som_num_clusters = '${num_clusters}'),
                                output_file='som_metaclusters_${num_clusters}_report.html')"
  """

}


process SEURATVCELESTA { // This script should output the final cluster files, so the previous ones don't need to output cluster files
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/seuratvcelesta",
        pattern: "*.html",
        mode: "copy"
    )
    
  publishDir(
        path: "${params.output_dir}/output_tables/seuratvcelesta",
        pattern: "*.csv",
        mode: "copy"
    )
  
  input:
  path seurat_vs_celesta_script
  tuple val(roi), path(seurat_clusters), path(celesta_classes)
  
  output:
  path "seurat_vs_celesta_report_${roi}.html"
  path "seurat_vs_celesta_clusters_*.csv"
  
  script:
  roi = celesta_classes.baseName.replace("celesta_classes_", "")
  
  """
  Rscript -e "rmarkdown::render('${seurat_vs_celesta_script}', 
                                output_file='seurat_vs_celesta_report_${roi}.html')"
  """
  
}

process SEURATVSCIMAP {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/seuratvscimap",
        pattern: "*.html",
        mode: "copy"
    )
  
  input:
  path seurat_vs_scimap_script
  tuple val(roi), path(seurat_clusters), path(scimap_clusters)
  
  output:
  path "seurat_vs_scimap_report_${roi}.html"
  
  script:
  roi = scimap_clusters.baseName.replace("scimap_clusters_", "")
  
  """
  Rscript -e "rmarkdown::render('${seurat_vs_scimap_script}', 
                                output_file='seurat_vs_scimap_report_${roi}.html')"
  """
  
}
