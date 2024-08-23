#!/usr/bin/env nextflow

params.celesta_prior_matrix = "${projectDir}/celesta_prior_matrix.csv"
params.qcscript = "${projectDir}/scripts/QC.Rmd"
params.collect_bin_density_script= "${projectDir}/scripts/collect_bin_density.Rmd"
params.collect_sigsum_script= "${projectDir}/scripts/collect_sigsum.Rmd"
params.seuratscript = "${projectDir}/scripts/seurat_clustering.Rmd"
params.celestascript = "${projectDir}/scripts/CELESTA_clustering.Rmd"
params.scimapscript = "${projectDir}/scripts/scimap_clustering.py"
params.scimap_report_script = "${projectDir}/scripts/scimap_report.Rmd"
params.seurat_metacluster_script = "${projectDir}/scripts/seurat_metaclustering.Rmd"
params.seurat_vs_celesta_script = "${projectDir}/scripts/seurat_vs_celesta.Rmd"
params.seurat_vs_scimap_script = "${projectDir}/scripts/seurat_vs_scimap.Rmd"

process WRITECONFIGFILE {
  input:
  val sigsum_quantile_high
  val sigsum_quantile_low
  val bin_size
  val density_cutoff
  val cluster_metric
  val clustering_res
  val min_clusters
  val min_res
  val max_res
  val res_step
  val scimap_resolution
  val min_metaclusters
  val max_metaclusters
  
  output:
  path "configs.csv", emit: configfile
  
  shell:
  """
  echo "object,value" > configs.csv
  
  echo "sigsum_quantile_high,!{sigsum_quantile_high}" >> configs.csv
  echo "sigsum_quantile_low,!{sigsum_quantile_low}" >> configs.csv
  echo "bin_size,!{bin_size}" >> configs.csv
  echo "density_cutoff,!{density_cutoff}" >> configs.csv
  echo "cluster_metric,!{cluster_metric}" >> configs.csv
  echo "clustering_res,!{clustering_res}" >> configs.csv
  echo "min_clusters,!{min_clusters}" >> configs.csv
  echo "min_res,!{min_res}" >> configs.csv
  echo "max_res,!{max_res}" >> configs.csv
  echo "res_step,!{res_step}" >> configs.csv
  echo "scimap_resolution,!{scimap_resolution}" >> configs.csv
  echo "min_metaclusters,!{min_metaclusters}" >> configs.csv
  echo "max_metaclusters,!{max_metaclusters}" >> configs.csv
  """
}

process WRITEMARKERFILE {
  input:
  val markers
  
  output:
  path "marker_configs.csv", emit: markerconfigfile
  
  shell:
  """
  list=!{markers}

  IFS=',' read -r -a array <<< "\$list"
  
  echo "marker" > marker_configs.csv

  for item in "\${array[@]}"; do
    echo \$item >> marker_configs.csv
  done
  """
}

process RUNQC {
  maxForks 10
  
  publishDir(
    path: "${params.output_dir}/output_tables/qc",
    pattern: "*.csv",
    mode: "copy"
  )

  input:
  path quantfile
  path QC_scriptpath
  path configs

	output:
	path "QC_report_${roi}.html"
	path "bin_density_${roi}.png", emit: bin_density
	path "sigsum_${roi}.png", emit: sigsum
  path "all_markers_clean_${roi}.csv", emit: all_markers
	
	script:
	roi = quantfile.baseName
	
	"""
	Rscript -e "rmarkdown::render('${QC_scriptpath}', output_file='QC_report_${quantfile.baseName}.html')"  $quantfile $configs
	"""
}

process COLLECTBINDENSITY {
  publishDir(
    path: "${params.output_dir}/output_reports/qc",
    pattern: "*.html",
    mode: "copy"
  )
  
  input:
  path collect_bin_density_script
  path bin_density_collected
  
  output:
  path "bin_density_report.html"
  
  script:
  """
  Rscript -e "rmarkdown::render('${collect_bin_density_script}', 
                                  output_file='bin_density_report.html')"
  """  
}

process COLLECTSIGSUM {
  publishDir(
        path: "${params.output_dir}/output_reports/qc",
        pattern: "*.html",
        mode: "copy"
  )

  
  input:
  path collect_sigsum_script
  path sigsum_collected
  
  output:
  path "sigsum_report.html"
  
  script:
  """
  Rscript -e "rmarkdown::render('${collect_sigsum_script}', 
                                  output_file='sigsum_report.html')"
  """  
}

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

process RUNMETACLUSTERS {

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

process RUNCOMPARISON { // This script should output the final cluster files, so the previous ones don't need to output cluster files

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


workflow {
  file_ch = Channel.fromPath(params.data_pattern)
  
  WRITEMARKERFILE(params.markers)
  WRITECONFIGFILE(params.sigsum_quantile_high,params.sigsum_quantile_low,
  params.bin_size,params.density_cutoff,params.cluster_metric,
  params.clustering_res,params.min_clusters,params.min_res,params.max_res,
  params.res_step,params.scimap_resolution,params.min_metaclusters,params.max_metaclusters)
  
	RUNQC(file_ch, params.qcscript, WRITECONFIGFILE.output.configfile)
	COLLECTBINDENSITY(params.collect_bin_density_script, RUNQC.output.bin_density.collect())
	COLLECTSIGSUM(params.collect_sigsum_script, RUNQC.output.sigsum.collect())
	
	if (!params.qc_only) {
	  if (params.run_seurat)	{
	    // Run Seurat with metaclustering
	    RUNSEURAT(params.seuratscript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile)
	    RUNMETACLUSTERS(params.seurat_metacluster_script, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile, RUNQC.output.all_markers.collect(), RUNSEURAT.output.seurat_clusters_noid.collect())
	  }  
	  
	  if (params.run_celesta) {
	    // Run CELESTA
	    RUNCELESTA(params.celestascript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, params.celesta_prior_matrix)
	  
	    if(params.run_seurat) {  	  // combine seurat and CELESTA output for comparison
  	  combined_output = RUNSEURAT.output.seurat_clusters \
        | combine(RUNCELESTA.output.celesta_classes, by:0)
  	  
  	  RUNCOMPARISON(params.seurat_vs_celesta_script, combined_output)
	    }
	  }
	
	  if (params.run_scimap) {
	    // Run Scimap and build report
	    RUNSCIMAP(params.scimapscript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile)
	    SCIMAPREPORT(params.scimap_report_script, RUNSCIMAP.output.matrixplot, RUNSCIMAP.output.spatialplot, RUNSCIMAP.output.umap)
	
	    if(params.run_seurat) {
	      // Combine seurat and scimap output for comparison
	      combined_output_seurat_scimap = RUNSEURAT.output.seurat_clusters \
          | combine(RUNSCIMAP.output.scimap_clusters, by:0)
	  
	     SEURATVSCIMAP(params.seurat_vs_scimap_script, combined_output_seurat_scimap)
	    }
	  }
	      
	}
	  
}



