#!/usr/bin/env nextflow

params.configfile = "${projectDir}/configs.csv"
params.markerconfigfile = "${projectDir}/marker_configs.csv"
params.qcscript = "${projectDir}/scripts/QC.Rmd"
params.clusterscript = "${projectDir}/scripts/clustering.Rmd"
params.metascript = "${projectDir}/scripts/metaclustering.Rmd"
params.data_pattern = "${projectDir}/data/*.[tc]sv"

process RUNQC {
  publishDir = 'output'

  input:
  path quantfile
  path QC_scriptpath
  path configs

	output:
	path "QC_report_${roi}.html"
  path "all_markers_clean_${roi}.csv", emit: all_markers
	
	script:
	roi = quantfile.baseName
	
	"""
	Rscript -e "rmarkdown::render('${QC_scriptpath}', output_file='QC_report_${quantfile.baseName}.html')"  $quantfile $configs
	"""
}

process RUNCLUSTERS {
  publishDir = 'output'

  input:
  path clustering_script
  path all_markers
  path configs
  path marker_configs
  
  output:
  path "clustering_report_${roi}.html"
  path "clusters_${roi}.csv", emit: clusters
    
  script:
  roi = all_markers.baseName.replace("all_markers_clean_", "")
  
  """
  Rscript -e "rmarkdown::render('${clustering_script}', output_file='clustering_report_${roi}.html')" $all_markers $configs $marker_configs
  """
}

process RUNMETACLUSTERS {

  publishDir = 'output'
  
  input:
  path metascript
  path configs
  path marker_configs
  path allmarkers_collected
  path clusters_collected
  
  output:
  path "metacluster_report.html"
  path "*_metaclusters.csv"
  
  script:
  """
  Rscript -e "rmarkdown::render('${metascript}', 
                                output_file='metacluster_report.html', 
                                params = list(allmarkers_collected='${allmarkers_collected.join(",")}',
                                              clusters_collected='${clusters_collected.join(",")}',
                                              configs_path='${configs}',
                                              marker_configs_path = '${marker_configs}'))"
  """
  
}

params.qc_only = false
params.qc_and_cluster=false

workflow {
  file_ch = Channel.fromPath(params.data_pattern)
  
	RUNQC(file_ch, params.qcscript, params.configfile)
	
	if (params.qc_and_cluster || !params.qc_only) {
	  RUNCLUSTERS(params.clusterscript, RUNQC.output.all_markers, params.configfile, params.markerconfigfile)
	}
	
	if (!params.qc_and_cluster & !params.qc_only) {
	  RUNMETACLUSTERS(params.metascript, params.configfile, params.markerconfigfile, RUNQC.output.all_markers.collect(), RUNCLUSTERS.output.clusters.collect())
	}	
}









