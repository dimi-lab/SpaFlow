#!/usr/bin/env nextflow

params.configfile = "${projectDir}/configs.csv"
params.markerconfigfile = "${projectDir}/marker_configs.csv"
params.qcscript = "${projectDir}/scripts/QC.Rmd"
params.collect_bin_density_script= "${projectDir}/scripts/collect_bin_density.Rmd"
params.collect_sigsum_script= "${projectDir}/scripts/collect_sigsum.Rmd"
params.clusterscript = "${projectDir}/scripts/clustering.Rmd"
params.metascript = "${projectDir}/scripts/metaclustering.Rmd"

params.input_dir = "${projectDir}/data"
params.data_pattern = "${params.input_dir}/*.[tc]sv"
params.output_dir = "${projectDir}"

process RUNQC {
  maxForks 10
  publishDir(
          path: "${params.output_dir}/output_tables"
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
        path: "${params.output_dir}/output_reports",
        pattern: "*.html"
  )
  
  input:
  path collect_bin_density_script
  path bin_density_collected
  
  output:
  path "bin_density_report.html"
  
  script:
  """
  Rscript -e "rmarkdown::render('${collect_bin_density_script}', 
                                  output_file='bin_density_report.html', 
                                params = list(bin_density_collected='${bin_density_collected.join(",")}'))"
  """  
}

process COLLECTSIGSUM {
  publishDir(
        path: "${params.output_dir}/output_reports",
        pattern: "*.html"
  )
  
  input:
  path collect_sigsum_script
  path sigsum_collected
  
  output:
  path "sigsum_report.html"
  
  script:
  """
  Rscript -e "rmarkdown::render('${collect_sigsum_script}', 
                                  output_file='sigsum_report.html', 
                                params = list(sigsum_collected='${sigsum_collected.join(",")}'))"
  """  
}

process RUNCLUSTERS {
  maxForks 10

  publishDir(
        path: "${params.output_dir}/output_reports",
        pattern: "*.html"
  )
    
  publishDir(
        path: "${params.output_dir}/output_tables",
        pattern: "*.csv"
  )
  
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

  publishDir(
        path: "${params.output_dir}/output_reports",
        pattern: "*.html"
  )
    
  publishDir(
        path: "${params.output_dir}/output_tables",
        pattern: "*.csv"
  )
  
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
                                output_file='metacluster_report.html')"
  """
  
}

params.qc_only = false
params.qc_and_cluster=false

workflow {
  file_ch = Channel.fromPath(params.data_pattern)
  
	RUNQC(file_ch, params.qcscript, params.configfile)
	COLLECTBINDENSITY(params.collect_bin_density_script, RUNQC.output.bin_density.collect())
	COLLECTSIGSUM(params.collect_sigsum_script, RUNQC.output.sigsum.collect())
	
	if (params.qc_and_cluster || !params.qc_only) {
	  RUNCLUSTERS(params.clusterscript, RUNQC.output.all_markers, params.configfile, params.markerconfigfile)
	}
	
	if (!params.qc_and_cluster & !params.qc_only) {
	  RUNMETACLUSTERS(params.metascript, params.configfile, params.markerconfigfile, RUNQC.output.all_markers.collect(), RUNCLUSTERS.output.clusters.collect())
	}	
}









