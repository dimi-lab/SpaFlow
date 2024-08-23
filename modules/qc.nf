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