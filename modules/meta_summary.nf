process SINGLE_SAMPLE_META_RANKING {
  cpus 8
  memory '16 GB'

  input:
  path singlerankingscript
  tuple val(sample), path(anndata_file)
  
  output:
  tuple val(sample), path("ranked_clusterings_${sample}.csv"), emit: ranking
  path("*.png"), emit: plots

  script:
  """
  python rank_individual_metaclustering.py $sample $anndata_file
  """ 
}


process SINGLE_SAMPLE_META_RANKING_REPORT {
  cpus 8
  memory '24 GB'

  publishDir(
        path: "${params.output_dir}/output_reports/meta_summary",
        pattern: "*.html",
        mode: "copy"
  )
  
  input:
  path single_meta_report_script
  tuple val(sample), path(rank_table)
  path plots

  output:
  path "*.html"
  
  """
  Rscript -e "rmarkdown::render('${single_meta_report_script}', output_file='single_meta_report_${sample}.html')"
  """  
}
