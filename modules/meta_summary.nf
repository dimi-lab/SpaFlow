process SINGLE_SAMPLE_META_RANKING {
  cpus 8
  memory '16 GB'
  publishDir(
        path: "${params.output_dir}/output_reports/meta_ranks",
        pattern: "*.pdf",
        mode: "copy"
  )
  
  input:
  path singlerankingscript
  tuple val(sample), path(anndata_file)
  
  output:
  tuple val(sample), path("ranking_${sample}.csv"), emit: ranks

  script:
  """
  python rank_individual_metaclustering.py $sample $anndata_file
  """
  
}
