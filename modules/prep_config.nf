process WRITECONFIGFILE {
  cpus 8
  memory '24 GB'

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
  val som_grid_x
  val som_grid_y
  val min_som_clusters
  val max_som_clusters
  val globals_maxsize_MB
  
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
  echo "som_grid_x,!{som_grid_x}" >> configs.csv
  echo "som_grid_y,!{som_grid_y}" >> configs.csv
  echo "min_som_clusters,!{min_som_clusters}" >> configs.csv
  echo "max_som_clusters,!{max_som_clusters}" >> configs.csv
  echo "globals_maxsize_MB,!{globals_maxsize_MB}" >> configs.csv
  """
}

process WRITEMARKERFILE {
  cpus 8
  memory '24 GB'

  input:
  val markers
  
  output:
  path "marker_configs.csv", emit: markerconfigfile
  
  shell:
  """
  list="${markers}"

  IFS=',' read -r -a array <<< "\$list"
  
  echo "marker" > marker_configs.csv

  for item in "\${array[@]}"; do
    echo \$item >> marker_configs.csv
  done
  """
}
