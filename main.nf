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
params.leiden_metacluster_script = "${projectDir}/scripts/leiden_metaclustering.Rmd"
params.kmeans_metacluster_script = "${projectDir}/scripts/kmeans_metaclustering.Rmd"
params.seurat_vs_celesta_script = "${projectDir}/scripts/seurat_vs_celesta.Rmd"
params.seurat_vs_scimap_script = "${projectDir}/scripts/seurat_vs_scimap.Rmd"


// Load modules
include { WRITECONFIGFILE; WRITEMARKERFILE } from './modules/prep_config.nf'

include { RUNQC; COLLECTBINDENSITY; COLLECTSIGSUM } from './modules/qc.nf'

include { RUNSEURAT; RUNCELESTA; RUNSCIMAP; SCIMAPREPORT } from './modules/clustering.nf'

include {RUNMETACLUSTERSSEURAT; RUNMETACLUSTERSLEIDEN; RUNMETACLUSTERSKMEANS; SEURATVCELESTA; SEURATVSCIMAP; RUNSOMCLUSTERS } from './modules/post_clustering.nf'


workflow {
  file_ch = Channel.fromPath(params.data_pattern)
  
  WRITEMARKERFILE(params.markers)
  WRITECONFIGFILE(params.sigsum_quantile_high,params.sigsum_quantile_low,
  params.bin_size,params.density_cutoff,params.cluster_metric,
  params.clustering_res,params.min_clusters,params.min_res,params.max_res,
  params.res_step,params.scimap_resolution,params.min_metaclusters,params.max_metaclusters)
  
	RUNQC(file_ch, params.qcscript, WRITECONFIGFILE.output.configfile, params.filter_column)
	COLLECTBINDENSITY(params.collect_bin_density_script, RUNQC.output.bin_density.collect())
	COLLECTSIGSUM(params.collect_sigsum_script, RUNQC.output.sigsum.collect())
	
	// Print a message if there are NAs in the dataset
	RUNQC.output.NACHECK
    .branch { v ->
        nacheck_ok: v == "0"
        nacheck_warn: v == "1"
    }
    .set { result }

  result.nacheck_warn.view { v -> "***WARNING: Rows containing NaN values have been removed from your dataset. Please double-check your quantification files.***" }

	
	if (!params.qc_only) {
	  if (params.run_seurat)	{
	    // Run Seurat with metaclustering
	    RUNSEURAT(params.seuratscript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile)
	    RUNMETACLUSTERSSEURAT(params.seurat_metacluster_script, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile, RUNQC.output.all_markers.collect(), RUNSEURAT.output.seurat_clusters_noid.collect())

		if(params.run_som)  {
                  somSplitList = Channel.from(params.min_som_clusters..params.max_som_clusters)
		  RUNSOMCLUSTERS(params.som_clustering_script, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile, RUNQC.output.all_markers.collect(), RUNSEURAT.output.seurat_centroids.collect(), RUNSEURAT.output.seurat_clusters_noid.collect(), somSplitList)
		}
	  
	  if (params.run_celesta) {
	    // Run CELESTA
	    RUNCELESTA(params.celestascript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, params.celesta_prior_matrix)
	  
	    if(params.run_seurat) {  	  // combine seurat and CELESTA output for comparison
  	  combined_output = RUNSEURAT.output.seurat_clusters \
        | combine(RUNCELESTA.output.celesta_classes, by:0)
  	  
  	  SEURATVCELESTA(params.seurat_vs_celesta_script, combined_output)
	    }
	
	  if (params.run_scimap) {
	    // Run Scimap and build report
	    RUNSCIMAP(params.scimapscript, RUNQC.output.all_markers, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile)
	    SCIMAPREPORT(params.scimap_report_script, RUNSCIMAP.output.matrixplot_K, RUNSCIMAP.output.matrixplot_L, RUNSCIMAP.output.optimalplot, \
	     RUNSCIMAP.output.spatialplot_K, RUNSCIMAP.output.spatialplot_L, RUNSCIMAP.output.umap_K, RUNSCIMAP.output.umap_L, RUNSCIMAP.output.clustermeterics)
	    RUNMETACLUSTERSLEIDEN(params.leiden_metacluster_script, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile, RUNQC.output.all_markers.collect(), RUNSCIMAP.output.scimap_clusters_noid.collect())
	    RUNMETACLUSTERSKMEANS(params.kmeans_metacluster_script, WRITECONFIGFILE.output.configfile, WRITEMARKERFILE.output.markerconfigfile, RUNQC.output.all_markers.collect(), RUNSCIMAP.output.scimap_clusters_noid.collect())
	
	    if(params.run_seurat) {
	      // Combine seurat and scimap output for comparison
	      combined_output_seurat_scimap = RUNSEURAT.output.seurat_clusters \
          | combine(RUNSCIMAP.output.scimap_clusters, by:0)
	  
	     SEURATVSCIMAP(params.seurat_vs_scimap_script, combined_output_seurat_scimap)
	    }
	  }
	      
	}
	  
}



