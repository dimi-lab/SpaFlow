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


// Load modules
include { WRITECONFIGFILE; WRITEMARKERFILE } from './modules/prep_config.nf'

include { RUNQC; COLLECTBINDENSITY; COLLECTSIGSUM } from './modules/qc.nf'

include { RUNSEURAT; RUNCELESTA; RUNSCIMAP; SCIMAPREPORT } from './modules/clustering.nf'

include {RUNMETACLUSTERS; SEURATVCELESTA; SEURATVSCIMAP } from './modules/post_clustering.nf'


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
  	  
  	  SEURATVCELESTA(params.seurat_vs_celesta_script, combined_output)
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



