<img src="https://github.com/dimi-lab/mxif_clustering_pipeline/blob/main/images/SpaFlow.png" width="1000"/>

Pipeline visualized below

<img src="https://github.com/dimi-lab/SpaFlow/blob/main/images/Spaflow_metro.png" width="1000"/>

## Requirements/Dependencies

-   Nextflow 23.04.2 (requires: bash, java 11 [or later, up to 21] git and docker)
-   pandoc 3.1.2
-   R 4.2.2
    -   knitr 1.44
    -   ggplot2 3.4.4
    -   data.table 1.14.8
    -   dplyr 1.1.3
    -   Seurat \>5.0
    -   progressr 0.14.0
    -   kableExtra 1.3.4
    -   ComplexHeatmap 2.15.4
    -   ggridges 0.5.4
    -   clustree 0.5.0
    -   pheatmap 1.0.12
    -   plyr 1.8.9
    -   pander 0.6.5
    -   CELESTA 0.0.0.9000
        -   Rmixmod 2.1.9

        -   spdep 1.3-3

        -   reshape2 1.4.4

        -   zeallot 0.1.0
-   Python 3.8
    -   scimap 1.3.2
    -   anndata 0.7.8
    -   pandas 1.5.3
    -   scanpy 1.9.6

------------------------------------------------------------------------

## Instructions

Note: This pipeline requires exported QuPath (0.4.3) measurement tables (quantification files) generated from segmented single cell MxIF images.

1.  Clone repository to your machine
2.  Place quantification files in your input directory
    i.  Files should be in the format `<fov_name>.tsv` or `<fov_name>.csv` (e.g. `region_001.tsv`)
3.  Adjust configuration values in `configs.csv` (see 'Configurable parameters')
4.  If running CELESTA, set up `celesta_prior_matrix.csv` with your markers and cell types according to the CELESTA specification (see "Inputs" here: <https://github.com/plevritis-lab/CELESTA>)
5.  Add desired markers to `marker_configs.csv` - all markers provided must be present in the dataset
    i.  If markers are not provided or if all markers provided are not present in the dataset, the quantification file will search for a default list of markers; if all default markers are not present, all markers in the dataset will be used, excluding DAPI.
6.  Set up `params.yaml` with your input and output directories, as well as the location of your configuration files.
    i.  qc_only=true will run only the QC step. It is a good idea to first run QC and check the results before moving on to clustering.

    ii. run_scimap=false or run_celesta=false will skip running scimap or CELESTA, respectively.

    iii. export_intermediates=true will create a new directory "intermediates" in your outputs directory, containing the centroids of the seurat clusters from each FOV using the CLR normalization method, and those using the arcsin/Z-score method (AKA the inputs for metaclustering)
7.  Call main pipeline script: `nextflow run main.nf -params-file=params.yaml`

------------------------------------------------------------------------

### Configurable parameters

| object               | value                                                                                                                                                                                                              |
|----------------|--------------------------------------------------------|
| sigsum_quantile_high | Upper quantile cutoff for sigsum filtering (default 0.99)                                                                                                                                                          |
| sigsum_quantile_low  | Lower quantile cutoff for sigsum filtering (default 0.05)                                                                                                                                                          |
| bin_size             | Size of bounding box for low-density cell search (default 50). Smaller bin size is more stringent (will remove more cells at the same density cutoff).                                                             |
| density_cutoff       | Cutoff number of cells defined as low-density (default 5). Higher density cutoff is more stringent (will remove more cells at the same bin size)                                                                   |
| cluster_metric       | Metric to use for Seurat clustering (default Median)                                                                                                                                                               |
| clustering_res       | If specified, this clustering resolution will be used for all ROIs and will override the clustree method. Set to NA or remove row to use clustree method. Larger clustering resolutions will create more clusters. |
| min_res              | Minimum clustering resolution to search with clustree (default 0.1)                                                                                                                                                |
| max_res              | Maximum clustering resolution to search with clustree (default 1.9)                                                                                                                                                |
| res_step             | Increment for searching clustering resolutions; functions as `by` argument in `seq()` (default 0.2)                                                                                                                |
| min_clusters         | Minimum number of clusters for per-ROI clustering in Seurat (default 6)                                                                                                                                            |
| min_metaclusters     | Starting number of metaclusters to create (default 5)                                                                                                                                                              |
| max_metaclusters     | Ending number of metaclusters to create (default 10)                                                                                                                                                               |

------------------------------------------------------------------------

## Analysis steps and outputs

**QC.Rmd**

-   `output_reports/bin_density_report.html` and `output_reports/sigsum_report.html`
    -   Reports contain one QC image for each ROI
    -   Plots for sigsum cutoffs and bin density flags
-   `output_tables/all_markers_clean_<roi>.csv`
    -   Quantification file in the same format as input files but with additional columns for QC metrics and QC flags

**seurat_clustering.Rmd**

-   `output_reports/clustering_report_<roi>.html`
    -   Clustree plot, selected resolution, marker vs cluster heatmaps and ridgeplots
-   `output_tables/seurat_clusters_<roi>.html`
    -   Clusters mapped to cell coordinates - includes artifacts

**metaclustering.Rmd**

-   `output_reports/metaclustering_report.html`
    -   Marker vs metacluster heatmaps; barplots for proportion of ROI per metacluster
-   `output_tables/<roi>_mapped_metaclusters_n_metaclusters.csv`
    -   Clusters and metaclusters mapped to cell coordinates - includes artifacts

**CELESTA_clustering.Rmd**

-   `output_reports/celesta_report_<roi>.html`

    -   CELESTA marker thresholding plots, final cell type counts and spatial plots for CELESTA classification.

-   `output_tables/CELESTA_classes_<roi>.csv`

    -   CELESTA class predictions mapped to coordinates, including results from each round of classification.

**seurat_vs_celesta.Rmd**

-   `output_reports/classification_comparison_<roi>.html`

    -   Tables and plots comparing seurat clusters to celesta class predictions.

-   `output_tables/<roi>_combined_classes.csv`

    -   Combined mapping of seurat clusters and CELESTA final round classes with coordinates.

**scimap_clustering.py and scimap_report.Rmd**

-   `output_reports/scimap_report_<roi>.csv`

    -   UMAP, spatial, and heatmap plots for scimap clustering.

-   `output_tables/scimap_clusters_<roi>.csv`

    -   Mapping of scimap leiden clusters to cell coordinates.

**seurat_vs_scimap.Rmd**

-   `output_reports/scimap_report_<roi>.html`

    -   Tables and plots comparing scimap and seurat clustering results.

------------------------------------------------------------------------

## Example Data

Within the data directory of this repository, there is a sample dataset with four quantification files generated using images from the [Multiplexed Imaging Mass Cytometry of Chemokine Milieus in Metastatic Melanoma](https://zenodo.org/records/6004986) dataset (Hoch et al. 2022). The templates for configs.csv and marker_configs.csv have been set up for this dataset, therefore you can run a test of the pipeline by cloning the repo, unzipping the files into the data directory, and running `nextflow run main.nf` from the top-level directory.

Hoch, T., Schulz, D., Eling, N., Martínez-Gómez, J., Levesque, M., & Bodenmiller, B. (2022). Multiplexed Imaging Mass Cytometry of Chemokine Milieus in Metastatic Melanoma - Raw Data. <https://doi.org/10.5281/zenodo.6004986>
