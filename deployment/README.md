# Deploy spaFlow on GCP using Google Batch

Run spaFlow on Google Cloud using Google Batch service.  

1. Verify Logging, Cloud Storage, Google Batch and GCE APIs are enabled on your GCP Project
2. Run a test job (hello NF) to verify your [GCP service account permissions](https://cloud.google.com/batch/docs/nextflow) and `nextflow.config` file is setup correctly to use Google Batch as a backend provider. 
3. Clone this repo to a local directory on a GCE instance (VM), unzip the sample data files and build the container (see below)
4. Follow the instructions to set up the spaFlow param files - [link](https://github.com/dimi-lab/SpaFlow?tab=readme-ov-file#instructions)
5. Setup your `nextflow.config` and `main.nf` files to use the appropriate Google Batch resources (see below)
6. Run spaFlow using the provided command `nextflow run main.nf -params-file=params.yaml -work-dir gs://my-bucket/some/path`

## Build Docker Image with Cloud Build 

Use GCP Cloud Build to build a spaFlow container.  

1. Create an Artifact Registry in your GCP project, build and push a spaFlow container using Cloud Build and the `Dockerfile`
2. Create a `cloudbuild.yaml` for use in building the container using Cloud Build
3. Use `gcloud builds submit --tag us-central1-docker.pkg.dev/<GCP-project>/images/spaflow`

Example `cloudbuild.yaml` shown below.

```
    steps:
    - name: 'gcr.io/cloud-builders/docker' 
    args: ['build', '-t', 'us-central1-docker.pkg.dev/<GCP-project>/images/spaflow', '.']
    images:
    - 'us-central1-docker.pkg.dev/<GCP-project>/images/spaflow'
```

## Nextflow Config for Google Batch

Your `nextflow.config` file must have the following section to run on [Google Batch](https://www.nextflow.io/docs/latest/google.html#configuration).  

```
    process {
    executor = 'google-batch'
    container = 'your/container:latest'
    }

    google {
        project = 'your-project-id'
        location = 'us-central1'
    }
```

Your `main.nf` file should have config information at the task level.  General syntax shown below.

```
process myTask {
    cpus 8
    memory '40 GB'

    """
    your_command --here
    """
}

process anotherTask {
    machineType 'n1-highmem-8'

    """
    your_command --here
    """
}
```
