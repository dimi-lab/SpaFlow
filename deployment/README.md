# Deploy spaFlow on GCP using Google Batch

Run spaFlow on Google Cloud using Google Batch service.  

## Perform Setup and Test Steps

Run a test Nextflow job to verify your service account permissions for Google Batch and GCP resources.  This will ensure your `nextflow.config` file is setup correctly to use Google Batch as a backend provider.  The test job will run a simple `hello` process to verify the setup.  


1. Verify Logging, Cloud Storage, Google Batch and GCE APIs are enabled on your GCP Project
2. Run a test job (hello NF) to verify your [GCP service account permissions](https://cloud.google.com/batch/docs/nextflow) and `nextflow.config` file is setup correctly to use Google Batch as a backend provider. 
3. Clone this repo to a local directory on a GCE instance (VM), unzip the sample data files and build the container (see below)
4. Follow the instructions to set up the spaFlow param files - [link](https://github.com/dimi-lab/SpaFlow?tab=readme-ov-file#instructions)
5. Setup your `nextflow.config` and `main.nf` files to use the appropriate Google Batch resources (see below)
6. Run spaFlow using the provided command `../nextflow run main.nf -profile gcb` from the spaFlow directory

## Build Docker Image with Cloud Build 

Use GCP Cloud Build to build a spaFlow container and push it to GCP Artifact Registry.

1. Create an Artifact Registry in your GCP project named `images`, note the region and project ID
2. To build and push a spaFlow container using Cloud Build and the `Dockerfile` & `.dockerignore` files, change to the directory containing the Dockerfile and spaFlow code
3. Create a `cloudbuild.yaml` for use in building the container using Cloud Build
4. Auth to cloud shell and auth to GCP Artifact Registry `gcloud auth configure-docker us-central1-docker.pkg.dev`
5. Use `gcloud builds submit --tag us-central1-docker.pkg.dev/<YOUR_PROJECT>/images/spaflow --machine-type=E2_HIGHCPU_8`

Example `cloudbuild.yaml` shown below.

```
    steps:
    - name: 'gcr.io/cloud-builders/docker' 
    args: ['build', '-t', 'us-central1-docker.pkg.dev/<YOUR_PROJECT>/images/spaflow', '.']
    images:
    - 'us-central1-docker.pkg.dev/<YOUR_PROJECT>/images/spaflow'
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

NOTES for enterprise customers:  
- specify a named (project) GCP service account in the nextflow.config file
- specify a named (project) GCP network and subnet in the nextflow.config file
- specify 'no external IP' for the GCP VM instances in the nextflow.config file
- specify both CPU and MEMORY requirements for **each process** in the main.nf file to request a specific machine type
