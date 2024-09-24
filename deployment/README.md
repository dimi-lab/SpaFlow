# Deploy spaFlow on GCP using Google Batch

Run spaFlow on Google Cloud using Google Batch service.  

## Perform Setup and Test Steps

Run a test Nextflow job to verify your service account permissions for Google Batch and GCP resources.  This will ensure your `nextflow.config` file is setup correctly to use Google Batch as a backend provider.  The test job will run a simple `hello` process to verify the setup.  


1. Clone this repo to a local directory eg on a GCE instance (VM).
2. Install Nextflow. Follow [Google's instructions to get version 23.04.1](https://cloud.google.com/batch/docs/nextflow#before-you-begin). Put the resulting `nextflow` binary on your path, or, in the repository root directory.
3. Verify Logging, Cloud Storage, Google Batch and GCE APIs are enabled on your GCP Project.
4. Make sure you have a bucket created for Nextflow job files.
5. Run a test job (hello NF) to verify your [GCP service account permissions](https://cloud.google.com/batch/docs/nextflow) and `nextflow.config` file is setup correctly to use Google Batch as a backend provider.

   1. You need at least these for your default GCE service account: logging, cloud storage, cloud build, batch.

6. Follow the instructions to set up the spaFlow param files - [link](https://github.com/dimi-lab/SpaFlow?tab=readme-ov-file#instructions)
7. Setup your `nextflow.config` and `main.nf` files to use the appropriate Google Batch resources (see below)
8. Run spaFlow using the provided command `./nextflow run main.nf -profile gcb` from the spaFlow repository root directory. (Or `nextflow` if you put it on your path.)

## Build Docker Image with Cloud Build 

Use GCP Cloud Build to build a spaFlow container and push it to GCP Artifact Registry.

1. Create an Artifact Registry in your GCP project named `images`, note the region and project ID.

2. Auth to GCP Artifact Registry (using your region):

   `gcloud auth configure-docker us-central1-docker.pkg.dev`

3. Update cloud-build.yaml with your region & project name.

4. Submit a Cloud Build job; from the root repository directory:

   ````sh
   gcloud builds submit --config cloud-build.yaml --machine-type=E2_HIGHCPU_8
   ````

## Nextflow Config for Google Batch

Update the `nextflow.config` file section `profiles` with the correct region and project id and bucket.

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

