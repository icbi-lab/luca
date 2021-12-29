# Lung cancer single cell atlas

This repository contains the source code to reproduce the analyses from 

> t.b.a

The analyses are wrapped into [nextflow](https://github.com/nextflow-io/nextflow/) pipelines. All dependencies are provided as [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) containers. Note that some computation steps have specialized hardware requirements (see below). 

For clarity, the project is split up into two separate workflows: 

 * `build_atlas`: Takes one `AnnData` object with UMI counts per dataset and integrates them into an atlas. 
 * `downstream_analyses`: Takes the integrated and annotated atlas and 

 To avoid the computationally expensive integration step, the second workflow can also operate on the integrated version downloaded from zenodo.

 ## Launching the workflows 

 The two workflows can be launched using the `run_build_atlas.sh` and `run_downstream_analyses.sh` scripts, respectively. You may need to adapt the
 nextflow calls according to your computational environment. 

 ## Structure of this repository

* `analyses`: Place for e.g. jupyter/rmarkdown notebooks, gropued by their respective (sub-)workflows. 
* `bin`: executable scripts called by the workflow
* `conf`: nextflow configuration files for all processes
* `lib`: custom libraries and helper functions
* `modules`: nextflow DSL2.0 modules
* `subworkflows`: nextflow subworkflows
* `tables`: contains static content that should be under version control (e.g. manually created tables) 
* `workflows`: the main nextflow workflows

## Data availability

 The input data for both workflows are available from zenodo 
 **TODO**

## Build atlas workflow

The `build_atlas` workflow comprises the following steps: 
  * QC of the individual datasets based on detected genes, read counts and mitochondrial fractions
  * Merging of all datasets into a single `AnnData` object. Harmonization of gene symbols. 
  * Annotation of two "seed" datasets as input for [scANVI](https://scarches.readthedocs.io/en/latest/scanvi_surgery_pipeline.html).
  * Integration of datasets with scANVI
  * Doublet removal with [SOLO](https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.SOLO.html)
  * Annotation of cell-types based on marker genes and unsupervised [leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html) clustering. 

### Notes on reproducibility

For the workflow to be reproducible (consistent neighborhood graph -> consistent leiden clusters -> consistent annotation), you need specific hardware: 
   * Nvidia Quadro RTX 8000 GPU (any Nvidia GPU of the same generation *should* work)
   * Intel(R) Xeon(R) CPU E5-2699A v4 @ 2.40GHz (any Intel CPU of the same generation *should* work)

Note that changing the number of cores per process will break reproducibility. Therefore, you'll need a CPU with at least 44 cores. 

 ## Downstream analysis workflow


