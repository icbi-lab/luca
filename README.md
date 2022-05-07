# LuCA - The single-cell **Lu**ng **C**ancer **A**tlas

The single cell lung cancer atlas is a resource integrating more than 1.2 million cells from 309 patients across 29 datasets. 
The atlas is publicly available for interactive exploration through a *cell-x-gene* instance. We also provide 
`h5ad` objects and a [scArches](https://scarches.readthedocs.io/en/latest/) model which allows to project custom datasets
onto the atlas. 

For more information, check out the project website (https://luca.icbi.at) and our preprint: 

> Salcher, Sturm, Horvath et al., Manuscript in preparation


This repository contains the source-code to reproduce the single-cell data analysis for the paper. 
The analyses are wrapped into [nextflow](https://github.com/nextflow-io/nextflow/) pipelines, all dependencies are 
provided as [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) containers, and input data are
available [from zenodo](https://doi.org/10.5281/zenodo.6411868).

For clarity, the project is split up into two separate workflows: 

 * `build_atlas`: Takes one `AnnData` object with UMI counts per dataset and integrates them into an atlas. 
 * `downstream_analyses`: Runs analysis tools on the annotated, integrated atlas and produces plots for the publication. 

The `build_atlas` step requires specific hardware (CPU + GPU) for exact reproducibility 
(see [notes on reproducibility](#notes-on-reproducibility)) and is relatively computationally 
expensive. Therefore the `downstream_analysis` step can also operate pre-computed results of the `build_atlas` step, 
which are available from zenodo. 

 ## Launching the workflows

### 1. Prerequisites

* [Nextflow](https://www.nextflow.io/index.html#GetStarted), version 21.10.6 or higher
* [Singularity/Apptainer](https://apptainer.org/), version 3.7 or higher (tested with 3.7.0-1.el7)
* A high performance cluster (HPC) or cloud setup. The whole analysis will consume several thousand CPU hours. 
 
 ### 2. Obtain data

 Before launching the workflow, you need to obtain input data and singularity containers from zenodo. 
 First of all, clone this repository:

 ```bash
git clone https://github.com/icbi-lab/luca.git
cd luca
 ```

Then, within the repository, download the data archives and extract then to the corresponding directories: 

 ```bash
 # singularity containers
curl TODO | tar xvf

# input data
curl TODO | tar xvf

# OPTIONAL: obtain intermediate results if you just want to run the `downstream_analysis` workflow
curl TODO | tar xvf
 ```

 Note, that some steps of the downstream analysis depend on an additional [cohort of checkpoint-inhibitor-treated patients](https://ega-archive.org/studies/EGAS00001005013), which is only available under protected access agreement. For obvious reasons, these data 
 are not included in our data archive. You'll need to obtain the dataset yourself and place it in the `data/14_ici_treatment/Genentech` folder. 
 The corresponding analysis steps are skipped by default. 

 ### 3. Configure nextflow

Depending on your HPC/cloud setup you likely will need to adjust the nextflow profile in `nextflow.config`, to tell 
nextflow how to submit the jobs. There is a `withLabel:gpu` directive, that can be used to assign special 
resources to gpu jobs. You can get an idea by checking out the `icbi_lung` profile - which we used to run the 
workflow on our on-premise cluster. Note that only the `build_atlas` workflow makes use of GPU processes. 

### 4. Launch the workflows

```bash
# Run `build_atlas` workflow
nextflow run main.nf --workflow build_atlas -resume -profile <YOUR_PROFILE> \
    --outdir "./data/20_build_atlas"

# Run `downstream_analysis` workflow
nextflow run main.nf --workflow downstream_analyses -resume -profile <YOUR_PROFILE> \
    --atlas "./data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad" \
    --outdir ./data/30_downstream_analyses 
```

As you can see, the `downstream_analysis` workflow requires an input file generated by the `build_atlas` workflow. 
By extracting the intermediate results from zenodo, the same directory structure will be created. 

 ## Structure of this repository

* `analyses`: Place for e.g. jupyter/rmarkdown notebooks, gropued by their respective (sub-)workflows. 
* `bin`: executable scripts called by the workflow
* `conf`: nextflow configuration files for all processes
* `containers`: place for singularity image files. Not part of the git repo and gets created by the download command. 
* `data`: place for input data and results in different subfolders. Gets populated by the download commands and by running the workflows. 
* `lib`: custom libraries and helper functions
* `modules`: nextflow DSL2.0 modules
* `preprocessing`: scripts used to preprocess data upstream of the nextflow workflows. The processed data are part of the archives on zenodo. 
* `subworkflows`: nextflow subworkflows
* `tables`: contains static content that should be under version control (e.g. manually created tables) 
* `workflows`: the main nextflow workflows


## Build atlas workflow

The `build_atlas` workflow comprises the following steps: 
  * QC of the individual datasets based on detected genes, read counts and mitochondrial fractions
  * Merging of all datasets into a single `AnnData` object. Harmonization of gene symbols. 
  * Annotation of two "seed" datasets as input for [scANVI](https://scarches.readthedocs.io/en/latest/scanvi_surgery_pipeline.html).
  * Integration of datasets with scANVI
  * Doublet removal with [SOLO](https://docs.scvi-tools.org/en/stable/api/reference/scvi.external.SOLO.html)
  * Annotation of cell-types based on marker genes and unsupervised [leiden](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html) clustering.
TODO update

## Downstream analysis workflow

TODO 

## Notes on reproducibility

For the workflow to be reproducible (consistent neighborhood graph -> consistent leiden clusters -> consistent annotation), you need specific hardware: 
   * Nvidia Quadro RTX 8000 GPU (any Nvidia GPU of the same generation *should* work)
   * Intel(R) Xeon(R) CPU E5-2699A v4 @ 2.40GHz (any Intel CPU of the same generation *should* work)

Note that changing the number of cores per process will break reproducibility. Therefore, you'll need a CPU with at least 44 cores. 

TODO link github issues



