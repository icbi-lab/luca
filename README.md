# LuCA - The single-cell **Lu**ng **C**ancer **A**tlas

The single cell lung cancer atlas is a resource integrating more than 1.2 million cells from 309 patients across 29 datasets.

The atlas is publicly available for interactive exploration through a *cell-x-gene* instance. We also provide
`h5ad` objects and a [scArches](https://scarches.readthedocs.io/en/latest/) model which allows to project custom datasets
into the atlas. For more information, check out the

 * [project website](https://luca.icbi.at) and
 * our [preprint](https://doi.org/10.1101/2022.05.09.491204).

> Salcher, Sturm, Horvath et al., Manuscript in preparation

This repository contains the source-code to reproduce the single-cell data analysis for the paper.
The analyses are wrapped into [nextflow](https://github.com/nextflow-io/nextflow/) pipelines, all dependencies are
provided as [singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) containers, and input data are
available [from zenodo](https://doi.org/10.5281/zenodo.6411867).

For clarity, the project is split up into two separate workflows:

 * `build_atlas`: Takes one `AnnData` object with UMI counts per dataset and integrates them into an atlas.
 * `downstream_analyses`: Runs analysis tools on the annotated, integrated atlas and produces plots for the publication.

The `build_atlas` step requires specific hardware (CPU + GPU) for exact reproducibility
(see [notes on reproducibility](#notes-on-reproducibility)) and is relatively computationally
expensive. Therefore, the `downstream_analysis` step can also operate on pre-computed results of the `build_atlas` step,
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
curl "https://zenodo.org/record/6997383/files/containers.tar.xz?download=1" | tar xvJ

# input data
curl "https://zenodo.org/record/6997383/files/input_data.tar.xz?download=1" | tar xvJ

# OPTIONAL: obtain intermediate results if you just want to run the `downstream_analysis` workflow
curl "https://zenodo.org/record/6997383/files/build_atlas_results.tar.xz?download=1" | tar xvJ
```

Note that some steps of the downstream analysis depend on an additional [cohort of checkpoint-inhibitor-treated patients](https://ega-archive.org/studies/EGAS00001005013), which is only available under protected access agreement. For obvious reasons, these data
are not included in our data archive. You'll need to obtain the dataset yourself and place it in the `data/14_ici_treatment/Genentech` folder.
The corresponding analysis steps are skipped by default. You can enable them by adding the `--with_genentech` flag to the `nextflow run` command.

### 3. Configure nextflow

Depending on your HPC/cloud setup you will need to adjust the nextflow profile in `nextflow.config`, to tell
nextflow how to submit the jobs. Using a `withName:...` directive, special
resources may be assigned to GPU-jobs. You can get an idea by checking out the `icbi_lung` profile - which we used to run the
workflow on our on-premise cluster. Only the `build_atlas` workflow makes use of GPU processes.

### 4. Launch the workflows

```bash
# Run `build_atlas` workflow
nextflow run main.nf --workflow build_atlas -resume -profile <YOUR_PROFILE> \
    --outdir "./data/20_build_atlas"

# Run `downstream_analysis` workflow
nextflow run main.nf --workflow downstream_analyses -resume -profile <YOUR_PROFILE> \
    --build_atlas_dir "./data/20_build_atlas" \
    --outdir "./data/30_downstream_analyses"
```

As you can see, the `downstream_analysis` workflow requires the output of the `build_atlas` workflow as input.
The intermediate results from zenodo contain the output of the `build_atlas` workflow.

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
  * Integration of additional datasets with transfer learning using [scArches](scarches.readthedocs.io/).

## Downstream analysis workflow

 * Patient stratification into immune phenotypes
 * Subclustering and analysis of the neutrophil cluster
 * Differential gene expression analysis using [pseudobulk + DESeq2](https://www.nature.com/articles/s41467-021-25960-2)
 * Differential analysis of transcription factors, cancer pathways and cytokine signalling using [Dorothea](https://github.com/saezlab/dorothea-py), [progeny](https://github.com/saezlab/progeny-py), and [CytoSig](https://github.com/data2intelligence/CytoSig).
 * Copy number variation analysis using [SCEVAN](https://github.com/AntonioDeFalco/SCEVAN)
 * Cell-type composition analysis using [scCODA](https://github.com/theislab/scCODA)
 * Association of single cells with phenotypes from bulk RNA-seq datasets with [Scissor](https://github.com/sunduanchen/Scissor)
 * Cell2cell communication based on differential gene expression and the [CellphoneDB database](https://github.com/ventolab/CellphoneDB).

## Contact

For reproducibility issues or any other requests regarding single-cell data analysis, please use the [issue tracker](https://github.com/icbi-lab/luca/issues). For anything else, you can reach out to the corresponding author(s) as indicated in the manuscript.

## Notes on reproducibility

We aimed at making this workflow reproducible by providing all input data, containerizing all software
dependencies and integrating all analysis steps into a nextflow workflow.
In theory, this allows to execute the workflow on any system that can run nextflow and singularity.
Unfortunately, some single cell analysis algorithms (in particular scVI/scANVI and UMAP) will yield
slightly different results on different hardware, trading off computational reproducibility for a
significantly faster runtime. In particular, results will differ when changing the number of cores, or
when running on a CPU/GPU of a different architecture. See also https://github.com/scverse/scanpy/issues/2014 for a discussion.

Since the cell-type annotation depends on clustering, and the clustering depends on the neighborhood graph,
which again depends on the scANVI embedding, running the `build_atlas` workflow on a different machine
will likely break the cell-type labels.

Below is the hardware we used to execute the `build_atlas` workflow. Theoretically,
any CPU/CPU of the same generation shoud produce identical results, but we did not have the chance to test this yet.

 * Compute node CPU: `Intel(R) Xeon(R) CPU E5-2699A v4 @ 2.40GHz` (2x)
 * GPU node CPU: `EPYC 7352 24-Core` (2x)
 * GPU node GPU: `Nvidia Quadro RTX 8000 GPU`

