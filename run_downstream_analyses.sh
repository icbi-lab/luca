#!/usr/bin/bash

nextflow run main.nf --workflow downstream_analyses \
    --additional_input ./tables/samplesheet_scrnaseq_preprocessing2.csv \
    --atlas /data/projects/2020/Pircher-scRNAseq-lung/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad \
    --reference_scanvi_h5ad /data/projects/2020/Pircher-scRNAseq-lung/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad \
    --reference_scanvi_model /data/projects/2020/Pircher-scRNAseq-lung/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_scanvi_model \
    --outdir ./data/30_downstream_analyses \
    -resume \
    -profile icbi \
    --publish_dir_mode symlink \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/downstream
