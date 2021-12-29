#!/usr/bin/bash

nextflow run main.nf --workflow downstream_analyses \
    --atlas /data/projects/2020/Pircher-scRNAseq-lung/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad \
    --outdir ./data/30_downstream_analyses \
    -resume \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/downstream
