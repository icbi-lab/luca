#!/usr/bin/bash

nextflow run main.nf --workflow downstream_analyses \
    --atlas "./data/30_downstream_analyses/03_update_annotation/artifacts/full_atlas_merged.h5ad" \ # TODO
    --outdir ./data/30_downstream_analyses \
    -resume \
    -profile icbi \
    --publish_dir_mode symlink \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/downstream
