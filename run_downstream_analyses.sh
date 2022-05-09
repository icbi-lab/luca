#!/usr/bin/bash

nextflow run main.nf --workflow downstream_analyses \
    --build_atlas_dir "./data/20_build_atlas" \
    --outdir "./data/30_downstream_analyses" \
    --with_genentech \
    -resume \
    -profile icbi_lung \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/downstream
