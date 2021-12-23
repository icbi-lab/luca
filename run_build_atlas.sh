#!/usr/bin/bash

nextflow run main.nf --workflow build_atlas \
    --input ./tables/samplesheet_scrnaseq_preprocessing.csv \
    --outdir ./data/20_build_atlas \
    -resume \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/atlas
