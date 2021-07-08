#!/usr/bin/bash

cd ../.. && \
nextflow run ./single-cell-analysis-nf/main.nf \
    --input tables/samplesheet_scrnaseq_preprocessing.csv \
    --outdir data/20_qc_norm_scrnaseq \
    -resume \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung \
    --publish_dir_mode=link \
    --skip_solo
