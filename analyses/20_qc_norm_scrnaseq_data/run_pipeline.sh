#!/usr/bin/bash

cd ../.. && \
nextflow run /home/sturm/projects/2020/single-cell-analysis-nf/main.nf \
    --input tables/samplesheet_scrnaseq_preprocessing.csv \
    --outdir data/20_qc_norm_scrnaseq \
    -resume \
    -profile icbi_long \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung

