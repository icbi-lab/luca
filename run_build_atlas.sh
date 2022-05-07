#!/usr/bin/bash

# HACK: run in separate hidden subdirectory - this makes it possible for the
# two workflows to run at the same time.
mkdir -p .build_atlas && cd .build_atlas && \
nextflow run .. --workflow build_atlas \
    --outdir ../data/20_build_atlas \
    -resume \
    -profile icbi_lung \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/atlas
