#!/usr/bin/bash

nextflow run main.nf \
    -resume \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung \
