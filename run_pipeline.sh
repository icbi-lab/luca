#!/usr/bin/bash

nextflow run integrate_single_cell.nf \
    -resume \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung \
    -c integrate_single_cell.config
