#!/bin/bash

nextflow run velocyto.nf \
    -c nf_velocyto.config \
    -profile icbi \
    -w /data/scratch/sturm/projects/2020/pircher-scrnaseq-lung/scvelo
