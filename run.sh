#!/bin/bash

~/.nextflow/nextflow run . --workflow build_atlas -resume --samplesheet example/sampleSheet.csv --samplesheet2 example/sampleSheet2.csv --outdir out