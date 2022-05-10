#!/bin/bash

# # singularity containers needed to run all processes in the workflow
tar cvf - --dereference containers | xz -T 32 -9 > data/zenodo/containers.tar.xz

# # input data required to run all workflows from scratch
tar cvf - --dereference data/10_references data/13_tcga/for_scissor data/12_input_adatas data/11_own_datasets/{batch1_3patients,batch2_5patients,velocyto} | xz -T 32 -9 > data/zenodo/input_data.tar.xz

# intermediate results (required to run the second part of the workflow that doesn't depend on specialized hardware)
tar cvf - --dereference data/20_build_atlas | xz -T 32 -9 > data/zenodo/build_atlas_results.tar.xz

# final results (all results produced by the second part of the workflow)
tar cvf - --dereference data/30_downstream_analyses | xz -T 32 -9 > data/zenodo/downstream_analyses_results.tar.xz

# copy the core atlas, extended atlas and scArches model as separate files
tar cvf - --dereference  \
  data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_scanvi_model \
  data/20_build_atlas/annotate_datasets/35_final_atlas/full_atlas_hvg_integrated_scvi_integrated_scanvi.h5ad | \
  pigz -p 32 > data/zenodo/core_atlas_scanvi_model.tar.gz

pigz -c -p 32 data/20_build_atlas/annotate_datasets/35_final_atlas/artifacts/full_atlas_annotated.h5ad > data/zenodo/core_atlas.h5ad.gz

pigz -c -p 32 data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad > data/zenodo/extended_atlas.h5ad.gz


