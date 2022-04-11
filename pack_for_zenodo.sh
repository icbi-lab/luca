#!/bin/bash

# # singularity containers needed to run all processes in the workflow
tar cvf - --dereference --hard-dereference containers | xz -9 -T 32 > data/zenodo/containers.tar.xz

# # input data required to run all workflows from scratch
tar cvf - --dereference --hard-dereference data/13_tcga/for_scissor data/12_input_adatas data/11_own_datasets/{batch1_3patients,batch2_5patients,velocyto} | xz -9 -T 32 > data/zenodo/input_data.tar.xz

# intermediate results (required to run the second part of the workflow that doesn't depend on specialized hardware)
# tar cvf - --dereference --hard-dereference data/20_build_atlas | pigz -p 32 > data/zenodo/build_atlas_results.tar.gz

# final results (all results produced by the second part of the workflow)
# tar cvf - --dereference --hard-dereference data/30_downstream_analyses | pigz -p 32 > data/zenodo/downstream_analyses_results.tar.gz
# tar cvf - --dereference --hard-dereference data/30_downstream_analyses | xz -9 -T 32 > data/zenodo/downstream_analyses_results.tar.xz
