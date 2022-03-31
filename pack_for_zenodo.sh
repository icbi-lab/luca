#!/bin/bash

# singularity containers needed to run all processes in the workflow
tar cvf - --dereference --hard-dereference containers | pigz -p 32 > data/zenodo/containers.tar.gz

# input data required to run all workflows from scratch
tar cvf - --dereference --hard-dereference data/13_tcga/for_scissor data/12_input_adatas data/11_own_datasets/{batch1_3patients,batch2_5patients,velocyto} | pigz -p 32 > data/zenodo/input_data.tar.gz

# intermediate results (required to run the second part of the workflow that doesn't depend on specialized hardware)

# final results (all results produced by the second part of the workflow)
