#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'

process SOLO {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // conda "/home/sturm/.conda/envs/pircher-sc-integrate2"
    // container "containers/sc-integrate2.sif"
    container = "containers/sc-integrate2_2021-11-16.sif"

    cpus 4

    input:
        tuple val(id), path(adata)
        tuple val(id), path(scvi_model)
        each batch

    output:
        path("solo*.csv"), emit: doublets

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    # from threadpoolctl import threadpool_limits
    # threadpool_limits(${task.cpus})

    import scvi

    def set_all_seeds(seed=0):
        import os
        import random
        import numpy as np
        import torch

        scvi.settings.seed = seed
        os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
        os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
        np.random.seed(seed)  # Numpy random
        random.seed(seed)  # Python random

        torch.manual_seed(seed)
        torch.use_deterministic_algorithms(True)
        if torch.cuda.is_available():
            torch.cuda.manual_seed(seed)
            torch.cuda.manual_seed_all(seed)  # For multiGPU

    set_all_seeds()

    adata = sc.read_h5ad("${adata}")

    scvi.data.setup_anndata(adata, batch_key="batch")
    scvi_model = scvi.model.SCVI.load("${scvi_model}", adata=adata)
    solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="${batch}")
    solo.train()
    res = solo.predict()
    res["label"] = solo.predict(False)

    res.to_csv("solo_${batch}.csv")
    """
}

