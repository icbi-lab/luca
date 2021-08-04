#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { initOptions; saveFiles; getSoftwareName } from './functions'

process SOLO {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda "/home/sturm/.conda/envs/pircher-sc-integrate2"

    cpus 4

    input:
        path adata
        path scvi_model
        each batch

    output:
        path("solo*.csv"), emit: doublets

    script:
    """
    #!/usr/bin/env python
    import scanpy as sc
    import scvi

    adata = sc.read_h5ad("${adata}")
    scvi_model = scvi.model.SCVI.load("${scvi_model}", adata=adata)
    solo = scvi.external.SOLO.from_scvi_model(scvi_model, restrict_to_batch="${batch}")
    solo.train()
    res = solo.predict()
    res["label"] = solo.predict(False)

    res.to_csv("solo_${batch}.csv")
    """
}

