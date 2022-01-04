#!/usr/bin/env python

import os
import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import sys

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor="white")
sc.settings.set_figure_params(figsize=(5, 5))

arg1 = sys.argv[1]
adata = sc.read_h5ad(arg1)

# TODO
cnv.io.genomic_position_from_gtf(
    "/data/genomes/hg38/annotation/gencode/gencode.v33.primary_assembly.annotation.gtf",
    adata,
)


immune_cells = [
    "B cell",
    "B cell dividing",
    "cDC1",
    "cDC2",
    "DC mature",
    "Granulocytes",
    "Macrophage",
    "Macrophage FABP4+",
    "Mast cell",
    "Monocyte",
    "NK cell",
    "pDC",
    "T cell CD4",
    "T cell CD8",
    "T cell dividing",
    "other (T assoc.)",
    "T cell regulatory",
]

immune_cells_in_set = set(immune_cells).intersection(adata.obs["cell_type"])
immune_cells_in_set = list(immune_cells_in_set)

cnv.tl.infercnv(
    adata,
    reference_key="cell_type",
    reference_cat=immune_cells_in_set,
    window_size=100,
    step=1,
)

pcorr = np.corrcoef(
    adata[adata.obs["cell_type"] == "Tumor cells", :].obsm["X_cnv"].todense(),
    rowvar=False,
)
q75, q25 = np.percentile(pcorr, [75, 25])
ithcna = q75 - q25
np.savetxt("ithcna.txt", [ithcna])

cnv.pl.chromosome_heatmap(adata, groupby="cell_type", save="cnv_cell_type.png")
cnv.pl.chromosome_heatmap(adata, groupby="origin", save="cnv_origin.png")

adata.write_h5ad(arg1.replace(".h5ad", ".infercnvpy.h5ad"))
