# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %%
import scvelo as scv
import scanpy as sc
from multiprocessing import Pool
import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from nxfvars import nxfvars
from glob import glob
from tqdm.contrib.concurrent import process_map

rcParams["axes.grid"] = False

# %%
adata_n_path = nxfvars.get(
    "adata_n_path",
    "/home/sturm/Downloads/adata_neutrophil_clusters.h5ad"
    # "../../data/30_downstream_analyses/04_neutrophil_subclustering/artifacts/adata_neutrophil_clusters.h5ad",
)
velocyto_dir = nxfvars.get("velocyto_dir", "../../data/11_own_datasets/velocyto/")
artifact_dir = nxfvars.get("artifact_dir", "/home/sturm/Downloads")
cpus = nxfvars.get("cpus", 8)

# %%
adata_n = sc.read_h5ad(adata_n_path)

# %%
filename_map = {
    "UKIM-V_P1": "Combined_HJNHMDRXX_111723_GCTACGCT_S1_Lanes_final",
    "UKIM-V_P2": "Combined_nsclc242_S1_Lanes_final",
    "UKIM-V_P3": "Combined_124223_final",
    "UKIM-V-2_P4": "Combined_182807_final",
    "UKIM-V-2_P5": "_1_Combined_173026_final",
    "UKIM-V-2_P6": "Combined_182808_final",
    "UKIM-V-2_P7": "Combined_182809_final",
    "UKIM-V-2_P8": "Combined_182810_final",
    "UKIM-V-2_P9": "Combined_200144",
    "UKIM-V-2_P10": "Combined_200146",
    "UKIM-V-2_P11": "Combined_202833",
    "UKIM-V-2_P12": "Combined_202834",
    "UKIM-V-2_P13": "Combined_202831",
    "UKIM-V-2_P14": "Combined_202832",
    "UKIM-V-2_P15": "Combined_202829",
    "UKIM-V-2_P16": "Combined_202830",
    "UKIM-V-2_P17": "Combined_202835",
}


# %%
def _read_scvelo(patient, filename):
    path = list(glob(f"{velocyto_dir}/{filename}_*.loom"))
    assert len(path) == 1
    adata = scv.read_loom(path[0])
    adata.obs["patient"] = patient
    adata.obs["filename"] = filename
    adata.obs_names = [patient + ":" + x.split(":")[1] for x in adata.obs_names]
    return adata


adatas = process_map(
    _read_scvelo, filename_map.keys(), filename_map.values(), max_workers=cpus
)

# %%
adata_ukim = adata_n[adata_n.obs["dataset"].str.contains("UKIM"), :].copy()
adata_ukim.obs_names = [
    f"{patient}:{bc.split('_')[0]}"
    for patient, bc in zip(adata_ukim.obs["patient"], adata_ukim.obs_names)
]

# %%
for ad in adatas:
    ad.var_names_make_unique()

# %%
adata_scvelo = anndata.concat(adatas)

# %%
assert not any(adata_scvelo.obs_names.duplicated())

# %%
# Number of common cells in transcriptomics and scvelo anndata objects.
len(set(adata_scvelo.obs_names)), len(set(adata_ukim.obs_names)), len(
    set(adata_scvelo.obs_names) & set(adata_ukim.obs_names)
)

# %% [markdown]
# ### Preprocess anndata scvelo object

# %%
adata_scvelo = scv.utils.merge(adata_scvelo, adata_ukim)

# %%
adata_ukim.obs["patient"].value_counts()

# %%
adata_scvelo.obs["patient"].value_counts()

# %%
adata_scvelo.shape

# %%
scv.pp.filter_and_normalize(adata_scvelo)

# %%
scv.pp.moments(adata_scvelo, n_pcs=30, n_neighbors=30)

# %%
scv.tl.velocity(adata_scvelo)

# %%
scv.tl.velocity_graph(adata_scvelo, n_jobs=cpus)

# %% [markdown]
# ## Plots

# %%
scv.pl.proportions(adata_scvelo, groupby="patient")

# %%
sc.set_figure_params(figsize=(10, 10))
rcParams["axes.grid"] = False

# %% [markdown]
# ### Only UKIM-V dataset

# %%
ax = sc.pl.embedding(
    adata_ukim,
    basis="umap",
    color="cell_type",
    show=False,
    legend_loc="None",
    size=70,
    alpha=0.4,
)
scv.pl.velocity_embedding_stream(
    adata_scvelo,
    basis="umap",
    color="cell_type",
    #     arrow_color="white",
    legend_loc="right margin",
    ax=ax,
    alpha=0,
    arrow_size=2,
)

# %% [markdown]
# ### All neutrophils, with velocity from UKIM-V dataset

# %%
ax = sc.pl.embedding(
    adata_n,
    basis="umap",
    color="cell_type",
    show=False,
    legend_loc="on data",
    legend_fontoutline=4,
    legend_fontsize=30,
    size=70,
    alpha=0.3,
)
scv.pl.velocity_embedding_stream(
    adata_scvelo,
    basis="umap",
    color="cell_type",
    #     arrow_color="white",
    legend_loc="right margin",
    ax=ax,
    alpha=0,
    arrow_size=2,
)
ax.get_figure().savefig(
    f"{artifact_dir}/umap_scvelo.svg", bbox_inches="tight", dpi=1200
)

# %%
scv.tl.paga(adata_scvelo, groups="cell_type", minimum_spanning_tree=False)

# %%
fig, ax = plt.subplots(figsize=(3, 3), dpi=150)
scv.pl.paga(
    adata_scvelo,
    layout="fr",
    figsize=(3, 3),
    dpi=150,
    dashed_edges=None,
    init_pos="umap",
    ax=ax,
)
fig.savefig(f"{artifact_dir}/velocyto_paga_graph.pdf", bbox_inches="tight")

# %%
df = scv.get_df(adata_scvelo, "paga/transitions_confidence", precision=2).T
df.index = df.columns = sorted(adata_scvelo.obs["cell_type"].unique())
df.style.background_gradient(cmap="Blues").format("{:.2g}")

# %%
ax = scv.pl.paga(
    adata_scvelo,
    basis="umap",
    size=100,
    alpha=0.2,
    legend_loc="right margin",
    min_edge_width=2,
    node_size_scale=7,
    dashed_edges=None,
    show=False,
    arrowsize=30,
)
ax.get_figure().savefig(
    f"{artifact_dir}/umap_paga_graph.pdf", bbox_inches="tight", dpi=1200
)

# %%

# %%
