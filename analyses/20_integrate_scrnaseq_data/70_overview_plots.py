# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-2020-pircher-downstream-vis]
#     language: python
#     name: conda-env-conda-2020-pircher-downstream-vis-py
# ---

# %%
import scanpy as sc
from scanpy_helpers.annotation import AnnotationHelper
from nxfvars import nxfvars
import altair as alt
import siuba as s
from siuba import _
from toolz.functoolz import pipe

# %%
alt.data_transformers.disable_max_rows()

# %%
sc.set_figure_params(figsize=(5, 5))

# %%
ah = AnnotationHelper()

# %%
main_adata = nxfvars.get(
    "main_adata",
    "/home/sturm/Downloads/adata_annotated_fine.h5ad",
)

# %%
adata_epi = sc.read_h5ad("../../data/zz_epi/adata_epithelial_cells.h5ad")

# %%
adata = sc.read_h5ad(main_adata)

# %%
sc.pl.umap(adata, color="cell_type")

# %%
sc.pl.umap(adata_epi, color="cell_type")

# %%
ah.integrate_back(adata, adata_epi)

# %%

# %%
plt_data = adata.obs.groupby(["patient", "cell_type"]).size().reset_index(name="count")

# %%
plt_data

# %%
adata.obs.columns

# %%

# %% [markdown]
# ## Plots by sample

# %%
sample_df = adata.obs.loc[
    :, ["sample", "tissue", "origin", "condition", "dataset"]
].drop_duplicates()
assert sample_df.shape[0] == adata.obs["sample"].nunique()

# %%
(
    sample_df
    >> s.group_by("dataset")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="dataset"))
)

# %%
(
    sample_df
    >> s.group_by("tissue")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="tissue"))
)

# %%
(
    sample_df
    >> s.group_by("condition")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="condition"))
)

# %%
(
    sample_df
    >> s.group_by("origin")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="origin"))
)


# %% [markdown]
# ## Plots by patient

# %%
@s.siu.symbolic_dispatch
def _group_origins(origins): 
    origins = set(origins)
    if len(origins) == 1:
        return next(iter(origins))
    elif "tumor_primary" in origins and "normal_adjacent" in origins:
        return "matched tumor/normal"
    else:
        print(origins)
        return "multiple"
    
patient_df = adata.obs >> s.group_by("patient", "dataset", "condition") >> s.summarize(origin = _group_origins(_.origin))
assert patient_df.shape[0] == adata.obs["patient"].nunique()

# %%
(
    patient_df
    >> s.group_by("origin")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="origin"))
)

# %%
(
    patient_df
    >> s.group_by("dataset")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="dataset"))
)

# %%
(
    patient_df
    >> s.group_by("condition")
    >> s.count()
    >> s.pipe(lambda _: alt.Chart(_).mark_bar().encode(x="n", color="condition"))
)

# %%
