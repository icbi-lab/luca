# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: SSH apollo-15 pircher-sc-integrate2
#     language: ''
#     name: rik_ssh_apollo_15_pircherscintegrate2
# ---

# %%
import scanpy_helpers as sh

# %%
cytosig = sh.compare_groups.compute_scores._cytosig_model()

# %%
cytosig

# %%
cytosig_df = cytosig.melt(ignore_index=False).loc[lambda x: x["value"] != 0].reset_index()

# %%
cytosig_df.loc[lambda x: x["value"] > 1].to_csv("/home/sturm/Downloads/cytosig_genes.csv")

# %%
cytosig.to_csv("/home/sturm/Downloads/cytosig_matrix.csv")
