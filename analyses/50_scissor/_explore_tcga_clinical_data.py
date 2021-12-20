# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import pandas as pd

# %%
df = pd.read_csv("/home/sturm/Downloads/media-2.tsv", sep="\t")

# %%
lung_filtered = df.loc[df["tumor"].isin(["LUAD", "LUSC"]), :].loc[lambda x: ~x["bcr_drug_barcode"].isnull(), :]

# %%
lung_filtered

# %%
lung_filtered["measure_of_response"].value_counts()

# %%
lung_filtered["drug_name"].value_counts().head(10)

# %%
