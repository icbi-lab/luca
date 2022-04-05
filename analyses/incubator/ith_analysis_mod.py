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
import os
import numpy as np
import pandas as pd
import scanpy as sc
import infercnvpy as cnv
import matplotlib.pyplot as plt
import anndata
import seaborn as sns
import statsmodels.api as sm
import scipy as sci
from scipy.stats import pearsonr
from IPython.display import display

# %%
sns.set_theme(style="darkgrid")
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# %%
pat_strat = pd.read_csv("/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/stratify_patients/artifacts/patient_stratification.csv")

# %%
pat_strat["patient"]=pat_strat.patient.str.lower().replace()

# %%
ithcna_list = list()
ithgex_list = list()
for i in range(0, len(pat_strat["patient"])):
    ithcna = pd.read_csv("/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/infercnv/infercnvpy/full_atlas_merged_"+pat_strat["patient"][i]+"/ithcna.txt")
    # ithcna.rename(columns={"Cell type": "Cell_type"}, inplace=True)
    ithgex = pd.read_csv("/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/infercnv/infercnvpy/full_atlas_merged_"+pat_strat["patient"][i]+"/ithgex.txt")
    # ithgex.rename(columns={"Cell type": "Cell_type"}, inplace=True)
    ithcna_list.append(ithcna.columns.tolist())
    ithgex_list.append(ithgex.columns.tolist())

ithcna = pd.DataFrame(ithcna_list)
ithcna.rename(columns = {0:'ITHCNA'}, inplace = True)
ithcna = ithcna.replace({'nan': 0})
ithcna['ITHCNA'] = pd.to_numeric(ithcna['ITHCNA'])

ithgex = pd.DataFrame(ithgex_list)
ithgex.rename(columns = {0:'ITHGEX'}, inplace = True)
ithgex = ithgex.replace({'nan': 0})
ithgex['ITHGEX'] = pd.to_numeric(ithgex['ITHGEX'])

# %%
pat_strat["immune_infiltration"].value_counts()

# %%
ithcna['immune_infiltration'] = pat_strat['immune_infiltration'].values
ithcna['patient'] = pat_strat['patient'].values
ithcna['index'] = pat_strat['Unnamed: 0'].values
ithcna = ithcna.set_index('index')

ithgex['immune_infiltration'] = pat_strat['immune_infiltration'].values
ithgex['patient'] = pat_strat['patient'].values
ithgex['index'] = pat_strat['Unnamed: 0'].values
ithgex = ithgex.set_index('index')


# %%
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(15, 5))

sns.scatterplot(data=ithcna, x=ithcna.index, y="ITHCNA", hue="immune_infiltration",legend=True, ax=ax1)
sns.regplot(x=ithcna.index, y="ITHCNA", data=ithcna, ax=ax2)

data1 = ithcna.index
data2 = ithcna["ITHCNA"]
corr, _ = pearsonr(data1, data2)
print("Pearson = %s" % corr)

# %%
ithgex = ithgex[(ithgex != 0).all(1)]
ithcna = ithcna[(ithcna != 0).all(1)]

# %%
ithcna["immune_infiltration"].value_counts()

# %%
ithgex["immune_infiltration"].value_counts()

# %%
sns.boxplot(data=ithgex, y="ITHGEX", x="immune_infiltration")

# %%
sns.boxplot(data=ithcna, y="ITHCNA", x="immune_infiltration")

# %%
ithcna_w = ithcna[(ithcna["ITHCNA"]>0)]

# %%
ithcna_w.to_csv("ithcna.csv", index=False)

# %%
ithgex_w = ithgex[(ithgex["ITHGEX"]>0)]

# %%
ithgex_w.to_csv("ithgex.csv", index=False)

# %%
ax = sns.boxplot(x="immune_infiltration", y="ITHCNA", data=ithcna, whis=np.inf)
ax = sns.stripplot(x="immune_infiltration", y="ITHCNA", data=ithcna,color=".3")

# %%
ax = sns.boxplot(x="immune_infiltration", y="ITHGEX", data=ithgex, whis=np.inf)
ax = sns.stripplot(x="immune_infiltration", y="ITHGEX", data=ithgex,color=".3")

# %%
sm.stats.diagnostic.kstest_normal(x=ithgex["immune_infiltration"]=="B", dist='norm', pvalmethod='table')

# %%
adata = anndata.read_h5ad(
    "/data/projects/2020/Pircher-scRNAseq-lung/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad"
)

# %%
# only on primary tumor samples;
# exclude datasets with only a single cell-type
frac_by_condition = (
    adata.obs.loc[
        lambda x: (x["origin"] == "tumor_primary")
        & ~x["dataset"].isin(["Guo_Zhang_2018"])
    ]
    .groupby(["dataset", "patient"])
    .apply(lambda x: x.value_counts("cell_type", normalize=True))
    .reset_index(name="n_cells")
)

# %%
frac_by_condition["patient"] = frac_by_condition["patient"].str.lower().replace()
result = pd.merge(ithcna, frac_by_condition, how='inner')


# %%
cell_types = np.unique(result['cell_type'])
for cell_t in cell_types:
    plt.figure()
    cell_type = result.loc[result['cell_type'] == "%s" % cell_t]
    data1 = cell_type["ITHCNA"]
    data2 = cell_type["n_cells"]
    corr, _ = pearsonr(data2, data1)
    sns.regplot(x="ITHCNA", y="n_cells", data=cell_type).set(title="%s\nPearson = %s" % (cell_t, corr))

plt.show()

# %%
