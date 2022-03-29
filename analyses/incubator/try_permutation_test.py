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
import numpy as np
import pandas as pd
import scipy.stats
import seaborn as sns
import statsmodels.formula.api as smf
from tqdm import tqdm

# %%
import numpy as np
import pandas as pd
import scipy.stats

np.random.seed(0)

# %%
lusc1 = np.random.normal(loc=4, scale=1, size=30)
luad1 = np.random.normal(loc=6, scale=1, size=3)
lusc2 = np.random.normal(loc=2, scale=1, size=3)
luad2 = np.random.normal(loc=3, scale=1, size=30)

# %% [markdown]
# `NKX2-1` is a marker for LUAD tumor cells. Let's assume it is has double the expression in LUAD than in LUSC. Let's also assume, that due to batch effects, in dataset 2 `NKX2-1`, due to batch effects, has a higher baseline expression 
#
# | dataset | # LUSC | # LUAD | NKX2-1 baseline expression  |
# | -- | -- | -- | -- |
# 1 | 30 | 3 | 4
# 2 | 3 | 30| 2
#

# %%
df = pd.DataFrame().assign(
    expr=np.hstack([lusc1, luad1, lusc2, luad2]),
    tumor_type=["LUSC"] * 30 + ["LUAD"] * 3 + ["LUSC"] * 3 + ["LUAD"] * 30,
    dataset=["batch1"] * 33 + ["batch2"] * 33,
)

# %%
np.mean(df.loc[lambda x: x["tumor_type"] == "LUAD", "expr"]) / np.mean(
    df.loc[lambda x: x["tumor_type"] == "LUSC", "expr"]
)

# %%
scipy.stats.ttest_ind(
    df.loc[lambda x: x["tumor_type"] == "LUSC", "expr"],
    df.loc[lambda x: x["tumor_type"] == "LUAD", "expr"],
)

# %%
sns.swarmplot(x="tumor_type", y="expr", data=df, hue="dataset")

# %%
model = smf.ols("expr ~ C(tumor_type, Treatment('LUSC'))", data=df)

# %%
res = model.fit()

# %%
res.summary()


# %% [markdown]
# ## Permutation test

# %%
def get_stat(df):
    return np.mean(df.loc[lambda x: x["tumor_type"] == "LUSC", "expr"]) - np.mean(
        df.loc[lambda x: x["tumor_type"] == "LUAD", "expr"]
    )


# %%
stats = []
df_rand = df.copy()
for i in tqdm(range(1000)):
    np.random.shuffle(df_rand["tumor_type"])
    stats.append(get_stat(df_rand))

# %%
stats = np.array(stats)

# %%
stats

# %%
stat = get_stat(df)

# %%
np.sum(stat > stats) / stats.size

# %%
np.sum(stat < stats) / stats.size

# %%
np.sum()

# %%
