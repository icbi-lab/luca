import altair as alt
import pandas as pd


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable):
    return alt.Scale(
        domain=list(getattr(COLORS, variable).keys()),
        range=list(getattr(COLORS, variable).values()),
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")


class COLORS:

    sex = {
        "male": "#80b1d3",
        "female": "#fccde5",
        "unknown": "#dddddd",
    }
    tumor_stage = {
        "early": "#998ec3",
        "late": "#f1a340",
    }
    TMIG = {  # tumor microenvironment infiltration group
        "-/-": "#CCCCCC",
        "-/S": "#999999",
        "T/S": "#1f78b4",
        "T/-": "#a6cee3",
        "M/S": "#33a02c",
        "M/-": "#b2df8a",
        "B/S": "#ff7f00",
        "B/-": "#fdbf6f",
    }
    immune_infiltration = {
        "B": "#ff7f00",
        "T": "#1f78b4",
        "M": "#33a02c",
        "-": "#999999",
    }
    infiltration_state = {
        "I/S": "#6a3d9a",
        "I/-": "#cab2d6",
        "-/S": "#999999",
        "-/-": "#CCCCCC",
    }
    tumor_type = {
        "LUAD": "#7fc97f",
        "NSCLC": "#999999",
        "LSCC": "#f0027f",
        "LUAD EMT": "#beaed4",
        "LUAD NE": "#fdc086",
        "LUAD dedifferentiated": "#ffff99",
        # TODO: Just show these as NSCLC NOS
        "PPC": "#bf5b17",
        "LCLC": "#bf5b17",
    }
    condition = {
        "healthy_control": "#7fc97f",
        "non-cancer": "#7fc97f",
        "LSCC": "#beaed4",
        "COPD": "#fdc086",
        "LUAD": "#ffff99",
        "NSCLC": "#999999",
        "NOS": "#999999",
        "other": "#DDDDDD",
        "LUAD EMT": "#386cb0",
        "LUAD NE": "#f0027f",
        "LUAD dedifferentiated": "#bf5b17",
    }
    origin = {
        "tumor_metastasis": "#a6761d",
        "metastases": "#a6761d",
        "other": "#DDDDDD",
        "normal": "#66a61e",
        "normal_adjacent": "#d95f02",
        "adj. normal": "#d95f02",
        "tumor_primary": "#e6ab02",
        "primary tumor": "#e6ab02",
        "nan": "#DDDDDD",
    }
    platform = {
        "10x": "#377eb8",
        "Smart-seq2": "#e41a1c",
        "DropSeq": "#984ea3",
        "BD-Rhapsody": "#4daf4a",
        "InDrop": "#ff7f00",
        "Singleron": "#ffff33",
    }
    cell_type_coarse = {
        "B cell": "#1f77b4",
        "Endothelial cell": "#ff7f0e",
        "Epithelial cell": "#279e68",
        "Macrophage/Monocyte": "#d62728",
        "Mast cell": "#aa40fc",
        "NK cell": "#8c564b",
        "Neutrophils": "#e377c2",
        "Plasma cell": "#b5bd61",
        "Stromal": "#17becf",
        "T cell": "#aec7e8",
        "cDC": "#ffbb78",
        "pDC": "#98df8a",
    }
    cell_type_major = {
        "Alveolar cell type 1": "#023fa5",
        "Alveolar cell type 2": "#7d87b9",
        "B cell": "#bec1d4",
        "Ciliated": "#d6bcc0",
        "Club": "#bb7784",
        "DC mature": "#8e063b",
        "Endothelial cell": "#4a6fe3",
        "Goblet": "#8595e1",
        "Macrophage": "#b5bbe3",
        "Macrophage FABP4+": "#e6afb9",
        "Mast cell": "#e07b91",
        "Monocyte": "#d33f6a",
        "NK cell": "#11c638",
        "Neutrophils": "#8dd593",
        "Plasma cell": "#c6dec7",
        "Stromal": "#ead3c6",
        "T cell CD4": "#f0b98d",
        "T cell CD8": "#ef9708",
        "T cell regulatory": "#0fcfc0",
        "Tumor cells": "#9cded6",
        "cDC1": "#d5eae7",
        "cDC2": "#f3e1eb",
        "other": "#f6c4e1",
        "pDC": "#f79cd4",
    }
    dataset = {
        "Adams_Kaminski_2020": "#023fa5",
        "Chen_Zhang_2020": "#7d87b9",
        "Goveia_Carmeliet_2020": "#bec1d4",
        "Guo_Zhang_2018": "#d6bcc0",
        "Habermann_Kropski_2020": "#bb7784",
        "He_Fan_2021": "#8e063b",
        "Kim_Lee_2020": "#4a6fe3",
        "Lambrechts_Thienpont_2018_6149v1": "#8595e1",
        "Lambrechts_Thienpont_2018_6149v2": "#b5bbe3",
        "Lambrechts_Thienpont_2018_6653": "#e6afb9",
        "Laughney_Massague_2020": "#e07b91",
        "Madissoon_Meyer_2020": "#d33f6a",
        "Maier_Merad_2020": "#11c638",
        "Maynard_Bivona_2020": "#8dd593",
        "Mayr_Schiller_2020": "#c6dec7",
        "Reyfman_Misharin_2018": "#ead3c6",
        "Travaglini_Krasnow_2020": "#f0b98d",
        "UKIM-V": "#ef9708",
        "Vieira_Teichmann_2019": "#0fcfc0",
        "Wu_Zhou_2021": "#9cded6",
        "Zilionis_Klein_2019": "#d5eae7",
    }

    scissor = {}
