#!/usr/bin/env python

"""
Requirements:
    * Python >= 3.6.2
    * pybedtools
    * pandas
    * numpy

Copyright (c) 2022 Dietmar Rieder <dietmar.rieder@i-med.ac.at>
MIT License <http://opensource.org/licenses/MIT>


"""

RELEASE = False
__version_info__ = (
    "0",
    "1",
)
__version__ = ".".join(__version_info__)
__version__ += "-dev" if not RELEASE else ""


import pandas as pd
import pybedtools as bt
import argparse
import sys
import os
import numpy as np


def is_valid(x):
    try:
        valid = (len(x) == 2) and (x[0] <= x[1])
    except:
        valid = False
    finally:
        return valid


def is_overlapping(x, y):
    return max(x[0], y[0]) <= min(x[1], y[1])


def overlapping_fraction(x, y):
    if (is_valid(x) and is_valid(y)) == False:
        raise ValueError("Invalid range")

    if is_overlapping(x, y):
        overlapping_distance = min(x[1], y[1]) - max(x[0], y[0]) + 1
        width_x = x[1] - x[0] + 1
        width_y = y[1] - y[0] + 1
        overlap_x = overlapping_distance / width_y
        overlap_y = overlapping_distance / width_x
        # return (overlap_x, overlap_y)
        return overlap_y
    return np.nan


def weighted_average(data):
    d = {}

    d["seg_mean_wa"] = np.average(data["Segment_Mean"], weights=data["weight"])

    return pd.Series(d)


if __name__ == "__main__":

    usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(
        description="Make intersected BED file of gene level CNV and segmented CNV files from GDC"
    )

    parser.add_argument(
        "--gene_cnv",
        required=True,
        type=str,
        help="Gene level CNV file",
    )
    parser.add_argument(
        "--seg_cnv",
        required=True,
        type=str,
        help="Segmented CNV file",
    )
    parser.add_argument(
        "--outfile",
        required=True,
        type=str,
        help="BED output",
    )

    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)

    args = parser.parse_args()

    gene_cnv_f = args.gene_cnv
    seg_cnv_f = args.seg_cnv
    bed_out_f = args.outfile

    gene_cnv = pd.read_table(gene_cnv_f)
    seg_cnv = pd.read_table(seg_cnv_f)

    gene_cnv = gene_cnv[["chromosome", "start", "end", "gene_name", "gene_id"]]
    gene_cnv_bed = bt.BedTool.from_dataframe(gene_cnv)

    seg_cnv.reset_index(inplace=True)
    seg_cnv = seg_cnv.rename(columns={"index": "Seg"})
    seg_cnv["Seg"] = "seq_" + seg_cnv["Seg"].astype(str)
    seg_cnv["Chromosome"] = "chr" + seg_cnv["Chromosome"].astype(str)

    seg_cnv = seg_cnv[["Chromosome", "Start", "End", "Seg", "Segment_Mean", "Num_Probes"]]
    seg_cnv_bed = bt.BedTool.from_dataframe(seg_cnv)

    gene_seg_cnv = gene_cnv_bed.intersect(seg_cnv_bed, loj=True).to_dataframe(
        disable_auto_names=True,
        names=[
            "chromosome",
            "start",
            "end",
            "gene_name",
            "gene_id",
            "Chromosome",
            "Start",
            "End",
            "Seg",
            "Segment_Mean",
            "Num_Probes",
        ],
    )
    gene_seg_cnv = gene_seg_cnv.replace(".", np.nan).astype({"Segment_Mean": "float"})

    gene_seg_cnv["weight"] = gene_seg_cnv.apply(lambda x: overlapping_fraction((x.start, x.end), (x.Start, x.End)), 1)

    gene_seg_cnv_agg = gene_seg_cnv.groupby(["chromosome", "start", "end", "gene_name", "gene_id"]).apply(
        weighted_average
    )

    gene_seg_cnv_agg.to_csv(bed_out_f, sep="\t")
