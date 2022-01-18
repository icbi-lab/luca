#!/usr/bin/env python

import anndata, argparse, os, pickle
import pandas as pd
import squidpy as sq
from tqdm import tqdm


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        usage="squidpy_cpdb.py [-h] -i {/path/to/counts_tpm/input/filename.h5ad}"
    )
    parser.add_argument("-i", "--InFile", help="H5AD sample file", required=True)
    parser.add_argument("-o", "--outDir", help="pickle output dir path", required=True)
    parser.add_argument(
        "-n", "--cpus", help="number of cpus for ligrec", required=True, type=int
    )
    parser.add_argument(
        "-c",
        "--cellTypeKey",
        help="column in adata.obs holding cell-type information",
        required=True,
    )

    args = parser.parse_args()

    sample = args.InFile
    ct_key = args.cellTypeKey
    s_name = sample.replace(".h5ad", "")
    outfile = args.outDir + os.path.basename(s_name) + ".pkl"
    adata_sample = anndata.read_h5ad(sample)
    cpus = args.cpus

    # squidpy expects log-norm data, which should be in adata.raw, see
    # https://github.com/theislab/squidpy/issues/446#issuecomment-1013244936
    res = sq.gr.ligrec(
        adata_sample,
        n_perms=0,
        cluster_key=ct_key,
        copy=True,
        use_raw=True,
        threshold=0,
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        interactions_params={"resources": "CellPhoneDB"},
        n_jobs=cpus,
    )

    with open(outfile, "wb") as handle:
        pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)
