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

    args = parser.parse_args()

    sample = args.InFile
    s_name = sample.replace(".h5ad", "")
    outfile = args.outDir + os.path.basename(s_name) + ".pkl"
    adata_sample = anndata.read_h5ad(sample)

    res = sq.gr.ligrec(
        adata_sample,
        n_perms=1000,
        cluster_key="cell_type",
        copy=True,
        use_raw=False,
        transmitter_params={"categories": "ligand"},
        receiver_params={"categories": "receptor"},
        interactions_params={"resources": "CellPhoneDB"},
        n_jobs=10,
    )

    with open(outfile, "wb") as handle:
        pickle.dump(res, handle, protocol=pickle.HIGHEST_PROTOCOL)
