import scanpy as sc
from module.config_local import dir_raw, dir_processed
import pandas as pd
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Import Xenium raw files and exctract gene expression for MMC.'
    )

parser.add_argument('--name_dir', type=str, help='Name of the experience')
parser.add_argument('--sample', type=str, required=True, help='Sample name to process')

args = parser.parse_args()
name_dir = args.name_dir
sample_id = args.sample

def import_MMC(sample_id:  str,
                  name_dir:     str,
                  dir_raw:      str = dir_raw,
                  dir_processed: str = dir_processed,
                  ):
    '''
    dir_raw (str) : folder containing raw Xenium files
    dir_processed (str)
    samples_ids (str)
    name_dir (str)
    remove_noise (bool) : remove genes below noise level from a list
    MMC (bool)
    '''

    print(" ")
    print("## Start importation ##")

    adata = sc.read_10x_h5(f"{dir_raw}/{name_dir}/{sample_id}/cell_feature_matrix.h5")
    df = pd.read_csv(f"{dir_raw}/{name_dir}/{sample_id}/cells.csv.gz")
    df.set_index(adata.obs_names, inplace=True)
    adata.obs = df.copy()
    adata.layers["counts"] = adata.X.copy()
    all_cells = adata.shape[0]
    sc.pp.filter_genes(adata, min_cells=5) ## Filter genes expressed in less than 5 cells
    adata.obs_names = [f"{sample_id}_{cell_id}" for cell_id in adata.obs_names]
    adata.obs['cell_id'] = adata.obs_names
    
    print(f"Sample {sample_id} done")
    print(" ")
    adata.write(f"{dir_processed}/{name_dir}/h5ad/{sample_id}_forMMC.h5ad")


import_MMC(sample_id=sample_id,
            name_dir=name_dir,
            dir_raw = dir_raw,
            dir_processed = dir_processed,
            )
