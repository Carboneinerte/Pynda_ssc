import scanpy as sc
from xenium_preprocessing import import_xenium
from config_local import dir_raw, dir_processed, sample_name_import
import pandas as pd
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description='Convert cell boundaries to GeoJSON.'
    )

parser.add_argument('--name_dir', type=str, help='Name of the experience')
args = parser.parse_args()

name_dir = args.name_dir
samples_ids = sample_name_import(name_dir)

def create_folders(name_dir:str,
                   dir_processed: str = dir_processed
                   ):
    if not os.path.exists(f"{dir_processed}/csv/{name_dir}/"):
        os.makedirs(f"{dir_processed}/csv/{name_dir}/")
        print("csv folder created")
    if not os.path.exists(f"{dir_processed}/h5ad/{name_dir}/"):
        os.makedirs(f"{dir_processed}/h5ad/{name_dir}/")
        print("h5ad folder created")
    if not os.path.exists(f"{dir_processed}/analysis/{name_dir}/"):
        os.makedirs(f"{dir_processed}/analysis/{name_dir}/")
        print('Analysis folder created')
    if not os.path.exists(f"{dir_processed}/plot/{name_dir}/"):
        os.makedirs(f"{dir_processed}/plot/{name_dir}/")
        print('Plotfolder created')

def undernoise_list(samples_ids:list,
                    name_dir:str,
                    dir_raw:str = dir_raw,
                    dir_processed:str = dir_processed,
                    ):

    for idx, sample in enumerate(samples_ids):
        print('Start Sample :', sample)
        print(idx+1," / ", len(samples_ids))
        df = pd.read_csv(f'{dir_raw}/{name_dir}/{sample}/transcripts.csv.gz')
        df = df[df['qv']>=20]
        
        data = pd.DataFrame({'feature_name': df.feature_name.value_counts().index,'count' : df.feature_name.value_counts()})
        data.sort_index(inplace=True)

        data['type'] = data['feature_name'].apply(lambda x: x.split('_')[0])

        percentile_threshold:float = 99.5
        threshold = np.percentile(data[data['type']=="NegControlProbe"]['count'].values,percentile_threshold)
        print('threshold = ', threshold)

        data['logfoldovernoise'] = data['count'].apply(lambda x: np.log(x / threshold))
        data_gen_only = data[~(data['feature_name'].str.contains('_'))]
        print('nb of gene under threshold : ', data_gen_only[data_gen_only['logfoldovernoise']<0].shape[0])
        if idx == 0:
            set_undernoise = set(data_gen_only[data_gen_only['logfoldovernoise']<0]['feature_name'].values)
        else:
            set_undernoise = set_undernoise.intersection(set(data_gen_only[data_gen_only['logfoldovernoise']<0]['feature_name'].values))

        print(" ")
        
    pd.Series(list(set_undernoise)).to_csv(f'{dir_processed}/analysis/{name_dir}/undernoise_{name_dir}.csv')
    list_noise = list(set_undernoise)
    return list_noise

def import_xenium(samples_ids:  list,
                  name_dir:     str,
                  dir_raw:      str = dir_raw,
                  dir_processed: str = dir_processed,
                  trans_min:    int = 40,
                  trans_max:    int = 4000,
                  remove_noise: bool = False,
                  MMC:          bool = False,
                  ):
    '''
    dir_raw (str) : folder containing raw Xenium files
    dir_processed (str)
    samples_ids (str)
    name_dir (str)
    remove_noise (bool) : remove genes below noise level from a list
    MMC (bool)
    '''
    create_folders(name_dir)

    adatas = []
    if remove_noise == True:
        print("## Noise evaluation ##")
        list_noise = undernoise_list(samples_ids, name_dir)
        print(f"Will exclude {len(list_noise)} genes")
    print(" ")
    print("## Start importation ##")

    for sample_id in samples_ids:
        adata = sc.read_10x_h5(f"{dir_raw}/{name_dir}/{sample_id}/cell_feature_matrix.h5")
        df = pd.read_csv(f"{dir_raw}/{name_dir}/{sample_id}/cells.csv.gz")
        df.set_index(adata.obs_names, inplace=True)
        adata.obs = df.copy()
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
        adata.layers["counts"] = adata.X.copy()
        all_cells = adata.shape[0]

        if remove_noise:
            mask = [gene not in list_noise for gene in adata.var_names]
            adata = adata[:, mask].copy()
  
        
        sc.pp.filter_cells(adata, max_counts=trans_max) ## Possible filter to remove cells with too many transcripts
        sc.pp.filter_cells(adata, min_counts=trans_min) ## Filter cells with less than 40 transcripts
        sc.pp.filter_genes(adata, min_cells=5) ## Filter genes expressed in less than 5 cells
        adata.obs_names = [f"{sample_id}_{cell_id}" for cell_id in adata.obs_names]
        adata.obs['cell_id'] = adata.obs_names
        print(f"Proportion of cells concerved after filtering = {adata.shape[0] / all_cells:.2%} ({adata.shape[0]} cells)")
        
        adatas.append(adata)
        print(f"Sample {sample_id} done")
        print(" ")
        if MMC:
            if not os.path.exists(f"{dir_processed}/Correlation_Mapping/{name_dir}/"):
                os.makedirs(f"{dir_processed}/Correlation_Mapping/{name_dir}/")
                print("Correlation_Mapping folder created")
            adata.write(f"{dir_processed}/h5ad/{name_dir}/{name_dir}_{sample_id}_forMMC.h5ad")

    print(f"Read all {len(samples_ids)} samples")

    ### merge all the anndata objects into a single object
    adata = adatas[0].concatenate(adatas[1:], index_unique=None)

    ### Add a sample column to the metadata
    adata.obs['sample'] = adata.obs_names.map(lambda name: name.split('_')[0])
    # samples = adata.obs['sample'].unique()
    adata.write(f"{dir_processed}/h5ad/{name_dir}/{name_dir}_import.h5ad.gz", compression = "gzip")
    return adata


import_xenium(samples_ids=samples_ids,
                  name_dir=name_dir,
                  dir_raw = dir_raw,
                  dir_processed = dir_processed,
                  trans_min = 40,
                  trans_max = 4000,
                  remove_noise = True,
                  MMC = True,
                  )