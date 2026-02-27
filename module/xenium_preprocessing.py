import pandas as pd
import scanpy as sc
import numpy as np
import os
from datetime import datetime
from matplotlib.pyplot import rc_context
from module.misc import list_annotations
from module.dataviz_analysis import cluster_plot

def create_folders(dir_notebook, name_dir):
    if not os.path.exists(f"{dir_notebook}/csv/{name_dir}/"):
        os.makedirs(f"{dir_notebook}/csv/{name_dir}/")
        print("csv folder created")
    if not os.path.exists(f"{dir_notebook}/h5ad/{name_dir}/"):
        os.makedirs(f"{dir_notebook}/h5ad/{name_dir}/")
        print("h5ad folder created")
    if not os.path.exists(f"{dir_notebook}/analysis/{name_dir}/"):
        os.makedirs(f"{dir_notebook}/analysis/{name_dir}/")
        print('Analysis folder created')
    if not os.path.exists(f"{dir_notebook}/plot/{name_dir}/"):
        os.makedirs(f"{dir_notebook}/plot/{name_dir}/")
        print('Plotfolder created')

def undernoise_list(dir:str, dir_notebook:str, samples_ids:list, name_dir:str):

    for idx, sample in enumerate(samples_ids):
        print('Start Sample :', sample)
        print(idx+1," / ", len(samples_ids))
        if sample.split('-')[0]=="SD1":
            df = pd.read_parquet(f'E:/Xenium_SD/{sample}/transcripts.parquet',
                            filters=[("qv",">=",20)]
                            )
        else:
            df = pd.read_parquet(f'{dir}/{sample}/transcripts.parquet',
                            filters=[("qv",">=",20)]
                            )

        
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
        
    pd.Series(list(set_undernoise)).to_csv(f'{dir_notebook}/analysis/{name_dir}/undernoise_{name_dir}.csv')
    list_noise = list(set_undernoise)
    return list_noise

def import_xenium(dir, dir_notebook, samples_ids, name_dir, trans_min: int = 40, trans_max: int = 4000, remove_noise=False, MMC = False):
    '''
    dir (str) : folder containing raw Xenium files
    dir_notebook (str)
    samples_ids (str)
    name_dir (str)
    remove_noise (bool) : remove genes below noise level from a list
    MMC (bool)
    '''
    create_folders(dir_notebook, name_dir)

    adatas = []
    
    if remove_noise == True:
        print("## Noise evaluation ##")
        list_noise = undernoise_list(dir, dir_notebook, samples_ids, name_dir)
        print(f"Will exclude {len(list_noise)} genes")
    print(" ")
    print("## Start importation ##")
    for sample_id in samples_ids:
        adata = sc.read_10x_h5(f"{dir}/{sample_id}/cell_feature_matrix.h5")
        df = pd.read_csv(f"{dir}/{sample_id}/cells.csv.gz")
        df.set_index(adata.obs_names, inplace=True)
        adata.obs = df.copy()
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
        adata.layers["counts"] = adata.X.copy()

        if remove_noise == True:
            mask = [gene not in list_noise for gene in adata.var_names]
            adata = adata[:, mask].copy()
            
        all_cells = adata.shape[0]
        sc.pp.filter_cells(adata, max_counts=trans_max) ## Possible filter to remove cells with too many transcripts
        sc.pp.filter_cells(adata, min_counts=trans_min) ## Filter cells with less than 40 transcripts
        sc.pp.filter_genes(adata, min_cells=5) ## Filter genes expressed in less than 5 cells
        adata.obs_names = [f"{sample_id}_{cell_id}" for cell_id in adata.obs_names]
        adata.obs['cell_id'] = adata.obs_names
        print(f"Proportion of cells concerved after filtering = {adata.shape[0] / all_cells:.2%} ({adata.shape[0]} cells)")
        
        adatas.append(adata)
        # adata.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_{sample_id}_forMMC.h5ad")
        print(f"Sample {sample_id} done")
        print(" ")
        if MMC == True:
            if not os.path.exists(f"{dir_notebook}/h5ad/{name_dir}/"):
                os.makedirs(f"{dir_notebook}/h5ad/{name_dir}/")
            adata.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_{sample_id}_forMMC.h5ad")

    print(f"Read all {len(samples_ids)} samples")

    ### merge all the anndata objects into a single object
    adata = adatas[0].concatenate(adatas[1:], index_unique=None)

    ### Add a sample column to the metadata
    adata.obs['sample'] = adata.obs_names.map(lambda name: name.split('_')[0])
    # samples = adata.obs['sample'].unique()
    adata.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_import.h5ad.gz", compression = "gzip")
    return adata

def mmc_merge(adata, dir_notebook: str, name_dir: str):
    '''
    Merge MMC correlation mapping from {dir_notebook}/Correlation_Mapping/ folder
    MMC files names should start with {name_dir} 
    adata : AnnData object
    dir_notebook : string
    name_dir : string
    '''

    import glob
    dir_corr = f'{dir_notebook}/Correlation_Mapping/'

    list_files = glob.glob(dir_corr + f'{name_dir}*')

    if len(list_files) != 0:

        for files in list_files:
            print(files)
            if 'corr' not in locals():
                corr = pd.read_csv(files, comment = '#')
            else:
                csv_temp = pd.read_csv(files, comment = '#')
                corr = pd.concat([corr, csv_temp], ignore_index=True)

        HC3_MMC = corr
        HC3_MMC.index = HC3_MMC['cell_id']
        HC3_MMC.index.name = None
        HC3_MMC.columns = [f"mmc:{i}" for i in HC3_MMC.columns]
        mmc_dict_class = dict(zip(HC3_MMC['mmc:cell_id'], HC3_MMC['mmc:class_name']))
        mmc_dict_classcoef = dict(zip(HC3_MMC['mmc:cell_id'], HC3_MMC['mmc:class_correlation_coefficient']))
        mmc_dict_subclass = dict(zip(HC3_MMC['mmc:cell_id'], HC3_MMC['mmc:subclass_name']))
        mmc_dict_supertype = dict(zip(HC3_MMC['mmc:cell_id'], HC3_MMC['mmc:supertype_name']))

        adata.obs['mmc:class_name'] = adata.obs['cell_id'].map(mmc_dict_class)
        adata.obs['mmc:class_correlation_coefficient'] = adata.obs['cell_id'].map(mmc_dict_classcoef)
        adata.obs['mmc:subclass_name'] = adata.obs['cell_id'].map(mmc_dict_subclass)
        adata.obs['mmc:supertype_name'] = adata.obs['cell_id'].map(mmc_dict_supertype)
    else:
        print('Empty input. Please check the names of files in Correlation_Mapping folder')
    return adata

def add_annotations(adata, df):
    '''
    Add annotations from adata.obs to matrix of gene expression (usually called 'df').
    List of annotations can be changed in module/misc.py
    '''

    if 'cell_id' not in df.columns:
        df['cell_id'] = df.index
    
    list_anno = list_annotations()
    list_anno = [anno for anno in list_anno if anno in adata.obs.columns]
    for anno in list_anno:
        df[anno] = df['cell_id'].map(dict(zip(adata.obs['cell_id'], adata.obs[anno])))
    
    return df

def add_annotations_unassigned(adata, df):
    '''
    Add annotations from adata.obs to matrix of gene expression (usually called 'df').
    List of annotations can be changed in module/misc.py
    '''

    if 'gridcell_id' not in df.columns:
        df['gridcell_id'] = df.index
    
    list_anno = list_annotations()
    list_anno = [anno for anno in list_anno if anno in adata.obs.columns]
    for anno in list_anno:
        df[anno] = df['gridcell_id'].map(dict(zip(adata.obs['gridcell_id'], adata.obs[anno])))
    
    return df


def normalization_scanpy(adata):
    print(f"Start")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    print(f"Normalize done")
    sc.pp.log1p(adata)
    print(f"log1p done")

    return adata

def clustering_scanpy(adata, pca_compo: int = 10, leiden_resolution: float = 0.7, ):
    print(f"Start Clustering")
    sc.pp.pca(adata, n_comps = pca_compo)
    print(f"PCA done")
    sc.pl.pca_variance_ratio(adata, n_pcs=pca_compo, log=False)
    sc.pl.pca(
        adata,
        color=["sample","sample"],
        dimensions=[(0, 1), (2, 3)],
        ncols=2,
        size=1,
    )
    sc.pp.neighbors(adata)
    print(f"Neighbors done")
    sc.tl.umap(adata, min_dist = 1)
    print(f"UMAP done")

    sc.tl.leiden(adata, resolution = leiden_resolution) ### Use a higher resolution value to obtain more clusters. They can be adjusted/merged/subclustered later
    print('End of clustering')
    with rc_context({"figure.figsize": (10, 10)}):
        sc.pl.umap(
            adata,
            color="leiden",
            add_outline=False,
            legend_loc="on data",
            legend_fontsize=12,
            legend_fontoutline=2,
            frameon=False,
            palette="tab20b",
            size = 1
        )

    

    return adata

def norm_to_cluster(adata, dir_notebook, name_dir,pca_compo:int = 10,leiden_resolution:float = 0.7):
    
    adata = normalization_scanpy(adata)
    
    adata = clustering_scanpy(adata,pca_compo,leiden_resolution)

    cluster_plot(adata, dir_notebook, name_dir, cluster_to_use = "leiden",cluster_to_map=['all'])

    return adata