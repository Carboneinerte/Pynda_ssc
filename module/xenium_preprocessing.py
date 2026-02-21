import pandas as pd
import scanpy as sc
import os
from module.misc import list_annotations

import pandas as pd
import scanpy as sc
import os
from module.misc import list_annotations

def import_xenium(dir, dir_notebook, samples, samples_ids, name_dir, remove_noise=False):
    adatas = []
    for sample, sample_id in zip(samples, samples_ids):
        adata = sc.read_10x_h5(f"{dir}/{sample}/cell_feature_matrix.h5")
        df = pd.read_csv(f"{dir}/{sample}/cells.csv.gz")
        df.set_index(adata.obs_names, inplace=True)
        adata.obs = df.copy()
        adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
        adata.layers["counts"] = adata.X.copy()

        if remove_noise == True:
            path_to_list = f'undernoise_{name_dir}.csv'
            list_noise = pd.read_csv(path_to_list)
            list_noise = list(list_noise['0'].values)
            list_noise[0:5]
            mask = [gene not in list_noise for gene in adata.var_names]
            adata = adata[:, mask].copy()
            
        all_cells = adata.shape[0]
        # sc.pp.calculate_qc_metrics(adata,  percent_top=(10, 20, 50, 150), inplace=True)
        sc.pp.filter_cells(adata, max_counts=4000) ## Possible filter to remove cells with too many transcripts
        sc.pp.filter_cells(adata, min_counts=40) ## Filter cells with less than 40 transcripts
        sc.pp.filter_genes(adata, min_cells=5) ## Filter genes expressed in less than 5 cells
        adata.obs_names = [f"{sample_id}_{cell_id}" for cell_id in adata.obs_names]
        adata.obs['cell_id'] = adata.obs_names
        print(f"Proportion of cells concerved after filtering = {adata.shape[0] / all_cells:.2%}")
        
        adatas.append(adata)
        # adata.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_{sample_id}_forMMC.h5ad")
        print(f"Sample {sample} done")

        if not os.path.exists(f"{dir_notebook}/h5ad/{name_dir}/"):
            os.makedirs(f"{dir_notebook}/h5ad/{name_dir}/")
        adata.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_{sample}_forMMC.h5ad")

    print(f"Read all {len(samples)} samples")

    ### merge all the anndata objects into a single object
    adata = adatas[0].concatenate(adatas[1:], index_unique=None)

    ### Add a sample column to the metadata
    adata.obs['sample'] = adata.obs_names.map(lambda name: name.split('_')[0])
    # samples = adata.obs['sample'].unique()

    return adata

def mmc_merge(adata, dir_notebook, name_dir):

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
    if 'cell_id' not in df.columns:
        df['cell_id'] = df.index
    
    list_anno = list_annotations()
    list_anno = [anno for anno in list_anno if anno in adata.obs.columns]
    for anno in list_anno:
        df[anno] = df['cell_id'].map(dict(zip(adata.obs['cell_id'], adata.obs[anno])))
    
    return df