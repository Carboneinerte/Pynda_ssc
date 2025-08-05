### Cluster cells using a partition algorithm
from banksy.cluster_methods import run_Leiden_partition
import pandas as pd
import pickle
from datetime import datetime
import numpy as np
import os
import scanpy as sc
from banksy_utils.cluster_utils import pad_clusters, create_spatial_nonspatial_adata
import scipy.sparse as sparse
from banksy_utils.load_data import load_adata, display_adata ### is it useful?
from banksy_utils.plot_utils import plot_qc_hist, plot_cell_positions
from banksy_utils.filter_utils import normalize_total, filter_hvg, print_max_min
from banksy.main import median_dist_to_nearest_neighbour
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy.main import concatenate_all

def print_with_elapsed_time(message):
    elapsed_time = datetime.now() - start_time
    elapsed_seconds = elapsed_time.total_seconds()
    print(f'[{elapsed_seconds:.2f} seconds] {message}')


start_time = datetime.now()
resolutions = [0.7]
dir_notebook = '/media/volume/volume_spatial/hugo/notebook'
name_dir = "MB_test"
seed = 1234


with open(f'{dir_notebook}/dict/banksy_dict_{name_dir}.pkl', 'rb') as f:
    banksy_dict = pickle.load(f)

# print(banksy_dict)

### Reduce dimensions of each data matrix
from banksy_utils.umap_pca import pca_umap
resolutions = [0.7] # clustering resolution for UMAP
pca_dims = [20] # Dimensionality in which PCA reduces to
lambda_list = [0.2] # list of lambda parameters

start_time = datetime.now()
print_with_elapsed_time(f"Start DimRed")
pca_umap(banksy_dict,
         pca_dims = pca_dims,
         add_umap = True,
         plt_remaining_var = False,
         )
print_with_elapsed_time(f"End DimRed")

if not os.path.exists(f"dict"):
   os.makedirs(f"dict")

with open(f'dict/banksy_dict_{name_dir}_dimred.pkl', 'wb') as f: 
    pickle.dump(banksy_dict, f)

print_with_elapsed_time(f'start Clustering')
                                        
results_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn = 50,
    num_iterations = -1,
    partition_seed = seed,
    match_labels = False,
)
 
print_with_elapsed_time('Finish Clustering')
np.save(f'{dir_notebook}/dict/max_num_labels', max_num_labels)

if not os.path.exists(f'{dir_notebook}/dict'):
   os.makedirs(f'{dir_notebook}/dict')

with open(f'{dir_notebook}/dict/banksy_dict_all-samples_results_df.pkl', 'wb') as f: 
    pickle.dump(results_df, f)

print_with_elapsed_time('Add needed columns')
coord_keys = ('xcoord', 'ycoord', 'coord_xy')
file_path = os.path.join("Banksy_py", "data", "slide_seq", "v1")
c_map =  'tab20b' # specify color map
weights_graph =  banksy_dict['scaled_gaussian']['weights'][0]

for params_name in results_df.index:

    label_index = 'labels'

    labels = results_df.loc[params_name, label_index]
    adata_temp = results_df.loc[params_name, "adata"]
    num_pcs = results_df.loc[params_name, 'num_pcs']

    pc_temp = adata_temp.obsm[f"reduced_pc_{num_pcs}"]
    umap_temp = adata_temp.obsm[f"reduced_pc_{num_pcs}_umap"]

    label_name = f"labels_{params_name}"
    adata_temp.obs[label_name] = np.char.mod('%d', labels.dense)
    adata_temp.obs[label_name] = adata_temp.obs[label_name].astype('category')

    adata_temp.obsm[coord_keys[2]] = np.vstack(
        (adata_temp.obs[coord_keys[0]].values,
            adata_temp.obs[coord_keys[1]].values)
    ).T

print_with_elapsed_time('Create adata_spatial')
# Here we manually assign clusters to their identity using a dictionary
cluster2annotation_spatial = {}
pad_clusters(cluster2annotation_spatial, list(range(max_num_labels)))
cluster2annotation_nonspatial = {}
print(cluster2annotation_spatial,"\n", cluster2annotation_nonspatial)

lambda_list = [0.2] # list of lambda parameters
pca_dims = [20] # Dimensionality in which PCA reduces to
resolutions = [0.7]

# save annotations in two different anndata objects (adata_spatial and adata_nonspatial)
adata_spatial, adata_nonspatial = create_spatial_nonspatial_adata(results_df,
                                    pca_dims,
                                    lambda_list, 
                                    resolutions,
                                    cluster2annotation_spatial,
                                    cluster2annotation_nonspatial
                                                                 )

### Automatic initial annotation
print_with_elapsed_time('Start Automatic annotation')
cont_tab = pd.crosstab(adata_spatial.obs['labels_scaled_gaussian_pc20_nc0.20_r0.70'], adata_spatial.obs['mmc:subclass_name'], normalize="index")
cont_tab_trans = cont_tab.T
max_col_dict = cont_tab_trans.idxmax(axis=0).to_dict()
adata_spatial.obs['cell_type_auto'] = adata_spatial.obs['labels_scaled_gaussian_pc20_nc0.20_r0.70'].map(max_col_dict)

all_cell_type = adata_spatial.obs['cell_type_auto'].unique()
list_cell_nb = range(0, len(all_cell_type))
mapping_dict = dict(zip(all_cell_type,list_cell_nb))
adata_spatial.obs['cell_type_newnum'] = adata_spatial.obs['cell_type_auto'].map(mapping_dict)
mapping_dict

print_with_elapsed_time('Start saving')

adata_spatial.write(f"{dir_notebook}/h5ad/{name_dir}/{name_dir}_Banksy_combined.h5ad.gz", compression='gzip')

print_with_elapsed_time(f"End")