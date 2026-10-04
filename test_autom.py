

import scanpy as sc
import anndata

anndata.settings.allow_write_nullable_strings = True

print("load file")

name_dir = "../data_processed/ret_data/"
path = "../data_raw/data/ret_annot_subsample.h5ad.gz"
adata = sc.read_h5ad(path)

print("Start")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
print(f"Normalize done")
sc.pp.log1p(adata)
print(f"log1p done")
sc.pp.pca(adata, n_comps = 10)
print("pca done")
sc.pp.neighbors(adata)
print('neighbors done')
sc.tl.umap(adata, min_dist=1)
print("umap done")
sc.tl.leiden(adata, resolution = 0.7, flavor = "igraph", n_iterations = 2, key_added = "test_autom")
print("clustering done, start saving")

adata.write(f"{name_dir}/test_automate_snake.h5ad.gz", compression='gzip')

print("done")
