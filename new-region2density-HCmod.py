import polars as pl
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import time
import anndata as ad
import scanpy as sc
import os

save_path = "/media/volume/volume_spatial/hugo/notebook/h5ad/AD_unassigned/"

ppath = "/media/volume/volume_spatial/hugo/data/"
#pfilelist = ['CTX_transcripts_circa4-IGM-ZT01.parquet','SCH_transcripts_circa4-IGM-ZT01.parquet']
pfilelist = ["3161-1","3159-2"
    # 'circa4-IGM-ZT13',
#     "SD1-ZT01","SD1-ZT05","SD1-ZT09","SD1-ZT13","SD1-ZT17","SD1-ZT21",
#  'circa4-IGM-ZT01','circa4-IGM-ZT05','circa4-IGM-ZT09','circa4-IGM-ZT13','circa4-IGM-ZT17','circa4-IGM-ZT21'
 ]
#out_path = "density_bins_by_region-test.csv"
#outpathcounter = 0 #lazy way to make unique outfiles. you can change this to a better naming convention
BIN = 50   # your x (same for width/height)
x0 = 0.0          # grid origin; set to 0 or use min below
y0 = 0.0

adatas = []

coordpath = '/media/volume/volume_spatial/hugo/notebook/coordinates/whole_brain/'

def bins_for_geom(geom, bin_size, x0=0.0, y0=0.0):
        minx, miny, maxx, maxy = geom.bounds

        gx0 = int(np.floor((minx - x0) / bin_size) * bin_size + x0)
        gx1 = int(np.floor((maxx - x0) / bin_size) * bin_size + x0)
        gy0 = int(np.floor((miny - y0) / bin_size) * bin_size + y0)
        gy1 = int(np.floor((maxy - y0) / bin_size) * bin_size + y0)
        #print(x0,gx0,minx,y0,gy0,miny) i think the floor function means that these will always align wih the parquet grid, but may need to do some testing
        xs = np.arange(gx0, gx1 + bin_size, bin_size, dtype=np.int64)
        ys = np.arange(gy0, gy1 + bin_size, bin_size, dtype=np.int64)

        keys, bx, by, cells = [], [], [], []
        for x in xs:
            for y in ys:
                bx.append(x)
                by.append(y)
                keys.append(f"{x}_{y}")
                cells.append(box(x, y, x + bin_size, y + bin_size))

        return keys, bx, by, cells


def clipped_bin_areas(
    gdf,
    name_col="names",
    bin_size=BIN,
):
    rows = []

    for _, r in gdf.iterrows():
        name = r[name_col]
        geom = r.geometry

        keys, bx, by, cells = bins_for_geom(geom, bin_size=bin_size)

        cells_gdf = gpd.GeoDataFrame(
            {
                "bin_key": keys,   # <-- FIRST COLUMN
                "bin_x": bx,
                "bin_y": by,
            },
            geometry=cells,
            crs=gdf.crs,
        )

        inter = cells_gdf.geometry.intersection(geom)
        area = inter.area

        mask = area > 0
        if mask.any():
            out = cells_gdf.loc[
                mask, ["bin_key", "bin_x", "bin_y"]
            ].copy()

            out.insert(0, "region", name)   # bin_key stays first
            out["area"] = area[mask].to_numpy()

            rows.append(out)

    if not rows:
        return pd.DataFrame(
            columns=["region", "bin_key", "bin_x", "bin_y", "area"]
        )
    
    return pd.concat(rows, ignore_index=True)

for pfile in pfilelist:
    sti = time.time()
    parquet_in = ppath+pfile+'/transcripts.parquet'#"/media/volume/volume_spatial/hugo/data/circa4-IGM-ZT01/splitregion/CTX_transcripts_circa4-IGM-ZT01.parquet"
    #parquet_out = "coords_binned-test.parquet"
    out_path = "data/density-bins-5k/"+pfile+"-densitybins.csv"
    print(pfile,'bin size = ',BIN)
    #this takes the parquet file and saves a new parquet file with bin columns. you could rename the out file if you want to save it in place
    df = (
        pl.scan_parquet(parquet_in)  # lazy scan (doesn't load all at once)
        .filter(~pl.col("feature_name").str.contains("_") & (pl.col("qv") > 20)) #should filter out the codewords and low quality counts
        .with_columns([
            (((pl.col("x_location") - x0) / BIN).floor()*BIN).cast(pl.Int64).alias("bin_x"),
            (((pl.col("y_location") - y0) / BIN).floor()*BIN).cast(pl.Int64).alias("bin_y"),
        ])
        # optional single key:
        .with_columns([
            (pl.col("bin_x").cast(pl.Utf8) + pl.lit("_") + pl.col("bin_y").cast(pl.Utf8)).alias("bin_key")
        ])
    )
    print('binned transcripts')
    #df.sink_parquet(parquet_out)  # writes efficiently



    #load in shapes. will have to loop through these too
    #dfg_r = gpd.read_file("/media/volume/volume_spatial/hugo/notebook/coordinates/Region_prediction/Xenium-data-coordinates-filtered_circa4-IGM-ZT01.geojson")
    dfg_r = gpd.read_file(coordpath+pfile+"_wholebrain_annotation.geojson")
    # this gets rid of the warnings Only do this if the CRS is currently geographic

    if dfg_r.crs is None or dfg_r.crs.is_geographic:
        dfg_r = dfg_r.set_crs("EPSG:3857", allow_override=True)
    #pctx = dfg_r.loc[dfg_r['cell_type_newnum_final']=="CTX",'geometry'].values[0]

    #get the density df
    print("checking region boundary/bin overlaps and getting area")

    # usage:
    # result = clipped_bin_areas(my_gdf, name_col="names", bin_size=50)
    # result columns: name, bin_x, bin_y, area
    cba = clipped_bin_areas(dfg_r, name_col="region", bin_size=50)
    #warings are assuming geographic corrdinates that need to take into account curvature. ok to ignore i think


    #parquet_path = "coords_binned-test.parquet"  # the parquet that already has bin_key (or bin_x/bin_y)
    print("adding U/A counts/density")
    bin_counts = (
        df#pl.scan_parquet(parquet_path)
        .with_columns([ #modify this to remove codewords
            (pl.col("cell_id") == "UNASSIGNED").cast(pl.Int64).alias("is_unassigned"),
            (pl.col("cell_id") != "UNASSIGNED").cast(pl.Int64).alias("is_assigned"),
        ])
        .group_by(["bin_key","feature_name"])
        .agg([
            pl.sum("is_unassigned").alias("n_unassigned"),
            pl.sum("is_assigned").alias("n_assigned"),
        ])
        .collect()
    )
    # Convert Polars -> pandas for a simple merge
    bin_counts_pd = bin_counts.to_pandas()

    # Merge onto result (keeps all bins in result)
    result = cba.merge(
    bin_counts_pd,
    on = ["bin_key"],
    how="left",
    )

    # bins that have no points in parquet will be NaN -> make them 0
    # bins that have no points in parquet will be NaN -> make them 0
    # bins that have no points in parquet will be NaN -> make them 0
    result[["n_unassigned", "n_assigned"]] = result[["n_unassigned", "n_assigned"]].fillna(0).astype(np.int64)

    # densities (count per unit area)
    # guard against area=0 just in case (shouldn’t happen if you filtered area>0)
    result["density_unassigned"] = np.where(result["area"] > 0, result["n_unassigned"] / result["area"], np.nan)
    # result["density_assigned"]   = np.where(result["area"] > 0, result["n_assigned"]   / result["area"], np.nan)
    result['x_centroid'] = result['bin_x'] + 25
    result['y_centroid'] = result['bin_y'] + 25
    result['gridcell_id'] = [f"grid_{cell_id}" for cell_id in result['bin_key']]
    result = result[result['area'] > 25 ]

    result.to_parquet(out_path,index=None)
    print('saved',out_path)
    X = (
    result
    .groupby(["gridcell_id", "feature_name"], as_index=False)["n_unassigned"]
    .mean()                     # or sum(), depending on biology
    .pivot(
        index="gridcell_id",
        columns="feature_name",
        values="n_unassigned"
    )
    .fillna(0)
    )

    adata = ad.AnnData(X)

    obs = (
        result.drop_duplicates("gridcell_id")
        .set_index("gridcell_id")[["x_centroid", "y_centroid", "area"]]
    )

    adata.obs = obs.loc[adata.obs_names]
    adata.obs['sample'] = pfile
    adata.obs['gridcell_id'] = adata.obs.index
    adata.obs_names = [f"{pfile}_{cell_id}" for cell_id in adata.obs['gridcell_id']]
    adata.obs['gridcell_id'] = adata.obs_names
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
    adata.layers["counts"] = adata.X.copy()

    adata.write(f"data/unassigned_{pfile}_forMMC.h5ad")

    adatas.append(adata)

    sto = time.time()
    ct = round(sto-sti,6)
    print(ct,'seconds')

adata = adatas[0].concatenate(adatas[1:], index_unique=None)

sc.pp.calculate_qc_metrics(adata,  percent_top=(10, 20, 50, 150), inplace=True)

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=10e4, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps = 50)

sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist = 1)

sc.tl.leiden(adata, resolution = 1)

if not os.path.exists(save_path):
   os.makedirs(save_path)

adata.write_h5ad(f'{save_path}AD_unassigned_adata.h5ad')


print('Analysis Done')