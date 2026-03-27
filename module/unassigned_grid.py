import polars as pl
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import box
import time
import anndata as ad
import scanpy as sc
import os
from IPython.display import clear_output
from module.config_local import dir_raw,dir_processed
from module.export import coordinates_to_Geojson


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
    bin_size=50,
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

def unassigned_density(samples_ids: list, name_dir:str):
    BIN = 50   # your x (same for width/height)
    x0 = 0.0          # grid origin; set to 0 or use min below
    y0 = 0.0
    adatas = []
    save_path = f"{dir_processed}/h5ad/{name_dir}_unassigned/"

    for pfile in samples_ids:
        print(pfile, name_dir)
        sti = time.time()
        parquet_in = f'{dir_raw}/{pfile}/transcripts_light.parquet'
        
        # out_path = "data/density-bins-5k/"+pfile+"-densitybins.parquet"
        print(pfile,'bin size = ',BIN)
        #this takes the parquet file and saves a new parquet file with bin columns. you could rename the out file if you want to save it in place
        df = (
            pl.scan_parquet(parquet_in)  # lazy scan (doesn't load all at once)
            .filter(~pl.col("feature_name").str.contains("_"))
                    #should filter out the codewords and low quality counts
            .with_columns([
                (((pl.col("x_location") - x0) / BIN).floor()*BIN).cast(pl.Int64).alias("bin_x"),
                (((pl.col("y_location") - y0) / BIN).floor()*BIN).cast(pl.Int64).alias("bin_y"),
            ])
            # optional single key:
            .with_columns([
                (pl.col("bin_x").cast(pl.Utf8) + pl.lit("_") + pl.col("bin_y").cast(pl.Utf8)).alias("bin_key")
            ])
        )
        print(df.collect_schema().names())
        print('binned transcripts')
        #df.sink_parquet(parquet_out)  # writes efficiently

        #load in shapes. will have to loop through these too
        try:
            dfg_r = gpd.read_file(f"{dir_processed}/coordinates/whole_section/{pfile}_whole_section_annotation.geojson")
        except:
            coordinates_to_Geojson(pfile,dir_raw)
            dfg_r = gpd.read_file(f"{dir_processed}/coordinates/whole_section/{pfile}_whole_section_annotation.geojson")
        # this gets rid of the warnings Only do this if the CRS is currently geographic
        print(dfg_r.head())
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
        print(cba.columns)
        print(df.collect_schema().names())
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
        print(bin_counts_pd.head())

        # Merge onto result (keeps all bins in result)
        result = cba.merge(
        bin_counts_pd,
        on = ["bin_key"],
        how="left",
        )

        # bins that have no points in parquet will be NaN -> make them 0
        result[["n_unassigned", "n_assigned"]] = result[["n_unassigned", "n_assigned"]].fillna(0).astype(np.int64)

        # densities (count per unit area)
        # guard against area=0 just in case (shouldn’t happen if you filtered area>0)
        result["density_unassigned"] = np.where(result["area"] > 0, result["n_unassigned"] / result["area"], np.nan)
        result['x_centroid'] = result['bin_x'] + 25
        result['y_centroid'] = result['bin_y'] + 25
        result['gridcell_id'] = [f"grid_{cell_id}" for cell_id in result['bin_key']]
        result = result[result['area'] > 25 ]
        result = result[(result['n_unassigned']!=0)&(result['n_assigned']!=0)]

        # result.to_parquet(out_path,index=None)
        # print('saved',out_path)
        
        X = (
        result
        .groupby(["gridcell_id", "feature_name"], as_index=False)["n_unassigned"]
        .mean()                    
        .pivot(
            index="gridcell_id",
            columns="feature_name",
            values="n_unassigned"
        )
        .fillna(0)
        )

        adata_temp = ad.AnnData(X)

        obs = (
            result.drop_duplicates("gridcell_id")
            .set_index("gridcell_id")[["x_centroid", "y_centroid", "area"]]
        )

        adata_temp.obs = obs.loc[adata_temp.obs_names]
        adata_temp.obs['sample'] = pfile
        adata_temp.obs['gridcell_id'] = adata_temp.obs.index
        adata_temp.obs_names = [f"{pfile}_{cell_id}" for cell_id in adata_temp.obs['gridcell_id']]
        adata_temp.obs['gridcell_id'] = adata_temp.obs_names
        adata_temp.obsm["spatial"] = adata_temp.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
        adata_temp.layers["counts"] = adata_temp.X.copy()

        try:
            adata_temp.write(f"{save_path}/unassigned_{pfile}_forMMC.h5ad")
        except:
            os.makedirs(save_path)
            adata_temp.write(f"{save_path}/unassigned_{pfile}_forMMC.h5ad")

        adata_temp.write(f"{save_path}/unassigned_{pfile}_forMMC.h5ad")

        adatas.append(adata_temp)

        sto = time.time()
        ct = round(sto-sti,6)
        print(ct,'seconds')
        clear_output()

    adata = adatas[0].concatenate(adatas[1:], index_unique=None)
    print('Analysis Done')

    return adata
    # try:
    #     adata.write_h5ad(f'{save_path}/{name_dir}_unassigned_adata.h5ad')
    # except:
    #     os.makedirs(save_path)
    #     adata.write_h5ad(f'{save_path}/{name_dir}_unassigned_adata.h5ad')
