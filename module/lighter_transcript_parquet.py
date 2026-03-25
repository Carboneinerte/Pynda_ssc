import pandas as pd
import geopandas as gpd
from shapely.geometry import box
from module.config_local import dir_raw

def lighter_transcript_parquet_resesg(samples_ids:list,
                               gene_list:list,
                               dir_raw:str = dir_raw,
                               ):

    for sample in samples_ids:
        print(f"{sample}")
        df_cell = gpd.read_file(f'{dir_raw}/{sample}/{sample}.geojson')
        df_cell = df_cell[df_cell['objectType']=='cell']
        
        df_trans = pd.read_parquet(f'{dir_raw}/{sample}/transcripts.parquet',
                                columns=['transcript_id','x_location','y_location','qv',"feature_name"],
                                filters=[('qv','>=',20), ("feature_name", "in", gene_list)])
        print("import done")             
        df_trans['x_loc'] = df_trans['x_location']/0.2125
        df_trans['y_loc'] = df_trans['y_location']/0.2125
        
        df_trans_gdf = gpd.GeoDataFrame(
            df_trans,
            geometry=gpd.points_from_xy(df_trans.x_loc, df_trans.y_loc)
        )
        
        df_trans_gdf.crs = 'EPSG:4326'
        df_cell.crs = 'EPSG:4326'
        
        print("Geopandas conversion done, start matching")

        matched_cells = gpd.sjoin(
            df_trans_gdf,
            df_cell,
            predicate='within',
            how='left'
        )
        print('matching done, start updating parquet')
        matched_cells['id'] = matched_cells['id'].fillna('UNASSIGNED')

        dict_resegment = dict(zip(matched_cells['transcript_id'], matched_cells['id']))
        df_trans['cell_id'] = df_trans['transcript_id'].map(dict_resegment)
        print("update done, start saving")
        df_trans.to_parquet(f'{dir_raw}/{sample}/transcripts_light.parquet')
        print(" ")

def lighter_transcript_parquet(samples_ids:list,
                               gene_list:list,
                               dir_raw:str = dir_raw,
                               
                               ):

    for sample in samples_ids:
        print(f"{sample}")


        df_trans = pd.read_parquet(f'{dir_raw}/{sample}/transcripts.parquet',
                                columns=['cell_id','feature_name','x_location','y_location',"qv"],
                                filters=[('qv','>=',20), ("feature_name", "in", gene_list)])
        

        df_trans.to_parquet(f'{dir_raw}/{sample}/transcripts_light.parquet')
        print(" ")