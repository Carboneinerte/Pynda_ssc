import os
import pandas as pd
import geojson
import geopandas as gpd

from module.config_local import dir_processed, dir_raw



def export_XE(df:pd.DataFrame,
              sample:str,
              group:str,
              name_dir:str
              ):
    
    df_temp = df[df['sample']==sample]
    df_temp['cell_id'] = df_temp.index
    df_temp['cell_id'] = df['cell_id'].apply(lambda x: x.split("_")[1])
    df_temp['group'] = df[group]
    df_temp = df_temp.filter(['cell_id','group'],axis=1)
    try:
        df_temp.to_csv(f'{dir_processed}/csv/{name_dir}/XE_export/{sample}_{group}_export_XE.csv', index=None)
    except:
        os.makedirs(f'{dir_processed}/csv/{name_dir}/XE_export/')
        df_temp.to_csv(f'{dir_processed}/csv/{name_dir}/XE_export/{sample}_{group}_export_XE.csv', index=None)



def coordinates_to_Geojson(sample:str,
                           dir_raw:str = dir_raw,
                           dir_processed:str = dir_processed):
    

    BR_df = pd.read_csv(f"{dir_raw}/{sample}/{sample}_whole_section_annotation.csv", comment = "#")

    # Group the dataframe by the "Selection" column
    grouped = BR_df.groupby('Selection')

    # List to hold GeoJSON features
    features = []

    for name, group in grouped:
        # Create a list of coordinates for each region
        coordinates = [(x, y) for x, y in zip(group['X'], group['Y'])]
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])
        
        # Create a GeoJSON polygon for the region
        polygon = geojson.Polygon([coordinates])
        feature = geojson.Feature(geometry=polygon, properties={"region": name})
        features.append(feature)

    # Create a GeoJSON FeatureCollection
    feature_collection = geojson.FeatureCollection(features)

    feature_collection['features'][0]['properties']['region'] = 'wholebrain'

    # Save the GeoJSON FeatureCollection to a file
    try:
        with open(f'{dir_processed}/coordinates/whole_section/{sample}_whole_section_annotation.geojson', 'w') as f:
            geojson.dump(feature_collection, f)
    except:
        os.makedirs(f'{dir_processed}/coordinates/whole_section/')
        with open(f'{dir_processed}/coordinates/whole_section/{sample}_whole_section_annotation.geojson', 'w') as f:
            geojson.dump(feature_collection, f)

    print("GeoJSON saved")