import pandas as pd
import geopandas as gpd
import geojson
import argparse
from config_local import dir_raw, dir_processed

parser = argparse.ArgumentParser(
    description='Convert cell boundaries to GeoJSON.'
    )

parser.add_argument('--sample', type=str, required=True, help='Sample name to process')
parser.add_argument('--name_dir', type=str, help='Name of the experience')
args = parser.parse_args()
sample_id = args.sample
name_dir = args.name_dir


def cell_boundaries_geojson(dir_raw, dir_processed, sample_id, name_dir):
    cells = pd.read_csv(f"{dir_raw}/{name_dir}/{sample_id}/cell_boundaries.csv.gz")
    
    cells.groupby("cell_id")[['vertex_x', 'vertex_y']].agg({"vertex_x": list, "vertex_y": list}).reset_index().rename(columns={"vertex_x": "xs", "vertex_y": "ys"})
    cells['cell_id'] = sample_id + '_' + cells['cell_id']
    
    # Group the dataframe by the "Selection" column
    grouped = cells.groupby('cell_id')
    features = []
    
    for cell_id, group in grouped:
        coordinates = [(x, y) for x, y in zip(group['vertex_x'], group['vertex_y'])]
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])
        
        polygon = geojson.Polygon([coordinates])
        feature = geojson.Feature(geometry=polygon, properties={"cell": cell_id})
        features.append(feature)
    
    feature_collection = geojson.FeatureCollection(features)
    
    with open(f'{dir_processed}/{name_dir}/polygons/{sample_id}_cells.geojson', 'w') as f:
        geojson.dump(feature_collection, f)
    
    print(f"{sample_id}: GeoJSON saved")

cell_boundaries_geojson(dir_raw, dir_processed, sample_id, name_dir)