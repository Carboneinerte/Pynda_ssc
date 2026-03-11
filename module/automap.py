import json
from shapely.geometry import Polygon, mapping, MultiPolygon
from shapely.ops import unary_union
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from sklearn.neighbors import KNeighborsClassifier
from datetime import datetime
from scipy.interpolate import CubicSpline
from skimage import measure
from scipy.interpolate import CubicSpline
from matplotlib.patches import Patch
import json
from shapely.geometry import Polygon, mapping, MultiPolygon
from shapely.ops import unary_union
import numpy as np



# data = df

# Function to print messages with elapsed time since the script started
# start_time = datetime.now()

###

def knn_mst_clustering(data:pd.DataFrame,
                       cell_type, k=10, std_threshold=2, min_component_size=100):
    # Filter the data for the selected cell type
    cell_data = data[data['knn_celltype'] == cell_type]
    coordinates = cell_data[['x_centroid', 'y_centroid']].values

    # Apply K-Nearest Neighbors
    knn = NearestNeighbors(n_neighbors=k)
    knn.fit(coordinates)
    distances, indices = knn.kneighbors(coordinates)

    # Create a graph for this specific cell type
    graph = nx.Graph()

    # Adding nodes with positions
    for i, (x, y) in enumerate(coordinates):
        graph.add_node(i, pos=(x, y))

    # Calculate mean and standard deviation of edge lengths
    all_distances = distances.flatten()
    mean_distance = np.mean(all_distances)
    std_distance = np.std(all_distances)
    threshold = mean_distance + std_threshold * std_distance

    # Adding edges based on nearest neighbors and pruning outliers
    for i, neighbors in enumerate(indices):
        for neighbor in neighbors:
            if i != neighbor:  # Avoid self-loops
                distance = distances[i][np.where(indices[i] == neighbor)[0][0]]
                if distance <= threshold:  # Prune edges longer than the threshold
                    graph.add_edge(i, neighbor, weight=distance)

    # Remove small components
    components = list(nx.connected_components(graph))
    large_components = [comp for comp in components if len(comp) >= min_component_size]

    for component in components:
        if len(component) < min_component_size:
            graph.remove_nodes_from(component)

    # Label nodes based on their cluster number and update the dataframe
    for cluster_num, component in enumerate(large_components, start=1):
        for node in component:
            original_index = cell_data.index[node]
            data.at[original_index, 'cluster_label'] = f"{cell_type}_cluster_{cluster_num}"

###

def assign_chunk_labels(df: pd.DataFrame,
                        chunk_size:int,
                        threshold: float = 0.5
                        ):
    # Multiply centroid coordinates by 10 to convert to pixel coordinates
    
    # Determine the grid size based on the pixel coordinates
    x_min, x_max = 0, df['x_pixel'].max() + 2 * chunk_size
    y_min, y_max = 0, df['y_pixel'].max() + 2 * chunk_size
    
    # Create bins for chunking the data
    x_bins = np.arange(x_min, x_max + chunk_size, chunk_size)
    y_bins = np.arange(y_min, y_max + chunk_size, chunk_size)
    
    # Initialize an array to store chunk labels
    grid = np.full((len(y_bins) - 1, len(x_bins) - 1), np.nan)
    
    # Iterate through each chunk and assign a label based on the majority cell type
    total_chunks = (len(x_bins) - 1) * (len(y_bins) - 1)
    chunk_count = 0
    
    for i in range(len(x_bins) - 1):
        for j in range(len(y_bins) - 1):
            chunk_count += 1
            # Define the bounds of the current chunk
            x_lower, x_upper = x_bins[i], x_bins[i+1]
            y_lower, y_upper = y_bins[j], y_bins[j+1]
            
            # Filter the dataframe to get cells within the current chunk
            chunk = df[(df['x_pixel'] >= x_lower) & (df['x_pixel'] < x_upper) &
                       (df['y_pixel'] >= y_lower) & (df['y_pixel'] < y_upper)]
            
            if len(chunk) > 0:
                # Calculate the proportion of each cell type in the chunk
                counts = chunk['knn_celltype'].value_counts(normalize=True)
                
                if counts.iloc[0] >= threshold:
                    grid[j, i] = counts.idxmax() + 1  # Adjust for 0-indexing
                else:
                    grid[j, i] = np.nan  # Below threshold, consider background
            else:
                grid[j, i] = np.nan  # No cells in chunk, consider background
            
            # Print progress
            print(f'Processing chunk {chunk_count}/{total_chunks} (x: {x_lower}-{x_upper}, y: {y_lower}-{y_upper})', end='\r')
    
    return grid, x_bins, y_bins

def fill_empty_chunks(grid, x_bins, y_bins, n=4):
    # Create a copy of the grid to store the updated values
    filled_grid = grid.copy()
    
    # Get the dimensions of the grid
    rows, cols = grid.shape
    
    # Iterate over the grid, skipping the borders to avoid index errors
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            if np.isnan(grid[i, j]):
                # Check the four neighbors (up, right, down, left)
                neighbors = [
                    grid[i-1, j],  # Up
                    grid[i+1, j],  # Down
                    grid[i, j-1],  # Left
                    grid[i, j+1]   # Right
                ]
                
                # Filter out NaN values from neighbors
                valid_neighbors = [neighbor for neighbor in neighbors if not np.isnan(neighbor)]
                
                # Check if there are at least `n` neighbors with the same value
                if len(valid_neighbors) >= n:
                    most_common_value = max(set(valid_neighbors), key=valid_neighbors.count)
                    if valid_neighbors.count(most_common_value) >= n:
                        filled_grid[i, j] = most_common_value
    
    return filled_grid, x_bins, y_bins

def fill_any_chunks(grid, x_bins, y_bins, n=4):
    # Create a copy of the grid to store the updated values
    filled_grid = grid.copy()
    
    # Get the dimensions of the grid
    rows, cols = grid.shape
    
    # Iterate over the grid, skipping the borders to avoid index errors
    for i in range(1, rows - 1):
        for j in range(1, cols - 1):
            # Check the four neighbors (up, right, down, left)
            neighbors = [
                grid[i-1, j],  # Up
                grid[i+1, j],  # Down
                grid[i, j-1],  # Left
                grid[i, j+1]   # Right
            ]
            
            # Filter out NaN values from neighbors
            valid_neighbors = [neighbor for neighbor in neighbors if not np.isnan(neighbor)]
            
            # Check if there are at least `n` neighbors with the same value
            if len(valid_neighbors) >= n:
                most_common_value = max(set(valid_neighbors), key=valid_neighbors.count)
                if valid_neighbors.count(most_common_value) >= n:
                    filled_grid[i, j] = most_common_value
    
    return filled_grid, x_bins, y_bins

def smooth_boundary(contour, x_bins, y_bins):
    # Convert contour points to the original coordinate space
    x = contour[:, 1]
    y = contour[:, 0]

    # Map grid coordinates to original pixel coordinates
    x_orig = np.interp(x, np.arange(len(x_bins) - 1), x_bins[:-1] + np.diff(x_bins) / 2)
    y_orig = np.interp(y, np.arange(len(y_bins) - 1), y_bins[:-1] + np.diff(y_bins) / 2)
    
    # Parameterization for Catmull-Rom spline
    t = np.arange(len(x_orig))
    
    # Create cubic splines for x and y
    cs_x = CubicSpline(t, x_orig, bc_type='clamped')
    cs_y = CubicSpline(t, y_orig, bc_type='clamped')
    
    # Generate smooth parameterization
    t_new = np.linspace(t.min(), t.max(), 10 * len(t))
    x_smooth = cs_x(t_new)
    y_smooth = cs_y(t_new)
    
    return x_smooth, y_smooth

def plot_cells_and_chunks(df: pd.DataFrame,
                          grid, x_bins, y_bins):
    # Create a scatter plot for all cells
    plt.figure(figsize=(16, 10))
    scatter = plt.scatter(df['x_pixel'], df['y_pixel'], c=df['cell_type_newnum_final'], cmap='tab20', s=0.5, alpha=0.6)
    
    # Get the colormap and normalize based on the number of cell types
    cmap = plt.cm.get_cmap('tab20')
    # Use MinMax normalization that accounts for the full range of knn_celltypes
    norm = plt.Normalize(vmin=df['knn_celltype'].min(), vmax=df['knn_celltype'].max())
    
    # Create a list to store patches for the legend
    legend_patches = []
    
    # Loop through each unique value in the grid (excluding NaN)
    for chunk_value in np.unique(grid[~np.isnan(grid)]):
        print(f"Plotting boundaries for chunk value: {chunk_value}")
        
        # Get the corresponding color for the chunk value (should match knn cell type colors)
        color = cmap(norm(chunk_value))  # Removed the -1 to properly handle 0
        
        # Find contours for the current chunk value
        contours = measure.find_contours(grid == chunk_value, level=0.5)
        
        print(f"Found {len(contours)} contours for chunk value: {chunk_value}")
        # Plot the smoothed contours and fill the interiors
        for i, contour in enumerate(contours):
            print(f"Plotting contour {i} for chunk value: {chunk_value}")
            x_smooth, y_smooth = smooth_boundary(contour, x_bins, y_bins)
            plt.plot(x_smooth, y_smooth, color=color, linewidth=2, alpha=0.5)  # Transparent contours
            
            # Fill the interior of the contour with the same color but more transparent
            plt.fill(x_smooth, y_smooth, color=color, alpha=0.2)
        
        # Find the corresponding label for the knn_celltype from the mapping
        matching_labels = df[df['knn_celltype'] == chunk_value]['knn_celltype_label']
        
        # Only proceed if there is at least one matching label
        if not matching_labels.empty:
            label = matching_labels.iloc[0]
            # Add patch for legend with the new label
            legend_patches.append(Patch(color=color, label=f"{label} (KNN {int(chunk_value)})"))
    
    # Add custom legend with knn cell type labels and their corresponding colors
    plt.legend(handles=legend_patches, title="KNN Cell Type Labels", loc='upper right', bbox_to_anchor=(1.15, 1))
    
    # Set plot labels and title
    plt.xlabel('X Pixel Coordinate')
    plt.ylabel('Y Pixel Coordinate')
    plt.title('Xenium Spatial Transcriptomics - Cells and Inferred Regions')
    plt.show()

def grid_to_geojson_with_scaling(grid, x_bins, y_bins, mapping_dict,
                                 scale_factor: int =10):
    # Dictionary to store polygons grouped by cell_type_newnum
    polygons_by_type = {}
    
    rows, cols = grid.shape
    
    for i in range(rows):
        for j in range(cols):
            if not np.isnan(grid[i, j]):
                # Map the grid cell value to the original cell_type_newnum
                cell_type_newnum = mapping_dict[int(grid[i, j])]
                
                # Define the corners of the polygon (rectangular area) for this grid cell
                # Scale the coordinates by the given scale factor
                polygon = Polygon([
                    (x_bins[j] / scale_factor, y_bins[i] / scale_factor),         # Bottom-left
                    (x_bins[j+1] / scale_factor, y_bins[i] / scale_factor),       # Bottom-right
                    (x_bins[j+1] / scale_factor, y_bins[i+1] / scale_factor),     # Top-right
                    (x_bins[j] / scale_factor, y_bins[i+1] / scale_factor),       # Top-left
                    (x_bins[j] / scale_factor, y_bins[i] / scale_factor)          # Closing the polygon
                ])
                
                # Add the polygon to the corresponding cell_type_newnum group
                if cell_type_newnum not in polygons_by_type:
                    polygons_by_type[cell_type_newnum] = []
                polygons_by_type[cell_type_newnum].append(polygon)
    
    # Merge polygons for each cell_type_newnum
    features = []
    for cell_type_newnum, polygons in polygons_by_type.items():
        # Combine all polygons into a single MultiPolygon and then merge them into one Polygon
        merged_polygon = unary_union(MultiPolygon(polygons))
        
        # Create a feature with the merged polygon and the corresponding cell_type_newnum
        feature = {
            "type": "Feature",
            "geometry": mapping(merged_polygon),
            "properties": {
                "cell_type_newnum_final": cell_type_newnum  # Store the cell type_newnum
            }
        }
        
        features.append(feature)
    
    # Create the GeoJSON structure
    geojson = {
        "type": "FeatureCollection",
        "features": features
    }
    
    return geojson