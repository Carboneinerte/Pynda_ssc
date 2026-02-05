### import necessary libraries
from datetime import datetime
import geopandas as gpd
from IPython.display import display
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns
import pytz
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def umap_plot_indi_multi(adata_to_plot,
                         name_dir : str,
                         dir_notebook : str,
                         cluster_to_use : str = 'cell_type_newnum_final',
                         individual_plot : bool = True,
                         save_plot : bool = False,
                         cmap_ : str = 'hls',
                        ):
    '''
    Plot UMAP with all samples combined or with individual samples.
    '''

    adata_to_plot.obsm['umap'] = adata_to_plot.obsm['reduced_pc_20_umap']
    adata_to_plot.obs['umap-1'] = adata_to_plot.obsm['umap'][:, 0]
    adata_to_plot.obs['umap-2'] = adata_to_plot.obsm['umap'][:, 1]
    adata_to_plot.obs['umap-3'] = adata_to_plot.obsm['reduced_pc_20_umap'][:, 0]
    adata_to_plot.obs['umap-4'] = adata_to_plot.obsm['reduced_pc_20_umap'][:, 1]

    adata_to_plot.obs[cluster_to_use] = adata_to_plot.obs[cluster_to_use].astype(str)
    num_clusters = len(adata_to_plot.obs[cluster_to_use].astype(int).unique())
    palette = sns.color_palette(cmap_, n_colors=num_clusters)
    adata_to_plot.obs['leiden_colors'] = adata_to_plot.obs[cluster_to_use].astype(int).apply(lambda x: palette[x])

    if individual_plot == True:
        ## Draw UMAP
        b = int(adata_to_plot.obs['sample'].nunique() / 2)
        _, axs = plt.subplots(b,2,
                                figsize=(15,25)
    
                                )
        axs = axs.flatten()

        def plot_umap(adata_to_plot, color_column, ax, title=None):
            scatter = ax.scatter(adata_to_plot.obs['umap-3'], adata_to_plot.obs['umap-4'], c=adata_to_plot.obs[color_column], s=0.05, alpha=0.8)
            ax.set_title(title)
            ax.axis('off')
        samples_ids = adata_to_plot.obs['sample'].unique()
        for i, sample in enumerate(samples_ids):
            sample_data = adata_to_plot[adata_to_plot.obs['sample'] == sample]
            plot_umap(sample_data, 'leiden_colors', axs[i], title=f"UMAP for {sample}")
            cluster_centroids = sample_data.obs.groupby(cluster_to_use)[['umap-3', 'umap-4']].median()
            
            for cluster_id, centroid in cluster_centroids.iterrows():
                axs[i].text(centroid['umap-3'], centroid['umap-4'], str(cluster_id), color='black', fontsize=12, ha = 'center')

        plt.show()
        
        if save_plot == True:
            plt.savefig(f"{dir_notebook}/plot/{name_dir}/{name_dir}_UMAP_{cluster_to_use}.png", dpi = 300, transparent = True)
        
    ####
    else:
        cell_type_unique = adata_to_plot.obs[cluster_to_use].unique()
        cluster_centroids = adata_to_plot.obs.groupby(cluster_to_use)[['umap-3', 'umap-4']].median()

        # Map each 'leiden' value to a color
        adata_to_plot.obs['leiden_colors'] = adata_to_plot.obs[cluster_to_use].astype(int).apply(lambda x: palette[x])

        fig, ax = plt.subplots(figsize=(15, 10))
        for idx, celltype in enumerate(cell_type_unique):
            adata_sel = adata_to_plot[(adata_to_plot.obs[cluster_to_use] == celltype)]
            if cluster_to_use == 'cell_type_newnum_final':
                celltype_name = adata_sel.obs['cell_type_final'].unique()[0]
            elif cluster_to_use == 'cell_type_newnum_auto_sub':
                celltype_name = adata_sel.obs['cell_type_auto_sub'].unique()[0]
            elif cluster_to_use == 'cell_class_newnum':
                celltype_name = adata_sel.obs['cell_class'].unique()[0]
            elif cluster_to_use == 'region_automap_num':
                celltype_name = adata_sel.obs['region_automap_name'].unique()[0]
            elif (cluster_to_use == 'leiden') or (cluster_to_use == 'kmeans'):
                celltype_name = 'leiden'
            celltype_combine = str(celltype) + '_' + celltype_name
            scat = ax.scatter(adata_sel.obs['umap-3'].values, adata_sel.obs['umap-4'].values, c=adata_sel.obs['leiden_colors'], s = 0.01, label = celltype_combine)
        for cluster_id, centroid in cluster_centroids.iterrows():
            ax.text(centroid['umap-3'], centroid['umap-4'], str(cluster_id), color='black', fontsize=12, ha = 'center')

        plt.legend(markerscale=20, scatterpoints=1000, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        
        if save_plot == True:
            plt.savefig(f"{dir_notebook}/plot/{name_dir}/{name_dir}_UMAP_all.png", dpi = 300, transparent = True)
    

def cluster_plot(adata_to_plot,
                 name_dir,
                 dir_notebook,
                 cluster_to_use : str = 'cell_type_newnum_final',
                 cluster_to_map : list = ['all'],
                 cmap_ : str = 'tab20b',
                 save_plot : bool = False
                 ):
    
    '''
    Spatial plot of individual cells for each sample separately.

    Cluster usable:
    'cell_type_newnum_auto_sub','cell_type_newnum_auto','cell_type_newnum_final',
    'cell_class_newnum','region_automap_num',"leiden","kmeans","circascore"

    Docstring for cluster_plot
    
    :param adata_to_plot: Description
    :param name_dir: Description
    :param dir_notebook: Description
    :param cluster_to_use: Description
    :type cluster_to_use: str
    :param cluster_to_map: Description
    :type cluster_to_map: list
    :param cmap_: Description
    :type cmap_: str
    :param save_plot: Description
    :type save_plot: bool
    '''

    label_to_use = cluster_to_use
    test_dict = {
    'cell_type_newnum_auto_sub':'cell_type_auto_sub',
    'cell_type_newnum_auto':'cell_type_auto',
    'cell_type_newnum_final':'cell_type_final',
    'cell_class_newnum': 'cell_class',
    'region_automap_num':'region_automap_name',
    "leiden":"leiden",
    "kmeans":"kmeans",
    "circascore":"circascore",
    }
    if (cluster_to_use not in test_dict):
        print('Unsupported cluster')
    else:
    
        for cluster, label in test_dict.items():
            label_to_use = label_to_use.replace(cluster, label)
        label_to_use

        adata_to_plot.obs['x_centroid'].astype('float')
        adata_to_plot.obs['y_centroid'].astype('float')
        ### Generate a color palette for the clusters - to make color stay consistent across samples
        adata_to_plot.obs[cluster_to_use] = adata_to_plot.obs[cluster_to_use].astype(str)
        num_clusters = len(adata_to_plot.obs[cluster_to_use].astype(int).unique())
        palette = sns.color_palette(cmap_, n_colors=num_clusters +1)
        adata_to_plot.obs['leiden_colors'] = adata_to_plot.obs[cluster_to_use].astype(int).apply(lambda x: palette[x])

        # Map all cells
        b = int(adata_to_plot.obs['sample'].nunique() / 3)
        fig, axs = plt.subplots(b,3,
                                figsize=(15,25))
        axs = axs.flatten()# Mapping of clusters

        if cluster_to_map != ['all']:
            cluster_to_map2 = cluster_to_map
            color_samples = ['red','green','blue',"black",'magenta','pink',"darkgreen",'coral','orchid','pink']
            while len(color_samples) < len(cluster_to_map):
                color_samples.extend(color_samples)
        
            clusters_plot = {}
            for l in range(0, len(cluster_to_map)):
                dict_temp = {cluster_to_map[l]:color_samples[l]}
                clusters_plot.update(dict_temp)    

        samples_ids = adata_to_plot.obs['sample'].unique()
        for idx, sample in enumerate(samples_ids):
            adata_sel = adata_to_plot[(adata_to_plot.obs['sample'] == sample)]
            cluster_to_map2 = adata_sel.obs[cluster_to_use].unique()
            for cluster_id in cluster_to_map2:
                cluster_data = adata_sel.obs[adata_sel.obs[cluster_to_use] == cluster_id]
                if cluster_to_map != 'all':
                    colors = clusters_plot[cluster_id] if cluster_id in clusters_plot else "none" ### for selected clusters in cluster_plot
                else:
                    colors = cluster_data['leiden_colors'].unique()[0] ### for all clusters
                axs[idx].scatter(cluster_data['x_centroid'], cluster_data['y_centroid'], color=colors, s=0.005, label=cluster_data[label_to_use].unique()[0])
                axs[idx].set_title(f"Sample {sample}")
                axs[idx].xaxis.set_tick_params(labelbottom=False)
                axs[idx].yaxis.set_tick_params(labelleft=False)
                axs[idx].set_aspect('equal', adjustable='box')
        plt.legend(markerscale=50, scatterpoints=1000, bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

        if save_plot == True:
            plt.savefig(f"{dir_notebook}/plot/{name_dir}/{name_dir}_map_{cluster_to_use}.png", dpi = 300, transparent = True)


def polygonplot_dataprep(adata_main,
                         sample_to_plot : str,
                         dir_notebook : str,
                         cluster_to_use : str = 'cell_type_newnum_final',
                         cmap_ : str = 'tab20b'
                         ):
    
    '''
    Docstring for polygonplot_dataprep
    
    Prepare data for polygon plot. Only ONE SAMPLE at the time.


    :param adata_main: main AnnData object
    :param sample_to_plot: Full name from 'sample' columns
    :type sample_to_plot: str
    :param dir_notebook: Description
    :type dir_notebook: str
    :param cluster_to_use: Description
    :type cluster_to_use: str
    :param cmap_: Will be used for the polygons colors
    :type cmap_: str
    '''

    ### Generate a color palette for the clusters - to make color stay consistent across samples
    num_clusters = len(adata_main.obs[cluster_to_use].astype(int).unique())
    palette = sns.color_palette(cmap_, n_colors=num_clusters)
    adata_main.obs['leiden_colors'] = adata_main.obs[cluster_to_use].astype(int).apply(lambda x: palette[x])

    all_samples = np.array(adata_main.obs['sample'].unique())
    sample_position = np.where(all_samples == sample_to_plot)
    sample_position = sample_position[0][0]

    adata_plot = adata_main[adata_main.obs['sample']==all_samples[sample_position]]
    
    cells_geo = gpd.read_file(f'{dir_notebook}/coordinates/polygons/{all_samples[sample_position]}_cells.geojson')
    cells_geo['centroid'] = cells_geo['geometry'].centroid
    cells_geo['x_coor'] = cells_geo['centroid'].x
    cells_geo['y_coor'] = cells_geo['centroid'].y

    if 'objectType' in cells_geo.columns:
        cells_geo = cells_geo[cells_geo['objectType']=='cell']


    # cluster_dict_region = dict(zip(adata_main.obs['cell_id'], adata_main.obs['region_manual_name']))
    cluster_dict_region_a = dict(zip(adata_main.obs['cell_id'], adata_main.obs['region_automap_name']))
    cluster_dict_leiden = dict(zip(adata_main.obs['cell_id'], adata_main.obs['leiden_colors']))
    # cluster_dict = dict(zip(adata_main.obs['cell_id'], adata_main.obs['cell_type_newnum_final']))
    cluster_dict_type = dict(zip(adata_main.obs['cell_id'], adata_main.obs['cell_type_final']))

    if 'circascore' in adata_main.obs.columns:
        cluster_dict_circascore = dict(zip(adata_main.obs['cell_id'], adata_main.obs['circascore']))
        cells_geo['circascore'] = cells_geo['cell'].map(cluster_dict_circascore)

    cells_geo['leiden_colors'] = cells_geo['cell'].map(cluster_dict_leiden)
    cells_geo['cell type'] = cells_geo['cell'].map(cluster_dict_type)
    # cells_geo['region_manual_name'] = cells_geo['cell'].map(cluster_dict_region)
    cells_geo['region_automap_name'] = cells_geo['cell'].map(cluster_dict_region_a)
    cells_geo = cells_geo.dropna(subset=['region_automap_name'])

    df = pd.DataFrame(data=adata_plot.X.toarray(), index=adata_plot.obs_names, columns=adata_plot.var_names)
    df['cell_id'] = df.index
    mapping_dict_region = dict(zip(adata_plot.obs['cell_id'], adata_plot.obs['region_automap_name']))
    mapping_dict_celltype = dict(zip(adata_plot.obs['cell_id'], adata_plot.obs['cell_type_final']))
    mapping_dict_manos = dict(zip(adata_plot.obs['cell_id'], adata_plot.obs['sample']))

    # # # Use .map() function to rename cell contents in 'col1' based on mapping dictionary
    df['region_automap'] = df['cell_id'].map(mapping_dict_region)
    df['cell_type_final'] = df['cell_id'].map(mapping_dict_celltype)
    df['sample'] = df['cell_id'].map(mapping_dict_manos)
    df.dropna(subset=['cell_type_final'], inplace=True)

    return df, cells_geo, cluster_to_use



def polygonplot_plot(df,
                     cells_geo,
                     cluster_to_use : str = 'cell_type_newnum_final',
                     gene_ : str = None,
                     region_ : str = None,
                     region_only : str = None,
                     coord_ : list = None,
                     save_plot : bool = False,
                     legend : bool = False
                     ):
    
    '''
    Docstring for polygonplot_plot
    
    :param df: Description
    :type df: Dataframe
    :param cells_geo: Description
    :type cells_geo: Dataframe
    :param cluster_to_use: Description
    :type cluster_to_use: str
    :param gene_: Description
    :type gene_: str
    :param region_: Description
    :type region_: str
    :param region_only: Description
    :type region_only: str
    :param coord_: Description
    :type coord_: list
    :param save_plot: Description
    :type save_plot: bool
    :param legend: Description
    :type legend: bool
    '''

    if gene_ != None:
        df_dict = dict(zip(df.index, df[gene_]))
        cells_geo[gene_] = cells_geo['cell'].map(df_dict)
    
    if region_only != None:
        cells_geo = cells_geo[cells_geo['region_automap_name']== region_only]
        xmin = cells_geo['x_coor'].min()
        xmax = cells_geo['x_coor'].max()
        ymin = cells_geo['y_coor'].min()
        ymax = cells_geo['y_coor'].max()

    elif region_ != None:
        xmin = cells_geo[cells_geo['region_automap_name']== region_]['x_coor'].min()
        xmax = cells_geo[cells_geo['region_automap_name']== region_]['x_coor'].max()
        ymin = cells_geo[cells_geo['region_automap_name']== region_]['y_coor'].min()
        ymax = cells_geo[cells_geo['region_automap_name']== region_]['y_coor'].max()
        
    elif coord_ != None:
        xmin, xmax, ymin, ymax = coord_

    else:
        print(region_)
        xmin = cells_geo['x_coor'].min()
        xmax = cells_geo['x_coor'].max()
        ymin = cells_geo['y_coor'].min()
        ymax = cells_geo['y_coor'].max()

    cells_geo_crop = cells_geo[(cells_geo['x_coor'] >= xmin) & (cells_geo['x_coor'] <= xmax)  
                            & (cells_geo['y_coor'] >= ymin) & (cells_geo['y_coor'] <= ymax)]

    all_cell_type = cells_geo_crop[cluster_to_use].unique()
    list_cell_nb = range(0, len(all_cell_type))
    mapping_dict = dict(zip(all_cell_type,list_cell_nb))
    cells_geo_crop['cell_type_new'] = cells_geo_crop[cluster_to_use].map(mapping_dict)

    # ##### Extract unique pairs of 'cell type' and 'leiden_colors'
    unique_cell_types = cells_geo_crop[[cluster_to_use, 'leiden_colors']].drop_duplicates()

    ##### Create a color map for the legend
    legend_patches = [
        mpatches.Patch(color=row['leiden_colors'], label=row[cluster_to_use]) 
        for _, row in unique_cell_types.iterrows()
    ]

    fig, ax = plt.subplots(
        figsize=(20,20)
    )
    cells_geo_crop.plot(ax=ax,
                        color = cells_geo_crop['leiden_colors'],
                        alpha=1,
                        aspect=1,
                        zorder=1,
                        edgecolor = "black", #cells_geo_crop['leiden_colors'],
                        linewidth= 1,
                    )
    ax.set_aspect('equal', adjustable='box')

    if gene_ != None:
        cells_geo_crop_bis = cells_geo_crop[cells_geo_crop[gene_] > 0]
        scatter = cells_geo_crop_bis.plot(ax=ax, kind="scatter", x="x_coor", y="y_coor",
                    # color = 'black',
                    c = cells_geo_crop_bis[gene_],
                    cmap = 'inferno',
                    marker = 'v',
                    zorder=3,
                    s = cells_geo_crop_bis[gene_]*25,
                    vmin=0.5,
                    legend=False,
                    edgecolor="black",
                    linewidth = 0.5,
                    # legend_kwds={"label": "Population in 2010", "orientation": "horizontal"},
        )
   


        if len(fig.axes) > 1:
            fig.delaxes(fig.axes[-1]) 

        cbar = fig.colorbar(
            scatter.collections[1],  # Pass the ScalarMappable from the scatter plot
            ax=ax,
            orientation='horizontal',
            fraction=0.06,  # Fraction of original axes size (height adjustment)
            pad=0,         # Padding between plot and colorbar
        )
        cbar.ax.set_position([
            ax.get_position().x0,                # Left coordinate matches the plot
            ax.get_position().y0-0.04,        # Lower the colorbar below the plot
            ax.get_position().width,            # Width matches the plot
            0.04,                               # Height of the colorbar (adjust this value)
        ])

        cbar.set_label(gene_, size=20)

    ##### Add the custom legend
    if legend:    
        ax.legend(handles=legend_patches,     loc='center left', 
            bbox_to_anchor=(1, 0.5), title='Cell Type')

    if save_plot:
        now_ = datetime.now(pytz.timezone('America/Los_Angeles'))
        tod_ = f'{now_.year}_{now_.month}_{now_.day}'
        plt.savefig(f"plot/{tod_}_plot_{region_}_{gene_}.png", format ="png", dpi=600, transparent=True)
    
    plt.show()



def polygonplot_plot_gradient(df, cells_geo, gene_ = None, region_ = None, region_only = None,
                               coord_ = None, cmap_ = 'inferno', save_plot = False):
    if gene_ != None:
        df_dict = dict(zip(df.index, df[gene_]))
        cells_geo[gene_] = cells_geo['cell'].map(df_dict)
    
    if region_only != None:
        cells_geo = cells_geo[cells_geo['region_automap_name']== region_only]
        xmin = cells_geo['x_coor'].min()
        xmax = cells_geo['x_coor'].max()
        ymin = cells_geo['y_coor'].min()
        ymax = cells_geo['y_coor'].max()

    elif region_ != None:
        xmin = cells_geo[cells_geo['region_automap_name']== region_]['x_coor'].min()
        xmax = cells_geo[cells_geo['region_automap_name']== region_]['x_coor'].max()
        ymin = cells_geo[cells_geo['region_automap_name']== region_]['y_coor'].min()
        ymax = cells_geo[cells_geo['region_automap_name']== region_]['y_coor'].max()
        
    elif coord_ != None:
        xmin, xmax, ymin, ymax = coord_

    else:
        print(region_)
        xmin = cells_geo['x_coor'].min()
        xmax = cells_geo['x_coor'].max()
        ymin = cells_geo['y_coor'].min()
        ymax = cells_geo['y_coor'].max()

    cells_geo_crop = cells_geo[(cells_geo['x_coor'] >= xmin) & (cells_geo['x_coor'] <= xmax)  
                            & (cells_geo['y_coor'] >= ymin) & (cells_geo['y_coor'] <= ymax)]

    all_cell_type = cells_geo_crop['cell type'].unique()
    list_cell_nb = range(0, len(all_cell_type))
    mapping_dict = dict(zip(all_cell_type,list_cell_nb))
    cells_geo_crop['cell_type_new'] = cells_geo_crop['cell type'].map(mapping_dict)

    # ##### Extract unique pairs of 'cell type' and 'leiden_colors'
    unique_cell_types = cells_geo_crop[['cell type', 'leiden_colors']].drop_duplicates()

    ##### Create a color map for the legend
    legend_patches = [
        mpatches.Patch(color=row['leiden_colors'], label=row['cell type']) 
        for _, row in unique_cell_types.iterrows()
    ]

    fig, ax = plt.subplots(
        figsize=(20,20)
    )

    cells_geo_crop.plot(ax=ax,
                    column = cells_geo_crop[gene_], 
                    cmap = cmap_, vmin = 0.25,
                    alpha=1,
                    aspect=1,
                    # edgecolor =cells_geo_crop[gene_],
                   )
    ax.set_aspect('equal', adjustable='box')


    # if len(fig.axes) > 1:
    #     fig.delaxes(fig.axes[-1]) 

    # cbar = fig.colorbar(
    #     scatter.collections[1],  # Pass the ScalarMappable from the scatter plot
    #     ax=ax,
    #     orientation='horizontal',
    #     fraction=0.06,  # Fraction of original axes size (height adjustment)
    #     pad=0,         # Padding between plot and colorbar
    # )
    # cbar.ax.set_position([
    #     ax.get_position().x0,                # Left coordinate matches the plot
    #     ax.get_position().y0-0.04,        # Lower the colorbar below the plot
    #     ax.get_position().width,            # Width matches the plot
    #     0.04,                               # Height of the colorbar (adjust this value)
    # ])

    # cbar.set_label(gene_, size=20)

    fig.colorbar()

    ##### Add the custom legend
    ax.legend(handles=legend_patches,     loc='center left', 
        bbox_to_anchor=(1, 0.5), title='Cell Type')

    if save_plot == True:
        now_ = datetime.now(pytz.timezone('America/Los_Angeles'))
        tod_ = f'{now_.year}_{now_.month}_{now_.day}'
        plt.savefig(f"plot/{tod_}_plot_{region_}_{gene_}.svg", dpi=300, transparent=True)
    
    plt.show()