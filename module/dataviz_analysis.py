### import necessary libraries
import os
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')
import progressbar
import geopandas as gpd
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from IPython.display import clear_output
import matplotlib.patches as mpatches
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import random
import scanpy as sc


import warnings

from anndata import ImplicitModificationWarning
warnings.simplefilter('ignore', ImplicitModificationWarning)

from module.misc import save_figure
from module.config_local import dir_processed

from pathlib import Path
import anndata as ad
import holoviews as hv
import panel as pn
import hvplot.pandas    # noqa
import numpy as np
import pooch

import scanpy as sc

import hv_anndata
from hv_anndata import *

hv_anndata.register()
hv.extension("bokeh")
pn.extension()
pn.config.throttled = True

def umap_plot_indi_multi(adata_to_plot: sc.AnnData,
                         name_dir : str,
                         dir_processed : str = dir_processed,
                         cluster_to_use : str = 'cell_type_newnum_final',
                         individual_plot : bool = True,
                         save_plot : bool = False,
                         cmap_ : str = 'hls',
                         size:float = 0.05,
                        ):
    '''
    Plot UMAP with all samples combined or with individual samples.
    '''

    adata_to_plot.obsm['umap'] = adata_to_plot.obsm['X_umap']
    adata_to_plot.obs['umap-1'] = adata_to_plot.obsm['umap'][:, 0]
    adata_to_plot.obs['umap-2'] = adata_to_plot.obsm['umap'][:, 1]
    adata_to_plot.obs['umap-3'] = adata_to_plot.obsm['umap'][:, 0]
    adata_to_plot.obs['umap-4'] = adata_to_plot.obsm['umap'][:, 1]

    adata_to_plot.obs[cluster_to_use] = adata_to_plot.obs[cluster_to_use].astype(str)
    num_clusters = len(adata_to_plot.obs[cluster_to_use].astype(int).unique())
    palette = sns.color_palette(cmap_, n_colors=num_clusters)
    adata_to_plot.obs['leiden_colors'] = adata_to_plot.obs[cluster_to_use].astype(int).apply(lambda x: palette[x])

    if individual_plot == True:
        ## Draw UMAP
        b = int(adata_to_plot.obs['sample'].nunique() / 2)
        fig, axs = plt.subplots(b,2,
                                figsize=(15,25)
    
                                )
        axs = axs.flatten()

        def plot_umap(adata_to_plot, color_column, ax, title=None):
            scatter = ax.scatter(adata_to_plot.obs['umap-3'], adata_to_plot.obs['umap-4'], c=adata_to_plot.obs[color_column], s=size, alpha=0.8)
            ax.set_title(title)
            ax.axis('off')
        samples_ids = adata_to_plot.obs['sample'].unique()
        for i, sample in enumerate(samples_ids):
            sample_data = adata_to_plot[adata_to_plot.obs['sample'] == sample]
            plot_umap(sample_data, 'leiden_colors', axs[i], title=f"UMAP for {sample}")
            cluster_centroids = sample_data.obs.groupby(cluster_to_use)[['umap-3', 'umap-4']].median()
            
            for cluster_id, centroid in cluster_centroids.iterrows():
                axs[i].text(centroid['umap-3'], centroid['umap-4'], str(cluster_id), color='black', fontsize=12, ha = 'center')
                axs[i].set_aspect('equal', adjustable='box')

        
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.3)

        plt.show()
        
        if save_plot == True:
            suffix_save = f'UMAP_indi_{cluster_to_use}' 
            save_figure(fig, suffix_save, name_dir, format='png')

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
        plt.show()
        if save_plot == True:
            suffix_save = f'UMAP_all_{cluster_to_use}' 
            save_figure(fig, suffix_save, name_dir, format='png')

def cluster_plot(adata_to_plot,
                 name_dir,
                 dir_processed: str = dir_processed,
                 cluster_to_use : str = 'cell_type_newnum_final',
                 cluster_to_map : list = ['all'],
                 cmap_ : str = 'tab20b',
                 save_plot : bool = False,
                 size: float = 0.01,
                 ):
    
    '''
    Spatial plot of individual cells for each sample separately.

    Cluster usable:
    'cell_type_newnum_auto_sub','cell_type_newnum_auto','cell_type_newnum_final',
    'cell_class_newnum','region_automap_num',"leiden","kmeans","circascore"

    Docstring for cluster_plot
    
    :param adata_to_plot: Description
    :param name_dir: Description
    :param dir_processed: Description
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
    'region_manual_map_num' : 'region_manual_map',
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
        samples_ids = adata_to_plot.obs['sample'].unique().sort_values()

        # Map all cells
        b = int(adata_to_plot.obs['sample'].nunique() / 3)
        if b == 0:
            b=1
        if len(samples_ids) <= 3:
            a = len(samples_ids)
            b = 1
        else:
            a=3
            
        fig, axs = plt.subplots(b,a,
                                figsize=(15,6))
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

        
        for idx, sample in enumerate(samples_ids):
            adata_sel = adata_to_plot[(adata_to_plot.obs['sample'] == sample)]
            cluster_to_map2 = adata_sel.obs[cluster_to_use].unique()
            for cluster_id in cluster_to_map2:
                cluster_data = adata_sel.obs[adata_sel.obs[cluster_to_use] == cluster_id]
                if cluster_to_map != ['all']:
                    colors = clusters_plot[cluster_id] if cluster_id in clusters_plot else "none" ### for selected clusters in cluster_plot
                else:
                    colors = cluster_data['leiden_colors'].unique()[0] ### for all clusters
                axs[idx].scatter(cluster_data['x_centroid'], cluster_data['y_centroid'], color=colors, s=size, marker = 'o', label=cluster_data[label_to_use].unique()[0])
                axs[idx].set_title(f"Sample {sample}")
                axs[idx].xaxis.set_tick_params(labelbottom=False)
                axs[idx].yaxis.set_tick_params(labelleft=False)
                axs[idx].set_aspect('equal', adjustable='box')
        # plt.legend(markerscale=5, scatterpoints=10,fontsize='small', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, mode="expand", ncol =3)
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1, hspace=0.3)
        plt.show()
        
        if save_plot == True:
            suffix_save = f'clusterplot_{cluster_to_use}' 
            save_figure(fig, suffix_save, name_dir, format='png')

def polygonplot_dataprep(adata_main: sc.AnnData,
                         sample_to_plot : str,
                         dir_processed : str = dir_processed,
                         cluster_to_use : str = 'cell_type_newnum_final',
                         cmap_ : str = 'tab20b'
                         ):
    
    '''
    Docstring for polygonplot_dataprep
    
    Prepare data for polygon plot. Only ONE SAMPLE at the time.


    :param adata_main: main AnnData object
    :param sample_to_plot: Full name from 'sample' columns
    :type sample_to_plot: str
    :param dir_processed: Description
    :type dir_processed: str
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
    
    cells_geo = gpd.read_file(f'{dir_processed}/coordinates/polygons/{all_samples[sample_position]}_cells.geojson')
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

def polygonplot_plot(df: pd.DataFrame,
                     cells_geo:gpd.GeoDataFrame,
                     name_dir: str,            
                     dir_processed: str = dir_processed,
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

    if save_plot == True:
        suffix_save = f'plot_{region_}_{gene_}' 
        save_figure(fig, suffix_save, name_dir, format='png')

    plt.show()

def polygonplot_plot_gradient(
        df, cells_geo,
        name_dir: str,
        dir_processed: str = dir_processed,
        gene_: str = None,
        region_: str = None,
        region_only: str = None,
        coord_: list = None,
        cmap_:str = 'inferno',
        save_plot: bool = False
        ):
    
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

    # fig.colorbar()

    ##### Add the custom legend
    ax.legend(handles=legend_patches,     loc='center left', 
        bbox_to_anchor=(1, 0.5), title='Cell Type')
    
    if save_plot == True:
        suffix_save = f'plot_{region_}_{gene_}' 
        save_figure(fig, suffix_save, name_dir, format='png')

    plt.show()

def DEG_one_condition(adata: sc.AnnData,
                      name_dir: str,
                      cluster_to_use: str,
                      group_col: str,
                      grp_ctrl: str,
                      filters_bool: bool,
                      filters_dict: dict,
                      dir_processed: str = dir_processed,):

    dfs = []
    dfs_filter = []
    all_groups = np.array(adata.obs[cluster_to_use].unique())
    all_groups_type_sheet = all_groups

    
    bar = progressbar.ProgressBar(maxval=len(all_groups)+1, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    for idx, cell_type_to_extract in enumerate(all_groups):
        adata_microglia = adata[adata.obs[cluster_to_use] == cell_type_to_extract]
        sc.pp.calculate_qc_metrics(adata_microglia, percent_top=[20, 50], inplace=True)
        print(f"Start analysis of {cell_type_to_extract}")
        ### Extract gene expression per cluster + log fold change + p-value
        
        clust_uniq = adata.obs[group_col].unique()
        dat = pd.DataFrame()

        if (len(adata_microglia[adata_microglia.obs[group_col] == clust_uniq[0]]) < 5) or (len(adata_microglia[adata_microglia.obs[group_col] == clust_uniq[1]]) < 5):
            all_groups_type_sheet = np.delete(all_groups_type_sheet, np.where(all_groups_type_sheet == cell_type_to_extract))
            idx +=1
            continue

        adata_microglia.obs[group_col] = adata_microglia.obs[group_col].astype(str)
        #sc.tl.dendrogram(adata, groupby = 'L04_newnum_subclassname')
        sc.tl.rank_genes_groups(adata_microglia, groupby=group_col, method="wilcoxon", tie_correct = True, pts = True,
                                #  layer = 'log1p_norm'
                                )
                
        for i in adata.obs[group_col].unique():
            dat1 = sc.get.rank_genes_groups_df(adata_microglia, group=i)
            dat1['group'] = i
            dat = pd.concat([dat, dat1])
            dat["mean_count"] = dat["names"].map(dict(zip(adata_microglia.var.index, adata_microglia.var.mean_counts)))

        if filters_bool:
            dat_filter = dat[ ### Choose filters here
            # (dat['pct_nz_group'] > filters_dic['percentage_pop']) & #Percentage of cell expressing the gene
            (dat['pvals_adj']<= filters_dict['pval_adj']) & # adjusted p-value
            (abs(dat['logfoldchanges']) > filters_dict['logfoldchanges']) & # logfoldchange
            (dat['mean_count'] >= filters_dict['mean_count']) &
            (dat['group'] != grp_ctrl)
            ]
            dfs_filter.append(dat_filter)

        idx +=1
        clear_output()
        bar.update(idx)

        dfs.append(dat)
    else:
        clear_output()
        print('Extraction done')

    bar.finish()

    if not os.path.exists(f"{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use}"):
        os.makedirs(f"{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use}")

    writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use}/DEG_{cluster_to_use}_no-filter.xlsx', engine='xlsxwriter')
    for j in range(0,len(dfs)):
        # print(j, " : ", all_cell_type_sheet[j])
        dfs[j].to_excel(writer, sheet_name=all_groups_type_sheet[j], index=False)
    writer.close()

    if filters_bool:
        writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use}/DEG_{cluster_to_use}_filter.xlsx', engine='xlsxwriter')
        for j in range(0,len(dfs_filter)):
            dfs_filter[j].to_excel(writer, sheet_name=all_groups_type_sheet[j], index=False)
        writer.close()

    return dfs

def DEG_two_conditions(adata: sc.AnnData,
                      name_dir: str,
                      cluster_to_use_1: str,
                      cluster_to_use_2: str,
                      group_col: str,
                      grp_ctrl: str,
                      filters_bool: bool,
                      filters_dict: dict,
                      dir_processed: str = dir_processed,
):
    dfs = []
    dfs_filter = []
    all_groups_C1 = np.array(adata.obs[cluster_to_use_1].unique())

    bar = progressbar.ProgressBar(maxval=len(all_groups_C1)+1, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()

    for idx, group_C1 in enumerate(all_groups_C1):
        
        clust_uniq = adata.obs[group_col].unique()
        
        adata3 = adata[adata.obs[cluster_to_use_1] == group_C1]

        
        if (len(adata3[adata3.obs[group_col] == clust_uniq[0]]) < 2) or (len(adata3[adata3.obs[group_col] == clust_uniq[1]]) < 2):
            continue

        dfs = []
        dfs_filter = []

        all_group_C2 = np.array(adata3.obs[cluster_to_use_2].unique())
        all_group_C2_sheet = all_group_C2
        
        
        for cell_type_to_extract in all_group_C2:
            print(f'Starting group: {cell_type_to_extract} in {group_C1}')

            adata_microglia = adata3[adata3.obs[cluster_to_use_2] == cell_type_to_extract]
            sc.pp.calculate_qc_metrics(adata_microglia, percent_top=[20, 50], inplace=True)

            ### Extract gene expression per cluster + log fold change + p-value
            
            dat = pd.DataFrame()

            if (len(adata_microglia[adata_microglia.obs[group_col] == clust_uniq[0]]) < 2) or (len(adata_microglia[adata_microglia.obs[group_col] == clust_uniq[1]]) < 2):
                all_group_C2_sheet = np.delete(all_group_C2_sheet, np.where(all_group_C2_sheet==cell_type_to_extract))
                clear_output()
                continue

            adata_microglia.obs[group_col] = adata_microglia.obs[group_col].astype(str)
            sc.tl.rank_genes_groups(adata_microglia, groupby=group_col, method="wilcoxon", tie_correct = True, pts = True)
            
            for i in adata3.obs[group_col].unique():
                # print(f"Cluster {cell_type_to_extract}_{i}")
                dat1 = sc.get.rank_genes_groups_df(adata_microglia, group=i)
                dat1['group'] = i
                dat = pd.concat([dat, dat1])
                dat["mean_count"] = dat["names"].map(dict(zip(adata_microglia.var.index, adata_microglia.var.mean_counts)))


            if filters_bool:
                dat_filter = dat[ ### Choose filters here
                    (dat['pct_nz_group'] > filters_dict['percentage_pop']) & #Percentage of cell expressing the gene
                    (dat['pvals_adj']<= filters_dict['pval_adj']) & # adjusted p-value
                    (abs(dat['logfoldchanges']) > filters_dict['logfoldchanges']) & # logfoldchange
                    (dat['mean_count'] >= filters_dict['mean_count']) & 
                    (dat['group'] != grp_ctrl)
                ]
                dfs_filter.append(dat_filter)

            dfs.append(dat)
            clear_output()
            
        if not os.path.exists(f"{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use_2}_in_{cluster_to_use_1}"):
            os.makedirs(f"{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use_2}_in_{cluster_to_use_1}")

        writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use_2}_in_{cluster_to_use_1}/{group_C1}_all_{cluster_to_use_2}_DEG.xlsx', engine='xlsxwriter')
        for j in range(0,len(dfs)):
            dfs[j].to_excel(writer, sheet_name=all_group_C2_sheet[j], index=False)
        writer.close()

        if filters_bool:
            writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges/{cluster_to_use_2}_in_{cluster_to_use_1}/{group_C1}_all_{cluster_to_use_2}_DEG_filter.xlsx', engine='xlsxwriter')
            for j in range(0,len(dfs_filter)):
                dfs_filter[j].to_excel(writer, sheet_name=all_group_C2_sheet[j], index=False)
            writer.close()

        clear_output()

        bar.update(idx)
        
    else:
        clear_output()
        print('Extraction done')

def deseq2_one_condition(adata:sc.AnnData,
                         name_dir:str,
                         cluster_to_use: str,
                         group_col: str,
                         filters_bool: bool,
                         filters_dict: dict,
                         pseudoreplicates: int = 3,
                         dir_processed: str = dir_processed,
                         ):

    list_celltypes = adata.obs[cluster_to_use].unique()
    list_celltypes_sheet = list_celltypes

    ddf = []
    ddf_filter = []
    list_ignored = []

    groups = adata.obs[group_col].unique()

    for idx, cell in enumerate(list_celltypes):  ### With replicates 
        print(cell, idx+1, "/", len(list_celltypes))
        adata_sub = adata[adata.obs[cluster_to_use]==cell]
        if len(adata_sub.obs[group_col].unique()) == 1:
            list_celltypes_sheet = np.delete(list_celltypes_sheet, np.where(list_celltypes_sheet == cell))
            list_ignored.append(cell)
            continue
        
        pbs = [] 
        for sample in adata_sub.obs['sample'].unique():
            print(sample)
            samp_adata_sub = adata_sub[adata_sub.obs['sample']==sample]
            


            samp_adata_sub.X = samp_adata_sub.layers['counts']

            random.seed(20150201)
            indices = list(samp_adata_sub.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), pseudoreplicates)

            for i, pseudo_rep in enumerate(indices):

                rep_adata = sc.AnnData(X = samp_adata_sub[indices[i]].X.sum(axis=0),
                                    var = samp_adata_sub[indices[i]].var[[]])
                
                rep_adata.obs_names = [sample + '_' + str(i)]
                rep_adata.obs[group_col] = samp_adata_sub.obs[group_col].iloc[0]
                rep_adata.obs['replicate'] = i

                pbs.append(rep_adata)

        pb = sc.concat(pbs)
        counts = pd.DataFrame(pb.X, columns = pb.var_names)
        dds = DeseqDataSet(counts = counts,
                    metadata = pb.obs,
                    design_factors = [group_col],
                    quiet = True)
        sc.pp.filter_genes(dds, min_cells = 1)
        dds.deseq2()
        stat_res = DeseqStats(dds, n_cpus = 32, contrast = (group_col,groups[0],groups[1]))

        stat_res.summary() 
        de = stat_res.results_df
        
        if filters_bool:
            de_filter = de[(abs(de['log2FoldChange'])>filters_dict['logfoldchanges']) &
                        (de['padj']<filters_dict['padj'])
                        ]
            ddf_filter.append(de_filter)
        
        ddf.append(de)
        clear_output()
        

    if not os.path.exists(f"{dir_processed}/analysis/{name_dir}/foldchanges_DeSeq2/{cluster_to_use}"):
        os.makedirs(f"{dir_processed}/analysis/{name_dir}/foldchanges_DeSeq2/{cluster_to_use}")

    writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges_DeSeq2/{cluster_to_use}/DEG_DeSeq2_celltype_no-filter.xlsx', engine='xlsxwriter')
    for j in range(len(list_celltypes_sheet)):
        ddf[j].to_excel(writer, sheet_name=list_celltypes_sheet[j], index=True)
    writer.close()

    if filters_bool:
        writer = pd.ExcelWriter(f'{dir_processed}/analysis/{name_dir}/foldchanges_DeSeq2/{cluster_to_use}/DEG_DeSeq2_celltype_filter.xlsx', engine='xlsxwriter')
        for j in range(len(list_celltypes_sheet)):
            ddf_filter[j].to_excel(writer, sheet_name=list_celltypes_sheet[j], index=True)
        writer.close()
    
    clear_output()

    print('Analysis Done')
    print(f'{len(ddf)} analyzed.')
    print(f'{len(list_ignored)} ignored.')

    return ddf, ddf_filter, list_celltypes_sheet

def interactive_volcano_plot(result_list:list,
                             deg_method:str,
                             key,
                             ctrl_grp:str,
                             test_grp:str,
                             pval_thshld:float = 0.05,
                             FC_thshld:float = 0.26,
                             ):
    
    '''
    Plot interractive volcano plot for quick exploration of the result. Cannot be saved at the moment.
    result_list = output from DEG analysis. Only work for one condition DEG for now
    deg_method = method of deg used. either wilcoxon or DeSeq2 for now
    key = element you want to plot (specific cell type for instance)
    ctrl_grp and test_grp = name of each group to avoid double plotting
    pval_thrshld = threshold for P-VALUE significance (default = 0.05)
    FC_thrshld = Threshold for logfoldchanges significance (default = 0.26 which ~= 1.2 FC)
    '''
    
    if deg_method == 'wilcoxon':
        logstr = 'logfoldchanges'
        pval_str = 'pvals'
        adjpval_str = 'pvals_adj'
        gene_name = 'names'

        for idx in range(len(result_list)):
            if ctrl_grp in result_list[idx]["group"].unique():
                result_list[idx] = result_list[idx][result_list[idx]['group'] == test_grp]

    elif deg_method == 'DeSeq2':
        logstr = 'log2FoldChange'
        pval_str = 'pvalue'
        adjpval_str = 'padj'
        gene_name = 'feature_name'
    else:
        print('Wrong method input. Either wilcoxon or DeSeq2.')
        

    # Prepare data for visualization
    min_thr = result_list[key][result_list[key][adjpval_str] != 0][adjpval_str].min()
    result_list[key]['neg_log10_p'] = -np.log10(result_list[key][pval_str] + min_thr)
    result_list[key]['neg_log10_padj'] = -np.log10(result_list[key][adjpval_str]+ min_thr)

    # Create significance categories based on both thresholds
    significance_threshold = -np.log10(pval_thshld)
    fold_change_threshold = FC_thshld

    result_list[key]['significance'] = 'Not-Significant'
    result_list[key].loc[
        (result_list[key]['neg_log10_padj'] > significance_threshold) & 
        (abs(result_list[key][logstr]) > fold_change_threshold), 
        'significance'
    ] = 'Significant'

    volcano_plot = result_list[key].hvplot.scatter(
        x=logstr, 
        y="neg_log10_padj",
        c='significance',
        cmap={'Not-Significant': 'lightgrey', 'Significant': 'black'},
        hover_cols=[gene_name, 'significance'],
        title=f"DEG {ctrl_grp} vs {test_grp} ({deg_method}): {key}",
        legend='bottom_right',
        alpha=0.6,
        size=20,
        responsive=True,
        height=500
    )

    # Add threshold lines
    (
        volcano_plot
        * hv.HLine(significance_threshold).opts(color='red', line_dash='dashed', line_width=2)
        * hv.VLine(-fold_change_threshold).opts(color='blue', line_dash='dashed', line_width=2) 
        * hv.VLine(fold_change_threshold).opts(color='blue', line_dash='dashed', line_width=2)
    )

    return volcano_plot