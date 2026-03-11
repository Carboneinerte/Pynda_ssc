## To do list before starting the analysis:
- Clone or copy repository from Github (main branch)
- Install Conda and create environment
```{python}
conda create -n "env_name" python==3.11

conda activate "env_name"

pip install -r requirements.txt
```
- Copy module/config.py and rename it as module/config_local.py
    - Define the folder containing raw files (**dir_raw**)
    - Define the folder containing processed files (**dir_notebook**)
    
- Define the name of the experiment (**name_dir**)
- Define the name of each samples
    - Rename folders in "name_dir" accordingly
    - Add the names as a list in the function "sample_name_import" in module/misc.py
```{python}
"name_dir" : ["sample1", "sample2", etc.],
```

- Draw region of interest in Xenium explorer and export the list of cells
    - save as : "{sample}_Whole-section_cells_stats.csv"
!["Xenium explorer export cell list"](module\image_readme\XE_cell_list.png)

### Optional
- Create "light" transcripts.parquet with only high quality gene transcripts (useful for unassigned transcripts analysis)
- Run [this notebook](./v11_Polygons_preprocessing.ipynb) if you plan to do plots with cell polygons
- Run [this notebook](./V11Z_Xenium_Quality_control.ipynb) to plot raw metrics and compare samples


## Folder organization

/root/\
&nbsp;|__data &nbsp;&nbsp;(="dir_raw")\
&nbsp;|__notebook &nbsp;&nbsp;(="dir_notebook")\
&nbsp;&nbsp;&nbsp;|__analysis\
&nbsp;&nbsp;&nbsp;|__coordinates\
&nbsp;&nbsp;&nbsp;|__csv &nbsp;&nbsp;(Will also contain parquet files)\
&nbsp;&nbsp;&nbsp;|__h5ad\
&nbsp;&nbsp;&nbsp;|__module\ (contains the .py files with the functions used in notebooks)
&nbsp;&nbsp;&nbsp;|__plot

## List of minimum files to have in data directory (for each sample)

- cell_boundaries.parquet
- cell_feature_matrix.h5
- cells.csv.gz
- metrics_summary.csv
- transcripts.parquet

## Gene lists

In module/misc, you can define your own list of genes to use in different functions (heatmap, violin plot, etc.)


<!--
## reload module

```{python}
import importlib
importlib.reload(xp)
```
-->