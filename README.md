## To do list before starting the analysis:

-   Clone or copy repository from Github (main branch) where you want to store it. It will be your "dir_main" folder.
-   Config git

``` python
git config user.name 'yourusername'
git config user.email rbiopandas2@gmail.com # or your own Github account
```

-   Install Conda and create environment

    -   On Jetstream, avoid installing the env in /home/, prefer to install it in /media/volume/...

``` python
conda create --prefix=/media/volume/.../yourEnvName python==3.11

conda activate yourEnvName

pip install -r requirements.txt

pip install leidenalg, xlsxwriter, geojson
```
- Copy module/config.py and rename it as module/config_local.py
    - Define the folder containing raw files (**dir_raw**)
    - Define the folder containing processed files (**dir_processed**)
    
- Define the name of the experiment (**name_dir**)
- Define the name of each samples
    - Rename folders in "name_dir" accordingly
    - Add the names as a list in the function "sample_name_import" in module/config_local.py
```{python}
"name_dir" : ["sample1", "sample2", etc.],
```

-   To start using the notebooks, **copy** one from "notebooks_blank" and rename it. It must ends with "...\_in_use.ipynb".

-   Please, do not modufy or use the notebooks in "notebooks_blank". You won't be able to import function from the module anyway.

### Optional

-   Draw region of interest in Xenium explorer and export the list of cells
    -   save as : "{dir_raw}/{sample_name}/{sample_name}\_ROI_cells_stats.csv"\
        !["Xenium explorer export cell list"](module\image_readme\XE_cell_list.png)
-   Create "light" transcripts.parquet with only high quality gene transcripts (useful for unassigned transcripts analysis)
-   Run [this notebook](./v11_Polygons_preprocessing.ipynb) if you plan to do plots with cell polygons
-   Run [this notebook](./V11Z_Xenium_Quality_control.ipynb) to plot raw metrics and compare samples

## Folder organization

/root/\
 \|-- main (contains git repo)   (="dir_main")\
   \|-- reference_files\
   \|-- module (contains the .py files with the functions used in notebooks)\
   \|-- notebooks_blank (contains blank notebook you will copy into your root folder to use)\
 \|-- raw_data   (="dir_raw")\
   \|-- "Sample1"\
   \|-- "Sample2"\
   \|-- "..."\
 \|-- processed_data   (="dir_processed")\
   \|-- analysis\
   \|-- coordinates\
   \|-- Coorelation_Mapping\
   \|-- csv   (Will also contain parquet files)\
   \|-- h5ad\
   \|-- plot

## List of minimum files to have in **data** directory (for each sample)

-   cell_boundaries.parquet
-   cell_feature_matrix.h5
-   cells.csv.gz
-   metrics_summary.csv
-   transcripts.parquet
-   {sample_name}\_ROI_cells_stats.csv (optional)

## Misc. info

-   In module/misc, you can define your own list of genes to use in different functions (heatmap, violin plot, etc.)
-   In module/misc, you can add your own annotation to export to a gene matrix (df)

## Recommended VSCode extensions

-   Data Wrangler
-   Rainbow Csv

## Useful links

-   <https://github.com/seandavi/awesome-single-cell?tab=readme-ov-file#experimental-design>
-   <https://www.w3schools.com/python/>
-   <https://www.sc-best-practices.org/preamble.html>
-   <https://scanpy.readthedocs.io/en/stable/>
-   <https://squidpy.readthedocs.io/en/stable/>

````{=html}
<!--
## reload module

```{python}
import importlib
importlib.reload(xp)
```
-->
````
