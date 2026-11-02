# Presentation
The present repository contains a pipeline designed for analysis of single-cell spatial transcriptomics dataset. It was developt in the "Regulatory Biology Laboratory" at Salk Institute for Biological Studies (San Diego, USA).

The workflow is designed to analyze Xenium dataset from scratch but other SC/SSC datasets can be imported and processed after the preprocessing part.

It is not as of today available as a package but feel free to follow to instructions below to install and use it.


## To do list before starting the analysis:

-   Clone or copy repository from Github (main branch) where you want to store it. It will be your "dir_main" folder.
``` python
git clone https://github.com/Carboneinerte/Pynda_ssc.git
```
-   Config git

``` python
git config user.name 'yourusername'
git config user.email rbiopandas2@gmail.com # or your own Github account
git init
```

-   Install Conda and create environment

    -   On Jetstream, avoid installing the env in /home/, prefer to install it in /media/volume/...
    -   Make sure you have still some space in /home/, conda apparently need some space (for cache ?)

``` python
conda create --prefix=/media/volume/.../yourEnvName python==3.11

conda activate yourEnvName

pip install -r requirements.txt

pip install leidenalg, xlsxwriter, geojson, goatools, polars 
```
(use: pip freeze > requirements.txt when changing packages version in the env. Here, I saved you yet another google search.)

- Copy module/config.py and rename it as module/config_local.py
    - Define the folder containing raw files (**dir_raw**)
    - Define the folder containing processed files (**dir_processed**)
    
- Define the name of the experiment (**name_dir**)
- Define the name of each samples
    - Rename folders in "name_dir" accordingly
    - Add the names as a list in the function "sample_name_import" in module/config_local.py
    - I would recommand to include an indication of the experimental group in the name of each samples (not required but can make things smoother downhill)
```{python}
"name_dir" : ["sample1", "sample2", etc.],
```

-   To start using the notebooks, **copy** one from "notebooks_blank" and rename it. It must ends with "...\_in_use.ipynb".

-   Please, do not modify or use the notebooks in "notebooks_blank". You won't be able to import function from the module anyway.

### Optional

-   Draw region of interest in Xenium explorer and export the list of cells
    -   save as : "{dir_raw}/{sample_name}/{sample_name}\_ROI_cells_stats.csv"\
        !["Xenium explorer export cell list"](image_readme\XE_cell_list.png)
-   Create "light" transcripts.parquet with only high quality gene transcripts (useful for unassigned transcripts analysis)
-   Run [this notebook](./v11_Polygons_preprocessing.ipynb) if you plan to do plots with cell polygons
-   Run [this notebook](./V11Z_Xenium_Quality_control.ipynb) to plot raw metrics and compare samples

## Project Structure
```
root/
│
├── Pynda_ssc/          # contains git repo (="dir_main")
│   ├── reference_files
│   ├── module          # contains the .py files with the functions used in notebooks
│   └── notebooks_blank  # contains blank notebooks you will copy into your root folder to use
├── raw_data            # contains raw files (="dir_raw")
│   ├──  "Sample1"
│   ├──  "Sample2"
│   └──  "..."
├──  processed_data     # contains project outputs (="dir_processed")
│   ├──  analysis
│   ├──  coordinates
│   ├──  Coorelation_Mapping
│   ├──  csv            # Will also contain parquet files
│   ├──  h5ad
│   └──  plot
├── .gitignore          # Git ignore file
├── README.md           # Project overview and instructions
└── requirements.txt    # Python package dependencies
```

## List of minimum files to have in **data** directory (Xenium ony, for each sample)

-   cell_boundaries.parquet
-   cell_feature_matrix.h5
-   cells.csv.gz
-   metrics_summary.csv
-   transcripts.parquet
-   {sample_name}\_whole_section_annotation.csv (optional, contain the coordinates of the ROI)
-   {sample_name}\_ROI_cells_stats.csv (optional, contain the list of cells in the ROI)

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


<!--
## reload module

```{python}
import importlib
importlib.reload(xp)
```
-->

