v. **6.0** (10/24/24)

*New feature* : Optional batch correction during Banksy (Harmony)

*Output*:
* Fold-change, pvalue can be automatically extracted for all cell types

*UMAP*:
- New UMAP to plot gene expression accross UMAP
- Only use UMAP from Banksy

*Sub-clustering*:
* improved workflow by automatically using cluster number when creating the dictionary to rename cell types
* added visualization of UMAP and PCA
* added option to choose between leiden and kmeans algorithms
* added visualization of markers genes for subclusters

*Analysis*:
* reorganized the analysis section to separate celltype, regions and both
* added plotting of violin plots, stacked violin plots
* Dotplot for markers genes


v. **5.31** (09/24/24)

*New feature*: Now supporting multi-sample Banksy
* Now all samples from a same experiment can be process at the same time and automatically with Banksy. The loop combine adatas and normalized/weighted expression into a dictionnary then used for dimension reduction and clustering.
* Plot step will raise an error but is still essential to get the 'cell type' columns in results_df
* Take really long time. 1.5 h for 3 samples; 3h for 6 samples.
* Dimension reduction saturate CPU (desktop) and clustering use around 30Go of RAM

*New feature*: importation and alignment of regions from Qupath
* tested with A-beta plaques, automatically detected in Qupath
* require key points if alignment is needed (other images than morphology.OME.tiff)
* Could be use for region registering if we improve the process with ABBA

*bug fix*
* to adapt to geoPandas 1.0 way to add a name to index, we drop said index now before mapping the cells and after Banksy clustering.
* Known issue, not reproduced: impossible to read h5ad file (something about obs not existing, I should have saved the error message...)
	* fix: added a section in "output file" to merge previous h5ad with cell type column from csv backup

*output*
* To reduce output file size, most h5ad and csv files are now compressed using gzip, it reduces the size of the output between $1/10^{th}$ and $1/2$ of original size. Good for now but might not be enough for 5K panel. 
	* /!\ You will have to read and resave your files to have them compressed. It is tedious but save space in the long term
* Gzip files can be extract like normal zip files but can also be directly read in python/R:

```{R}
read_csv("/path/filename.csv.gz")
```


v. **5.2** (2024-09-04)

*New feature*: Automap (temporary name)
*  <span style = "color:red">Require Module file 'automap.py' in same folder as the notebook </span>
* Automapping of regions based on cluster annotation
* Simplify clusters by merging similar region (CTX, HIPP, etc.; can be modified in the dictionary in Data pre-processing)
* Do not include glial and epithelial cells except VLMC and Ependymal cells (can be modified in the dictionary in Data pre-processing)
* Make sure to run all your samples before continuing to 'Match cells'

*New feature*: Match cells in automatically drawn regions
* Depend on Automap-generated geojson files
* Automatic region annotation based on cluster name from Automap
* Can probably be improved, some steps might be useless now

*misc.
* added timer for some long functions, help estimating time when working on multiple samples in a row
* added "more visible" warning when going back to process other samples/cluster/etc. is needed. Hope it's not too aggressive...

---- 
v. **5.1** (2024-08-30)

*output
- Removed saving of …norm_MMC.csv (redundant with other csv files)
- Changed save file format for scatter (UMAP and sections) from svg to png (svg where impossible to open)

*visualization
- Changed annotations in UMAP from cluster number to cell type name
- UMAP for individual samples is now Banksy UMAP (previously, the position of clusters annotation was unreliable)
- UMAP grouping all samples from the post-Banksy UMAP
- Added legend with cell type names
- Changed color palette for Banksy plotting (purely aesthetical reasons)
- Fixed visualizations for ‘total_transcript’, ‘number of genes’ and ‘MMC correlation coefficient’
- need to change the ‘to_use’ parameter (for now, I’ll find a fix)

*misc.
- Moved Output files after analysis to list all possible output there

  

- To do:
	- [ ]  find a way to save editable scatter plot (downsampling of clusters?)
	- [ ]  fix ‘to_use’ in data visualization
	- [ ]  order in legend (alphabetical, class?)
	- [x] Streamline Mapping from drawn regions
		- [x] resegmented runs
		- [x] other runs
	- [x] Make functions for plot of UMAP and sections
	- [x] Automatically relevant cells information and normalized gene count
		- [x] These CSV are already of 1Go for 6 samples / 247 genes, we need a better solution to be ready to work with the 5k panel
			- [x] Zarr? Gzip?
	- [ ] only keep useful columns to reduce size of files