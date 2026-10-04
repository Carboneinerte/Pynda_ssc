# rule preprocess:
# 	input:
# 		"../data_raw/data/ret_annot_subsample.h5ad.gz"
# 	output:
# 		"../data_processed/ret_data/test_automate_snake.h5ad.gz"
# 	shell:
# 		"python test_autom.py"

name_dir = "AD_10x"
samples = ["TgCRND8_2_5_months", "TgCRND8_5_7_months","wildtype_2_5_months", "wildtype_5_7_months"]

rule all:
    input:
        #expand("../data_processed/{name_dir}/polygons/{sample}_cells.geojson", name_dir=name_dir, sample=samples)
		expand("../data_processed/{name_dir}/h5ad/{name_dir}_{sample}_forMMC.h5ad", name_dir=name_dir, sample=samples)
	
rule MMC:
	input:
		"../data_raw/{name_dir}/{sample}/cell_feature_matrix.h5"
	output:
		"../data_processed/{name_dir}/h5ad/{name_dir}_{sample}_forMMC.h5ad"
	params:
		sample_name = samples,
		name_dir = name_dir
	shell:
		"python module/MMC.py --samples {params.sample_name} --name_dir {params.name_dir}"

rule polygons:
	input:
		"../data_raw/{name_dir}/{sample}/cell_boundaries.csv.gz"
	output:
		"../data_processed/{name_dir}/polygons/{sample}_cells.geojson"
	params:
		sample_name = "{sample}",
		name_dir = name_dir
	shell:
		"python module/polygons.py --sample {params.sample_name} --name_dir {params.name_dir}"