name_dir = 'circa-only'
samples = ['circa4-IGM-ZT01','circa4-IGM-ZT05','circa4-IGM-ZT09','circa4-IGM-ZT13','circa4-IGM-ZT17','circa4-IGM-ZT21','circa4-Salk-ZT01','circa4-Salk-ZT05','circa4-Salk-ZT09','circa4-Salk-ZT13','circa4-Salk-ZT17','circa4-Salk-ZT21']

rule all:
	input:
		expand("../data_processed/{name_dir}/folders_created.txt", name_dir = name_dir),
		expand("../data_processed/{name_dir}/QC/undernoise_{sample}.csv", name_dir=name_dir, sample=samples),
		expand("../data_processed/{name_dir}/polygons/{sample}_cells.geojson", name_dir=name_dir, sample=samples),
		expand("../data_raw/{name_dir}/{sample}/transcripts_light.parquet", name_dir=name_dir, sample=samples),
		expand("../data_processed/{name_dir}/h5ad/{sample}_forMMC.h5ad", name_dir=name_dir, sample=samples),
		expand("../data_processed/{name_dir}/Correlation_Mapping/{sample}_hierarchical_mapping.csv", name_dir=name_dir, sample=samples)



rule QC:
	input:
		"../data_raw/{name_dir}/{sample}/transcripts_light.parquet"
	output:
		"../data_processed/{name_dir}/QC/undernoise_{sample}.csv"
	params:
		name_dir = name_dir,
		sample_name = "{sample}",
		env = "jupnot_snake"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_QC.py --name_dir {params.name_dir} --sample {params.sample_name}"

rule MMC:
	output:
		"../data_processed/{name_dir}/h5ad/{sample}_forMMC.h5ad"
	params:
		name_dir = name_dir,
		sample_name = "{sample}",
		env = "jupnot_snake"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_MMC.py --name_dir {params.name_dir} --sample {params.sample_name}"

rule polygons:
	input:
		"../data_raw/{name_dir}/{sample}/cell_boundaries.parquet"
	output:
		"../data_processed/{name_dir}/polygons/{sample}_cells.geojson"
	params:
		sample_name = "{sample}",
		name_dir = name_dir,
		env = "jupnot_snake"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_polygons.py --sample {params.sample_name} --name_dir {params.name_dir}"

rule transcript_lighter:
	input:
		"../data_raw/{name_dir}/{sample}/transcripts.parquet"
	output:
		"../data_raw/{name_dir}/{sample}/transcripts_light.parquet"
	params:
		sample_name = "{sample}",
		name_dir = name_dir,
		env = "jupnot_snake"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_lighter_transcripts_files.py --sample {params.sample_name} --name_dir {params.name_dir}"

rule folders:
	priority:
		100
	output:
		"../data_processed/{name_dir}/folders_created.txt"
	params:
		name_dir = name_dir,
		env = "jupnot_snake"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_folders.py --name_dir {params.name_dir}"

rule mapping:
	input:
		"../data_processed/{name_dir}/h5ad/{sample}_forMMC.h5ad"
	output:
		"../data_processed/{name_dir}/Correlation_Mapping/{sample}_hierarchical_mapping.csv"
	params:
		name_dir = name_dir,
		sample_name = "{sample}",
		env = "cell_mapper_cpython"
	shell:
		"conda run --no-capture-output -n {params.env} python module/script_mapping.py --name_dir {params.name_dir} --sample {params.sample_name}"