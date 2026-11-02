from module.config_local import dir_processed
import os
import argparse

parser = argparse.ArgumentParser(
    description='Create folders for processed files'
    )

parser.add_argument('--name_dir', type=str, help='Name of the experience')
args = parser.parse_args()

name_dir = args.name_dir

def create_folders(name_dir:str,
                   dir_processed: str = dir_processed
                   ):
    if not os.path.exists(f"{dir_processed}/{name_dir}/csv/"):
        os.makedirs(f"{dir_processed}/{name_dir}/csv/")
        print("csv folder created")
    if not os.path.exists(f"{dir_processed}/{name_dir}/QC/"):
        os.makedirs(f"{dir_processed}/{name_dir}/QC/")
        print("QC folder created")
    if not os.path.exists(f"{dir_processed}/{name_dir}/h5ad/"):
        os.makedirs(f"{dir_processed}/{name_dir}/h5ad/")
        print("h5ad folder created")
    if not os.path.exists(f"{dir_processed}/{name_dir}/analysis/"):
        os.makedirs(f"{dir_processed}/{name_dir}/analysis/")
        print('Analysis folder created')
    if not os.path.exists(f"{dir_processed}/{name_dir}/plot/"):
        os.makedirs(f"{dir_processed}/{name_dir}/plot/")
        print('Plot folder created')
    if not os.path.exists(f"{dir_processed}/{name_dir}/Correlation_Mapping/"):
        os.makedirs(f"{dir_processed}/{name_dir}/Correlation_Mapping/")
        print('Correlation_Mapping folder created')
    with open(f"{dir_processed}/{name_dir}/folders_created.txt", "w") as f:
        f.write("folders created")

create_folders(name_dir = name_dir)