### Config file
# modify to fit your need then rename the file as "config_local.py". Do not move it.

## Path to the raw Xenium files
dir_raw = ""

## Path to where you want to save your processed files
dir_processed = ''

## Path containing the repository
dir_main = ""

## Name of your experiment
# name_dir = ""

# order of cell types for plotting
# order-list = 


def sample_name_import(name_experiment:str):
    
    dict_exp_name = {
        # 'circa-SD' : ['circa4-IGM-ZT01','circa4-IGM-ZT05','circa4-IGM-ZT09','circa4-IGM-ZT13','circa4-IGM-ZT17','circa4-IGM-ZT21', "SD1-ZT01","SD1-ZT05","SD1-ZT09","SD1-ZT13","SD1-ZT17","SD1-ZT21"], #EXAMPLE

        }
    
    return dict_exp_name[name_experiment]

#Default colors?