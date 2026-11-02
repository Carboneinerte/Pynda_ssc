from module.config_local import dir_raw, dir_processed
from module.quality_control import undernoise_list_indi
import argparse

parser = argparse.ArgumentParser(description='Initial QC investigations')

parser.add_argument('--name_dir', type=str, help='Name of the experience')
parser.add_argument('--sample', type=str, required=True, help='Sample name to process')

args = parser.parse_args()
name_dir = args.name_dir
sample_id = args.sample

undernoise_list_indi(sample = sample_id,
                     name_dir=name_dir,
                     dir_raw = dir_raw,
                     dir_processed = dir_processed
                    )