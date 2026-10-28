import argparse
import os

from lighter_transcript_parquet import lighter_transcript_parquet_indi
from config_local import dir_raw

parser = argparse.ArgumentParser(
    description='Import Xenium raw files and exctract gene expression for MMC.'
    )

parser.add_argument('--sample', type=str, required=True, help='Sample name to process')
parser.add_argument('--name_dir', type=str, help='Name of the experience')
args = parser.parse_args()
sample_id = args.sample
name_dir = args.name_dir

lighter_transcript_parquet_indi(sample = sample_id,
                                gene_list = None,
                                dir_raw = dir_raw,
                                name_dir = name_dir
                                )