import pathlib
import os
from config_local import dir_processed
import argparse

parser = argparse.ArgumentParser(description='Initial QC investigations')

parser.add_argument('--name_dir', type=str, help='Name of the experience')
parser.add_argument('--sample', type=str, required=True, help='Sample name to process')

args = parser.parse_args()
name_dir = args.name_dir
sample_id = args.sample

os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

from cell_type_mapper.cli.from_specified_markers import (
    FromSpecifiedMarkersRunner
)

ref_dir = pathlib.Path('reference_files')
processed_dir = pathlib.Path(dir_processed)

baseline_marker_path = ref_dir / 'mouse_markers_230821.json'
baseline_precompute_path = ref_dir / 'precomputed_stats_ABC_revision_230821.h5'
test_h5ad_path = processed_dir / f'{name_dir}/h5ad/{sample_id}_forMMC.h5ad'

## Define output files
baseline_json_output_path = processed_dir / f'{name_dir}/Correlation_Mapping/{sample_id}_hierarchical_mapping.json'
baseline_csv_output_path = processed_dir / f'{name_dir}/Correlation_Mapping/{sample_id}_hierarchical_mapping.csv'

baseline_mapping_config = {
    'query_path': str(test_h5ad_path),
    'extended_result_path': str(baseline_json_output_path),
    'csv_result_path': str(baseline_csv_output_path),
    'tmp_dir': str(processed_dir),
    'max_gb': 10,
    'cloud_safe': False,
    'verbose_stdout': False,
    "map_to_ensembl": True,
    'type_assignment': {
        'normalization': 'raw',
        'n_processors': 4,
        'chunk_size': 10000,
        'bootstrap_iteration': 100,
        'bootstrap_factor': 0.5,
        'rng_seed': 511020
    },
    'precomputed_stats': {
        'path': str(baseline_precompute_path)
    },
    'query_markers': {
        'serialized_lookup': str(baseline_marker_path)
    },
    

}

mapping_runner = FromSpecifiedMarkersRunner(
    args=[],
    input_data=baseline_mapping_config
)

mapping_runner.run()

os.delete(baseline_json_output_path)