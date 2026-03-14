import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from module.config_local import dir_processed, dir_raw, dir_main
from module.misc import save_figure
from datetime import datetime
today = datetime.today().strftime('%Y-%m-%d')

def plot_mcc_density (sample_1: list,
                      name_dir: str,
                      save_plot: bool = False,
                      save_name : str = 'CC_mmc.png'
                      ):
    
    for sample in sample_1:
        print(sample)
        df = pd.read_csv(f'{dir_processed}/Correlation_Mapping/{name_dir}/{name_dir}_{sample}_CorrelationMapping.csv', comment = "#")
        df['sample'] = df['cell_id'].map(lambda name: name.split('_')[0])
        df_temp = df.filter(['sample', 'subclass_correlation_coefficient'])
        if "df_all" not in locals():
            df_all = df_temp
        else:
            df_all = pd.concat([df_all, df_temp])

    dict_ = dict(zip(df_all['sample'],df_all['sample']))
    dict_.update({"Region1": "february-test","S1" : "march-test"})
    df_all['sample'] = df_all['sample'].map(dict_)

    order_sample = df_all['sample'].unique() ## All samples

    ax = sns.histplot(data = df_all, x = 'subclass_correlation_coefficient', hue='sample', hue_order= order_sample,
             element="step", cumulative= True, fill= False, common_norm=False,
             stat='density')
    plt.xlim(0,0.85)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    if save_name == 'CC_mmc.png':
        save_name = f'CC_{sample_1[0]}.png'

    if save_plot:
        plt.savefig(f'{dir_processed}/plot/{name_dir}/{save_name}')

def desc_metrics(samples_ids:list,
                 name_dir:str,
                 path_to_data:str = dir_raw,
                 reference:bool = False,
                 plot_average:bool = True,
                 save_plot: bool = False
                 ):
    
    parameters_to_plot = ['region_area', 'total_high_quality_decoded_transcripts','fraction_transcripts_decoded_q20', 'decoded_transcripts_per_100um2','estimated_number_of_false_positive_transcripts_per_cell',
                      'num_cells_detected', 'fraction_transcripts_assigned', 'median_genes_per_cell', 'median_transcripts_per_cell' ]
    

    for sample in samples_ids:
        with open(f"{path_to_data}\{sample}\metrics_summary.csv", 'r', encoding='utf-8') as file:
            file_content = pd.read_csv(file)
            
            if 'files_content' in locals():
                files_content = pd.concat([files_content,file_content])
            else:
                files_content = pd.DataFrame(file_content)

    parameters_to_plot = ['region_area', 'total_high_quality_decoded_transcripts','fraction_transcripts_decoded_q20',
                        'decoded_transcripts_per_100um2','estimated_number_of_false_positive_transcripts_per_cell','num_cells_detected',
                        'fraction_transcripts_assigned', 'median_genes_per_cell', 'median_transcripts_per_cell' ]

    if reference:
        try:
            reference_dataset = pd.read_csv(f'{dir_main}/reference_files/reference_dataset.csv')
            min_1 = 1
        except:
            print('No reference dataset found')
    else:
        min_1 = 0
    max_1 = min_1 + len(files_content) - 1 


    fig, axes = plt.subplots(3,3, figsize=(15,15))
    axes = axes.flatten()

    for n, ax in enumerate(axes):
        parameter = parameters_to_plot[n]
        if reference:
            ax.bar(x = reference_dataset['region_name'], height = reference_dataset[parameters_to_plot[n]], color = (0.32,0.13,0.102))
        ax.bar(x = files_content['region_name'], height = files_content[parameters_to_plot[n]], color = (0.898,0.603,0.32))
        if plot_average:
            ax.hlines(y = files_content[parameters_to_plot[n]].mean(), xmin = min_1, xmax = max_1, linestyles= 'dashed', colors = 'black')
        ax.set_title(parameter)
        ax.tick_params(axis = 'x', rotation = 90, direction = 'in', pad = -80)

    if save_plot:
        fig.savefig(f'{dir_processed}/plot/{name_dir}/{name_dir}_QC.svg')

def desc_metrics_double(samples_ids_1: list,
                        samples_ids_2: list,
                        name_dir: str,
                        dir_raw: str = dir_raw,
                        reference: bool = False,
                        plot_average: bool = True,
                        save_plot: bool = False
                        ):

    parameters_to_plot = ['region_area', 'total_high_quality_decoded_transcripts','fraction_transcripts_decoded_q20',
                            'decoded_transcripts_per_100um2','estimated_number_of_false_positive_transcripts_per_cell','num_cells_detected',
                            'fraction_transcripts_assigned', 'median_genes_per_cell', 'median_transcripts_per_cell' ]
    
    for sample in samples_ids_1:
        with open(f"{dir_raw}/{sample}/metrics_summary.csv", 'r', encoding='utf-8') as file:
            file_content_1 = pd.read_csv(file)
            
            if 'files_content_1' in locals():
                files_content_1 = pd.concat([files_content_1,file_content_1])
            else:
                files_content_1 = pd.DataFrame(file_content_1)

    for sample in samples_ids_2:
        with open(f"{dir_raw}/{sample}/metrics_summary.csv", 'r', encoding='utf-8') as file:
            file_content_2 = pd.read_csv(file)
            
            if 'files_content_2' in locals():
                files_content_2 = pd.concat([files_content_2,file_content_2])
            else:
                files_content_2 = pd.DataFrame(file_content_2)

    if reference:
        try:
            reference_dataset = pd.read_csv(f'{dir_main}/reference_files/reference_dataset.csv')
            min_1 = 1
        except:
            print('No reference dataset found')
    else:
        min_1 = 0
    max_1 = min_1 + len(files_content_1) - 1 
    min_2 = max_1 + 1
    max_2 = min_2 + len(files_content_2) - 1

    fig, axes = plt.subplots(3,3, figsize=(15,15))
    axes = axes.flatten()

    for n, ax in enumerate(axes):
        parameter = parameters_to_plot[n]
        if reference:
            ax.bar(x = reference_dataset['region_name'], height = reference_dataset[parameters_to_plot[n]], color = (0.32,0.13,0.102))
        ax.bar(x = files_content_1['region_name'], height = files_content_1[parameters_to_plot[n]], color = (0.898,0.603,0.32))
        ax.bar(x = files_content_2['region_name'], height = files_content_2[parameters_to_plot[n]], color = "lightblue")
        if plot_average:
            ax.hlines(y = files_content_1[parameters_to_plot[n]].mean(), xmin = min_1, xmax = max_1, linestyles= 'dashed', colors = 'black')
            ax.hlines(y = files_content_2[parameters_to_plot[n]].mean(), xmin = min_2, xmax = max_2, linestyles= 'dashed', colors = 'black')
        ax.set_title(parameter)
        ax.tick_params(axis = 'x', rotation = 90, direction = 'in', pad = -60)
    if save_plot:
        fig.savefig(f'{dir_processed}/plot/{name_dir}/{name_dir}_QC.svg')

def add_jitter(x, scale=0.5):
    return x + np.random.uniform(-scale, scale, size=len(x))

def noise_threshold_ploting(sample:str,
                            save_plot:bool = False):
    print('Start Sample :', sample)

    df = pd.read_parquet(f'D:/Xenium/{sample}/transcripts.parquet',
                            filters=[("qv",">=",20)]
                            )       

    data = pd.DataFrame({'feature_name': df.feature_name.value_counts().index,'count' : df.feature_name.value_counts()})
    data.sort_index(inplace=True)

    data['type'] = data['feature_name'].apply(lambda x: x.split('_')[0])

    percentile_threshold:float = 99.5
    threshold = np.percentile(data[(data['type']=="NegControlProbe") | (data['type']=="NegControlCodeword")]['count'].values,percentile_threshold)
    print('threshold = ', threshold)

    data['logfoldovernoise'] = data['count'].apply(lambda x: np.log(x / threshold))

    type2 = {"DeprecatedCodeword" : "DeprecatedCodeword",
            "NegControlCodeword" : 'NegControlCodeword',
            "NegControlProbe" : "NegControlProbe",
            "UnassignedCodeword" : "UnassignedCodeword"}

    data['type2'] = data['type'].map(type2)
    data['type2'] = data['type2'].fillna('Gene')

    categories = ['Gene', "DeprecatedCodeword","NegControlCodeword","NegControlProbe","UnassignedCodeword"]

    def add_jitter(x, scale=0.5):
        return x + np.random.uniform(-scale, scale, size=len(x))

    data['Jittered_Category'] = data['type2'].apply(lambda x: categories.index(x))
    data['Jittered_Category'] = add_jitter(data['Jittered_Category'])
    fig = plt.figure()
    plt.scatter(data['Jittered_Category'], data['logfoldovernoise'], alpha=0.5, s=1, color = 'black')
    plt.xticks(ticks=range(len(categories)), labels=categories)
    plt.ylabel('LogFold Over Noise')
    plt.xticks(rotation=45)
    plt.hlines(y=0, xmin=-0.5,xmax=4.5, linestyles='dashed', color = 'red', label='99.5 percentile')
    plt.title(f'LFG over noise - {sample}')
    plt.legend()
    plt.show()

    if save_plot:
        save_figure(fig,'LF_over_noise', name_dir=sample, dir_processed=dir_processed, format = "svg")