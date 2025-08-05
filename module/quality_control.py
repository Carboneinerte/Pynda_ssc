import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt



def plot_mcc_density (sample_1, save_plot = False, path_to_plot = 'plot', save_name = 'CC_mmc.png'):
    for sample in sample_1:
        df = pd.read_csv(f'D:\Jupyter_Notebook\Xenium_Jupyter_notebook\Correlation_Mapping\{sample}_CorrelationMapping.csv', comment = "#")
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
        plt.savefig(f'{path_to_plot}/{save_name}')


def desc_metrics(samples_ids, path_to_data, run_name = 'run_name', reference = True, plot_average = True, path_to_plot = 'plot', save_plot = False):
    if reference:
        reference_dataset = pd.read_csv('data/reference_dataset.csv')
    
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
        min_1 = 1
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
        fig.savefig(f'{path_to_plot}/{run_name}_QC.svg')


def desc_metrics_double(samples_ids_1, samples_ids_2, path_to_data_1, path_to_data_2, run_name = 'run_name', reference = True, plot_average = True,  path_to_plot = 'plot', save_plot = False):
    reference_dataset = pd.read_csv('data/reference_dataset.csv')
    
    parameters_to_plot = ['region_area', 'total_high_quality_decoded_transcripts','fraction_transcripts_decoded_q20',
                            'decoded_transcripts_per_100um2','estimated_number_of_false_positive_transcripts_per_cell','num_cells_detected',
                            'fraction_transcripts_assigned', 'median_genes_per_cell', 'median_transcripts_per_cell' ]
    
    for sample in samples_ids_1:
        with open(f"{path_to_data_1}\{sample}\metrics_summary.csv", 'r', encoding='utf-8') as file:
            file_content_1 = pd.read_csv(file)
            
            if 'files_content_1' in locals():
                files_content_1 = pd.concat([files_content_1,file_content_1])
            else:
                files_content_1 = pd.DataFrame(file_content_1)

    for sample in samples_ids_2:
        with open(f"{path_to_data_2}\{sample}\metrics_summary.csv", 'r', encoding='utf-8') as file:
            file_content_2 = pd.read_csv(file)
            
            if 'files_content_2' in locals():
                files_content_2 = pd.concat([files_content_2,file_content_2])
            else:
                files_content_2 = pd.DataFrame(file_content_2)

    if reference:
        min_1 = 1
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
        fig.savefig(f'{path_to_plot}/{run_name}_QC.svg')


