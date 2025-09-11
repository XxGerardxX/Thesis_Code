import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path



def find_scan_iteration_file(scan_results_path_str, iteration_number, file_to_find):
    """
    Finds a specific data file inside a specific parameter scan iteration folder.
    """
    scan_results_path = Path(scan_results_path_str)
    iteration_folder = scan_results_path / f"scan_iteration_{iteration_number}"

    if not iteration_folder.is_dir():
        print(f"WARNING: Iteration folder not found: {iteration_folder}")
        return None

    try:

        data_file_path = next(iteration_folder.glob(f'**/{file_to_find}'))
        print(f"Found data file for iteration {iteration_number}: {data_file_path}")
        return data_file_path
    except StopIteration:
        print(f"WARNING: Could not find '{file_to_find}' in {iteration_folder}")
        return None




def find_all_scan_files(scan_results_path_str, file_to_find):
    """
    Scans a parameter scan results folder and returns a list of paths
    to the specified data file in EVERY iteration subfolder.
    """
    scan_results_path = Path(scan_results_path_str)
    if not scan_results_path.is_dir():
        print(f"ERROR: Scan results directory not found at '{scan_results_path}'")
        return []


    file_paths = list(scan_results_path.glob(f'**/{file_to_find}'))

    if not file_paths:
        print(f"ERROR: No '{file_to_find}' files found in '{scan_results_path}' or its subdirectories.")
        return []

    print(f"Found {len(file_paths)} copies of '{file_to_find}'.")
    return file_paths


def create_master_cytokine_plot(file_path):
    """
    Loads cytokine summary data and creates a single, clear plot of all
    cytokine mean concentrations over time, with distinct colors.
    """
    print(f"Loading cytokine data from: {file_path}")
    df = pd.read_csv(file_path)


    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]



    color_map = {
        'ccl2': '#1f77b4',
        'damps': '#ff7f0e',
        'pamps': '#2ca02c',
        'tgf_b': '#d62728',
        'pdgf': '#9467bd',
        'fgf': '#8c564b',
        'tnf_a': '#e377c2',
        'il1a': '#7f7f7f',
        'il1b': '#bcbd22',
        'il6': '#17becf',
        'il8': '#aec7e8',
        'il10': '#ffbb78',
        'il1ra': '#98df8a'
    }

    sns.set_theme(style="whitegrid")

    PLOT_SMOOTHED_LINE = True
    smoothing_window = 5
    plt.figure(figsize=(16, 9))


    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'

        if mean_col in df.columns:
            color = color_map.get(cytokine, 'black')

            if PLOT_SMOOTHED_LINE:

                smoothed_col = f'{cytokine}_smoothed'
                df[smoothed_col] = df[mean_col].rolling(window=smoothing_window, center=True).mean()


                plt.plot(df['mcs'], df[mean_col], color=color, alpha=0.2)


                plt.plot(df['mcs'], df[smoothed_col], label=cytokine.upper(), color=color, linewidth=2.0)

            else:

                plt.plot(df['mcs'], df[mean_col], label=cytokine.upper(), color=color, linewidth=1.5)

        else:
            print(f"WARNING: Column '{mean_col}' not found. Skipping.")

    title_text = 'Mean Cytokine Concentrations Over Time'
    if PLOT_SMOOTHED_LINE:
        title_text += ' (Smoothed)'

    plt.title(title_text, fontsize=20)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Concentration (mg/ml)', fontsize=14)
    plt.yscale('log')
    plt.legend(title="Cytokine", ncol=2)
    plt.tight_layout()



    output_plot_path = file_path.parent / 'master_cytokine_plot.png'
    plt.savefig(output_plot_path)
    plt.show()
    plt.close()

    print(f" Successfully saved master cytokine plot to: {output_plot_path}")


def create_master_cytokine_plot_v_2(file_path):
    """
    Loads cytokine summary data from a SINGLE RUN, converts units, and creates
    a smoothed time-series plot with consistent colors against MCS.
    """
    print(f" Generating Single-Run Smoothed Cytokine Plot from: {file_path} ")
    df = pd.read_csv(file_path)


    VOXEL_VOLUME_ML = 8.0e-6

    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]
    color_map = {
        'ccl2': '#1f77b4', 'damps': '#ff7f0e', 'pamps': '#2ca02c', 'tgf_b': '#d62728',
        'pdgf': '#9467bd', 'fgf': '#8c564b', 'tnf_a': '#e377c2', 'il1a': '#7f7f7f',
        'il1b': '#bcbd22', 'il6': '#17becf', 'il8': '#aec7e8', 'il10': '#ffbb78', 'il1ra': '#98df8a'
    }




    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        if mean_col in df.columns:
            df[mean_col] /= VOXEL_VOLUME_ML


    sns.set_theme(style="whitegrid")
    PLOT_SMOOTHED_LINE = True
    smoothing_window = 5
    plt.figure(figsize=(16, 9))

    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        if mean_col in df.columns:
            color = color_map.get(cytokine, 'black')

            if PLOT_SMOOTHED_LINE:
                smoothed_col = f'{cytokine}_smoothed'
                df[smoothed_col] = df[mean_col].rolling(window=smoothing_window, center=True).mean()


                plt.plot(df['mcs'], df[mean_col], color=color, alpha=0.2)
                plt.plot(df['mcs'], df[smoothed_col], label=cytokine.upper(), color=color, linewidth=2.0)
            else:
                plt.plot(df['mcs'], df[mean_col], label=cytokine.upper(), color=color, linewidth=1.5)
        else:
            print(f"WARNING: Column '{mean_col}' not found. Skipping.")


    title_text = 'Cytokine Concentrations Over Time (Single Run, Smoothed)'
    plt.title(title_text, fontsize=20)

    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Concentration (pg/ml)', fontsize=14)
    plt.yscale('log')
    plt.legend(title="Cytokine", ncol=2)
    plt.tight_layout()


    output_plot_path = file_path.parent / 'single_run_cytokine_plot_mcs.png'
    plt.savefig(output_plot_path)
    plt.show()
    plt.close()

    print(f"Successfully saved single-run cytokine plot to: {output_plot_path} ")




def create_master_cytokine_plot_with_SD(file_path):
    """
    Loads cytokine summary data from a SINGLE RUN, converts units, and creates
    a time-series plot showing the mean concentration with a shaded region
    representing the spatial standard deviation (within-run variability).
    """
    print(f"Generating Single-Run Cytokine Plot with Spatial SD from: {file_path} ")
    df = pd.read_csv(file_path)


    VOXEL_VOLUME_ML = 8.0e-6

    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]
    color_map = {
        'ccl2': '#1f77b4', 'damps': '#ff7f0e', 'pamps': '#2ca02c', 'tgf_b': '#d62728',
        'pdgf': '#9467bd', 'fgf': '#8c564b', 'tnf_a': '#e377c2', 'il1a': '#7f7f7f',
        'il1b': '#bcbd22', 'il6': '#17becf', 'il8': '#aec7e8', 'il10': '#ffbb78', 'il1ra': '#98df8a'
    }

    # Convert BOTH mean and std deviation columns to pg/ml
    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        std_col = f'{cytokine}_std'
        if mean_col in df.columns:
            df[mean_col] /= VOXEL_VOLUME_ML
        if std_col in df.columns:
            df[std_col] /= VOXEL_VOLUME_ML


    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(16, 9))

    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        std_col = f'{cytokine}_std'


        if mean_col in df.columns and std_col in df.columns:
            color = color_map.get(cytokine, 'black')


            plt.plot(df['mcs'], df[mean_col], label=cytokine.upper(), color=color, linewidth=2.0)


            plt.fill_between(
                df['mcs'],
                df[mean_col] - df[std_col],
                df[mean_col] + df[std_col],
                color=color,
                alpha=0.2
            )
        else:
            print(f"WARNING: Columns for '{cytokine}' not found. Skipping.")




    title_text = 'Cytokine Concentrations with Spatial Standard Deviation (Single Run)'
    plt.title(title_text, fontsize=20)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Concentration (pg/ml)', fontsize=14)
    plt.yscale('log')
    plt.legend(title="Cytokine", ncol=2)
    plt.tight_layout()


    output_plot_path = file_path.parent / 'single_run_cytokine_plot_with_spatial_std.png'
    plt.savefig(output_plot_path)
    plt.show()
    plt.close()

    print(f" Successfully saved single-run cytokine plot to: {output_plot_path}")



def create_cell_count_plot(file_path):
    """
    Loads cell count data from a CSV and creates a clear plot of
    population dynamics over time with distinct, consistent colors.
    """

    CELL_COLOR_MAP = {
        "endothelial": "darkgreen",
        "neutrophil": "blue",
        "monocyte": "purple",
        "fibroblast": "saddlebrown",
        "neutrophila": "cyan",
        "neutrophilnec": "darkblue",
        "mast": "magenta",
        "macrophage1": "red",
        "macrophage2": "orange",
        "keratino": "gold",
        "platelet": "lightcoral",
        "necrotic": "black",
        "temp": "lime",
        "nectemp": "dimgray"
    }

    CELL_TYPE_MAP = {
        1: "endothelial", 2: "neutrophil", 3: "monocyte", 4: "fibroblast",
        5: "neutrophila", 6: "neutrophilnec", 7: "mast", 8: "macrophage1",
        9: "macrophage2", 10: "keratino", 11: "platelet", 12: "necrotic",
        13: "temp", 14: "nectemp"
    }

    print(f"Loading cell count data from: {file_path}")
    df = pd.read_csv(file_path)

    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(14, 8))


    cell_types_in_data = sorted(df['CellType'].unique())

    custom_palette = {}
    for cell_type_id in cell_types_in_data:

        cell_name = CELL_TYPE_MAP.get(cell_type_id, "Unknown")

        color = CELL_COLOR_MAP.get(cell_name, "gray")

        custom_palette[cell_type_id] = color


    sns.lineplot(
        data=df,
        x="MCS",
        y="Count",
        hue="CellType",
        palette=custom_palette,
        linewidth=2.5
    )


    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()

    new_labels = [CELL_TYPE_MAP.get(int(l), f"Unknown {l}") for l in labels[1:]]
    ax.legend(handles=handles[1:], labels=new_labels, title="Cell Type", ncol=2)


    plt.title('Cell Population Dynamics Over Time', fontsize=18)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Total Cell Count (Log Scale)', fontsize=14)

    plt.tight_layout()

    output_plot_path = file_path.with_suffix('.png')
    plt.savefig(output_plot_path)
    plt.show()
    plt.close()

    print(f" Successfully saved cell count plot to: {output_plot_path}")



def plot_aggregated_cell_counts(scan_results_path_str):
    """
    Loads cell count data from all scan iterations, calculates the mean and
    confidence interval, and plots the aggregated results versus MCS.
    (This version excludes progenitor cells and has an enhanced color map for clarity).
    """
    print("\n Processing Aggregated Scan Results for: Cell Counts ")

    all_file_paths = find_all_scan_files(scan_results_path_str, "cell_counts.csv")
    if not all_file_paths: return

    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)



    CELL_COLOR_MAP = {
        "endothelial": "darkgreen", "neutrophil": "blue", "monocyte": "purple",
        "fibroblast": "lime", "neutrophila": "cyan", "neutrophilnec": "darkblue",
        "mast": "magenta", "macrophage1": "red", "macrophage2": "orange",
        "keratino": "gold", "platelet": "lightcoral", "necrotic": "black",
        "temp": "gray", "nectemp": "gray"
    }
    CELL_TYPE_MAP = {
        1: "endothelial", 2: "neutrophil", 3: "monocyte", 4: "fibroblast",
        5: "neutrophila", 6: "neutrophilnec", 7: "mast", 8: "macrophage1",
        9: "macrophage2", 10: "keratino", 11: "platelet", 12: "necrotic",
        13: "temp", 14: "nectemp"
    }
    combined_df['CellName'] = combined_df['CellType'].map(CELL_TYPE_MAP)


    cell_types_to_plot = [
        "endothelial", "neutrophil", "monocyte", "fibroblast", "neutrophila",
        "neutrophilnec", "mast", "macrophage1", "macrophage2", "keratino",
        "platelet", "necrotic"
    ]
    plotting_df = combined_df[combined_df['CellName'].isin(cell_types_to_plot)]


    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")

    sns.lineplot(
        data=plotting_df,
        x='MCS',
        y='Count',
        hue='CellName',
        palette=CELL_COLOR_MAP
    )


    plt.title('Mean Cell Population Dynamics (with Confidence Interval)', fontsize=20)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Cell Count', fontsize=14)
    # plt.yscale('log')
    plt.legend(title='Cell Type', ncol=2)
    plt.tight_layout()


    output_path = Path(scan_results_path_str) / "Aggregated_Cell_Counts_Plot_MCS.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved aggregated cell count plot to: {output_path} ")



def plot_aggregated_cytokines(scan_results_path_str):
    """
    Loads cytokine data from all scan iterations, converts to pg/ml, calculates
    mean and std dev, and plots the aggregated results.
    """
    print("\n Processing Aggregated Scan Results for: Cytokines ")

    all_file_paths = find_all_scan_files(scan_results_path_str, "simulation_data.csv")
    if not all_file_paths: return

    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)


    VOXEL_VOLUME_ML = 8.0e-6


    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]
    color_map = {
        'ccl2': '#1f77b4', 'damps': '#ff7f0e', 'pamps': '#2ca02c', 'tgf_b': '#d62728',
        'pdgf': '#9467bd', 'fgf': '#8c564b', 'tnf_a': '#e377c2', 'il1a': '#7f7f7f',
        'il1b': '#bcbd22', 'il6': '#17becf', 'il8': '#aec7e8', 'il10': '#ffbb78', 'il1ra': '#98df8a'
    }


    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        std_col = f'{cytokine}_std'
        if mean_col in combined_df.columns:
            combined_df[mean_col] = combined_df[mean_col] / VOXEL_VOLUME_ML
            combined_df[std_col] = combined_df[std_col] / VOXEL_VOLUME_ML


    MCS_PER_DAY = 24000.0
    combined_df['Days'] = combined_df['mcs'] / MCS_PER_DAY



    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()


    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        if mean_col not in combined_df.columns: continue

        agg_df = combined_df.groupby('Days')[mean_col].agg(['mean', 'std']).reset_index()

        color = color_map.get(cytokine, 'black')

        ax.plot(agg_df['Days'], agg_df['mean'], label=cytokine.upper(), color=color, linewidth=2)
        ax.fill_between(
            agg_df['Days'],
            agg_df['mean'] - agg_df['std'],
            agg_df['mean'] + agg_df['std'],
            color=color,
            alpha=0.2
        )


    plt.title('Mean Cytokine Concentrations (with Std. Dev. Across Runs)', fontsize=20)
    plt.xlabel('Time (Days)', fontsize=14)
    plt.ylabel('Mean Concentration (pg/ml)', fontsize=14)
    plt.yscale('log')
    plt.legend(title="Cytokine", ncol=2)
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_Cytokine_Plot_pg_ml.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved aggregated cytokine plot (in pg/ml) to: {output_path} ")


def plot_aggregated_il1ra(scan_results_path_str):
    """
    Loads IL-1RA data from all scan iterations, calculates the mean and
    standard deviation, and plots only the aggregated IL-1RA result.
    """
    print("\nProcessing Aggregated Scan Results for: IL-1RA")


    all_file_paths = find_all_scan_files(scan_results_path_str, "simulation_data.csv")
    if not all_file_paths:
        print("ERROR: No data files found. Cannot generate plot.")
        return


    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)


    cytokine_to_plot = 'il1ra'
    mean_col = f'{cytokine_to_plot}_mean'


    if mean_col not in combined_df.columns:
        print(f"ERROR: Column '{mean_col}' not found in the data files.")
        return


    agg_df = combined_df.groupby('mcs')[mean_col].agg(['mean', 'std']).reset_index()


    plt.figure(figsize=(14, 8))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()

    il1ra_color = '#98df8a'


    ax.plot(agg_df['mcs'], agg_df['mean'], label='Mean IL-1RA', color=il1ra_color, linewidth=3)


    ax.fill_between(
        agg_df['mcs'],
        agg_df['mean'] - agg_df['std'],
        agg_df['mean'] + agg_df['std'],
        color=il1ra_color,
        alpha=0.3,
        label='Std. Dev.'
    )


    plt.title('Mean IL-1RA Concentration (with Std. Dev. Across Runs)', fontsize=20)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Concentration (mg/ml)', fontsize=14)
    plt.yscale('log')
    plt.legend()
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_IL1RA_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved aggregated IL-1RA plot to: {output_path} ")





def plot_aggregated_cell_movement(scan_results_path_str):
    """
    Loads cell movement data from all scan iterations, calculates the mean and
    confidence interval for distance traveled, and plots the results.
    """
    print("\n Processing Aggregated Scan Results for: Cell Movement ")

    all_file_paths = find_all_scan_files(scan_results_path_str, "movement_data.csv")
    if not all_file_paths:
        print("INFO: No movement_data.csv files found. Skipping plot.")
        return


    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)


    CELL_COLOR_MAP = {
        "endothelial": "darkgreen", "neutrophil": "blue", "monocyte": "purple",
        "fibroblast": "saddlebrown", "neutrophila": "cyan", "neutrophilnec": "darkblue",
        "mast": "magenta", "macrophage1": "red", "macrophage2": "orange",
        "keratino": "gold", "platelet": "lightcoral", "necrotic": "black",
        "temp": "lime", "nectemp": "dimgray"
    }
    CELL_TYPE_MAP = {
        1: "endothelial", 2: "neutrophil", 3: "monocyte", 4: "fibroblast",
        5: "neutrophila", 6: "neutrophilnec", 7: "mast", 8: "macrophage1",
        9: "macrophage2", 10: "keratino", 11: "platelet", 12: "necrotic",
        13: "temp", 14: "nectemp"
    }
    combined_df['CellName'] = combined_df['CellType'].map(CELL_TYPE_MAP)


    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")


    sns.lineplot(
        data=combined_df,
        x='MCS',
        y='Distance',
        hue='CellName',
        palette=CELL_COLOR_MAP
    )


    plt.title('Mean Cell Motility (with Confidence Interval)', fontsize=20)
    plt.xlabel('Monte Carlo Steps (MCS)', fontsize=14)
    plt.ylabel('Mean Distance Traveled per Step', fontsize=14)
    plt.legend(title='Cell Type', ncol=2)
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_Movement_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved aggregated cell movement plot to: {output_path} ")


def create_il1ra_comprehensive_validation_plot(scan_results_path_str):
    """
    Loads all IL-1RA scan data, calculates mean & std dev, and plots it against
    a comprehensive set of literature data points for validation.

    """
    print("\n Creating Comprehensive IL-1RA Validation Plot ")


    all_file_paths = find_all_scan_files(scan_results_path_str, "simulation_data.csv")
    if not all_file_paths:
        print("ERROR: No data files found for the scan.")
        return

    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)


    mean_col = 'il1ra_mean'
    if mean_col not in combined_df.columns:
        print(f"ERROR: Column '{mean_col}' not found.")
        return

    agg_df = combined_df.groupby('mcs')[mean_col].agg(['mean', 'std']).reset_index()


    agg_df['Days'] = agg_df['mcs'] / 24000.0
    VOXEL_VOLUME_ML = 8.0e-6
    agg_df['Mean_pg_ml'] = agg_df['mean'] / VOXEL_VOLUME_ML
    agg_df['Std_Dev_pg_ml'] = agg_df['std'] / VOXEL_VOLUME_ML


    literature_data = {
        'Vindenes Fig 2 (Day 0)': {'type': 'MeanSEM', 'days': 0, 'mean': 3199, 'error': 400},
        'Vindenes Fig 2 (Day 2)': {'type': 'MeanSEM', 'days': 2, 'mean': 2500, 'error': 200},
        'Vindenes Fig 2 (Day 5)': {'type': 'MeanSEM', 'days': 5, 'mean': 2600, 'error': 200},
        'Vindenes Fig 2 (Day 10)': {'type': 'MeanSEM', 'days': 10, 'mean': 2100, 'error': 400},
        'Vindenes Text (Day 0)': {'type': 'MeanMinMax', 'days': 0, 'mean': 3199, 'min': 1018, 'max': 9436},
        'Hur (Day 1)': {'type': 'MinMax', 'days': 1, 'min': 120.5, 'max': 10067.7},
        'Hur (Day 3)': {'type': 'MinMax', 'days': 3, 'min': 64.7, 'max': 1552.5}
    }


    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()


    model_color = 'steelblue'
    sem_color = 'salmon'
    range_color = 'mediumpurple'


    ax.plot(agg_df['Days'], agg_df['Mean_pg_ml'], label='Model Mean', color=model_color, linewidth=3, zorder=5)
    ax.fill_between(
        agg_df['Days'],
        agg_df['Mean_pg_ml'] - agg_df['Std_Dev_pg_ml'],
        agg_df['Mean_pg_ml'] + agg_df['Std_Dev_pg_ml'],
        color=model_color, alpha=0.2, label='Model Std. Dev.'
    )


    ax.errorbar([], [], yerr=[], fmt='s', color=sem_color, markersize=8, capsize=5, label='Literature (Mean ± SEM)')
    ax.errorbar([], [], yerr=[], fmt='D', color=range_color, markersize=8, capsize=5, label='Literature (Full Range)')

    for point_name, data in literature_data.items():
        if data['type'] == 'MeanSEM':
            ax.errorbar(
                x=data['days'], y=data['mean'], yerr=data['error'],
                fmt='s', color=sem_color, markersize=8, capsize=5, elinewidth=2,
                markeredgecolor='black', zorder=10
            )
        elif data['type'] in ['MeanMinMax', 'MinMax']:
            center_point = data.get('mean', (data['max'] + data['min']) / 2)
            lower_error = center_point - data['min']
            upper_error = data['max'] - center_point
            ax.errorbar(
                x=data['days'], y=center_point, yerr=[[lower_error], [upper_error]],
                fmt='D', color=range_color, markersize=8, capsize=5, elinewidth=2,
                markeredgecolor='black', zorder=10
            )




    ax.set_title('Model Validation: IL-1Ra Dynamics vs. Published Clinical Data', fontsize=20)
    ax.set_xlabel('Time (Days post-injury)', fontsize=14)
    ax.set_ylabel('IL-1Ra Concentration (pg/ml)', fontsize=14)
    ax.set_yscale('log')
    ax.set_xlim(-0.5, 12.5)
    ax.legend()
    plt.tight_layout()


    output_path = Path(scan_results_path_str) / "Aggregated_IL1RA_Validation_Plot_Pastel.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved comprehensive validation plot to: {output_path}")



def plot_spatial_variability(scan_results_path_str):
    """
    Plots the mean within-run (spatial) standard deviation for ALL cytokines
    over time to show how field heterogeneity evolves, using dashed lines.
    """
    print("\nGenerating Spatial Variability Plot ")

    all_file_paths = find_all_scan_files(scan_results_path_str, "simulation_data.csv")
    if not all_file_paths: return

    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)


    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]


    color_map = {
        'ccl2': '#1f77b4', 'damps': '#ff7f0e', 'pamps': '#2ca02c', 'tgf_b': '#d62728',
        'pdgf': '#9467bd', 'fgf': '#8c564b', 'tnf_a': '#e377c2', 'il1a': '#7f7f7f',
        'il1b': '#bcbd22', 'il6': '#17becf', 'il8': '#aec7e8', 'il10': '#ffbb78', 'il1ra': '#98df8a'
    }


    VOXEL_VOLUME_ML = 8.0e-6
    for cytokine in cytokines_to_plot:
        std_col = f'{cytokine}_std'
        if std_col in combined_df.columns:
            combined_df[std_col] /= VOXEL_VOLUME_ML

    combined_df['Days'] = combined_df['mcs'] / 24000.0


    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()

    for cytokine in cytokines_to_plot:
        std_col = f'{cytokine}_std'
        if std_col in combined_df.columns:

            agg_df = combined_df.groupby('Days')[std_col].agg('mean').reset_index()
            color = color_map.get(cytokine, 'black')


            ax.plot(agg_df['Days'], agg_df[std_col], label=cytokine.upper(),
                    linestyle='--', linewidth=2.5, color=color)


    ax.set_title('Spatial Heterogeneity of Cytokine Fields Over Time', fontsize=20)
    ax.set_xlabel('Time (Days)', fontsize=14)
    ax.set_ylabel('Mean Spatial Std. Dev. (pg/ml)', fontsize=14)
    ax.set_yscale('log')
    ax.legend(title='Cytokine', ncol=2)
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_Spatial_Variability_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved spatial variability plot to: {output_path}")


def plot_macrophage_polarization(scan_results_path_str):
    """
    Calculates and plots the M2/M1 macrophage ratio over time to show polarization dynamics.
    """
    print("\n Generating Macrophage Polarization Plot")

    cell_files = find_all_scan_files(scan_results_path_str, "cell_counts.csv")
    if not cell_files: return

    list_of_dfs = [pd.read_csv(f) for f in cell_files]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)
    combined_df['Days'] = combined_df['MCS'] / 24000.0


    m1_counts = combined_df[combined_df['CellType'] == 8].groupby('Days')['Count'].mean()
    m2_counts = combined_df[combined_df['CellType'] == 9].groupby('Days')['Count'].mean()
    mac_ratio = m2_counts / (m1_counts + 1e-9)


    plt.figure(figsize=(12, 7))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()

    ax.plot(mac_ratio.index, mac_ratio.values, color='purple', linewidth=2.5)
    ax.axhline(y=1.0, color='grey', linestyle='--', label='Balance Point')

    ax.set_title('Macrophage Polarization Dynamics', fontsize=20)
    ax.set_xlabel('Time (Days)', fontsize=14)
    ax.set_ylabel('Mean M2 / M1 Ratio', fontsize=14)
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_Macrophage_Polarization_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved macrophage polarization plot to: {output_path}")


def plot_inflammatory_balance(scan_results_path_str):
    """
    Calculates and plots the composite pro-/anti-inflammatory cytokine ratio over time.
    """
    print("\n Generating Inflammatory Balance Plot")

    cytokine_files = find_all_scan_files(scan_results_path_str, "simulation_data.csv")
    if not cytokine_files: return

    list_of_dfs = [pd.read_csv(f) for f in cytokine_files]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)
    combined_df['Days'] = combined_df['mcs'] / 24000.0


    VOXEL_VOLUME_ML = 8.0e-6
    for col in combined_df.columns:
        if '_mean' in col or '_std' in col:
            combined_df[col] = combined_df[col] / VOXEL_VOLUME_ML


    pro_inflammatory_cols = ['ccl2_mean', 'tnf_a_mean', 'il8_mean', 'il6_mean', 'il1a_mean', 'il1b_mean', 'damps_mean',
                             'pamps_mean']
    anti_inflammatory_cols = ['il10_mean', 'tgf_b_mean', 'pdgf_mean', 'fgf_mean', 'il1ra_mean']
    pro_inflammatory_sum = combined_df.groupby('Days')[pro_inflammatory_cols].mean().sum(axis=1)
    anti_inflammatory_sum = combined_df.groupby('Days')[anti_inflammatory_cols].mean().sum(axis=1)
    cytokine_ratio = pro_inflammatory_sum / (anti_inflammatory_sum + 1e-9)


    plt.figure(figsize=(12, 7))
    sns.set_theme(style="whitegrid")
    ax = plt.gca()

    ax.plot(cytokine_ratio.index, cytokine_ratio.values, color='darkturquoise', linewidth=2.5)
    ax.axhline(y=1.0, color='grey', linestyle='--', label='Balance Point')

    ax.set_title('Composite Inflammatory Balance', fontsize=20)
    ax.set_xlabel('Time (Days)', fontsize=14)
    ax.set_ylabel('Pro- / Anti-inflammatory Ratio', fontsize=14)
    ax.set_yscale('log')
    ax.legend()
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Aggregated_Inflammatory_Balance_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved inflammatory balance plot to: {output_path} ")



def plot_estimated_biological_counts(scan_results_path_str):
    """
    Loads cell count data, and using the principle of coarse-graining, converts
    the agent counts into an estimated number of real biological cells based on volume.
    """
    print("\nGenerating Estimated Biological Cell Count Plot ")


    all_file_paths = find_all_scan_files(scan_results_path_str, "cell_counts.csv")
    if not all_file_paths: return

    list_of_dfs = [pd.read_csv(f) for f in all_file_paths]
    combined_df = pd.concat(list_of_dfs, ignore_index=True)




    VOXEL_VOLUME_UM3 = 8.0e6


    AGENT_TARGET_VOLUMES = {
        1: 25, 2: 20, 3: 30, 4: 25, 5: 20, 6: 20, 7: 23,
        8: 35, 9: 35, 10: 25, 11: 18, 12: 25, 13: 25, 14: 25
    }



    necrotic_mean_vol = np.mean([1600, 8, 2500, 300])
    REAL_CELL_VOLUMES = {
        1: 1000, 2: 300, 3: 3000, 4: 1600, 5: 300, 6: 300, 7: 1400,
        8: 4990, 9: 4990, 10: 2500, 11: 8, 12: necrotic_mean_vol,
        13: 1000, 14: necrotic_mean_vol
    }


    combined_df['AgentVolume_voxels'] = combined_df['CellType'].map(AGENT_TARGET_VOLUMES)
    combined_df['RealCellVolume_um3'] = combined_df['CellType'].map(REAL_CELL_VOLUMES)


    combined_df['TotalOccupiedVolume_um3'] = combined_df['Count'] * combined_df['AgentVolume_voxels'] * VOXEL_VOLUME_UM3


    combined_df['EstimatedBiologicalCount'] = combined_df['TotalOccupiedVolume_um3'] / combined_df['RealCellVolume_um3']


    combined_df['Days'] = combined_df['MCS'] / 24000.0
    CELL_COLOR_MAP = {
        "endothelial": "darkgreen", "neutrophil": "blue", "monocyte": "purple",
        "fibroblast": "lime", "neutrophila": "cyan", "neutrophilnec": "darkblue",
        "mast": "magenta", "macrophage1": "red", "macrophage2": "orange",
        "keratino": "gold", "platelet": "lightcoral", "necrotic": "black",
        "temp": "gray", "nectemp": "gray"
    }
    CELL_TYPE_MAP = {
        1: "endothelial", 2: "neutrophil", 3: "monocyte", 4: "fibroblast",
        5: "neutrophila", 6: "neutrophilnec", 7: "mast", 8: "macrophage1",
        9: "macrophage2", 10: "keratino", 11: "platelet", 12: "necrotic",
        13: "temp", 14: "nectemp"
    }
    combined_df['CellName'] = combined_df['CellType'].map(CELL_TYPE_MAP)

    cell_types_to_plot = [
        "endothelial", "neutrophil", "monocyte", "fibroblast", "neutrophila",
        "neutrophilnec", "mast", "macrophage1", "macrophage2", "keratino",
        "platelet", "necrotic"
    ]
    plotting_df = combined_df[combined_df['CellName'].isin(cell_types_to_plot)]


    plt.figure(figsize=(16, 9))
    sns.set_theme(style="whitegrid")

    sns.lineplot(
        data=plotting_df,
        x='Days',
        y='EstimatedBiologicalCount',
        hue='CellName',
        palette=CELL_COLOR_MAP
    )


    plt.title('Estimated Biological Cell Population Dynamics', fontsize=20)
    plt.xlabel('Time (Days)', fontsize=14)
    plt.ylabel('Estimated Number of Biological Cells', fontsize=14)
    plt.yscale('log')
    plt.legend(title='Cell Type', ncol=2)
    plt.tight_layout()

    output_path = Path(scan_results_path_str) / "Estimated_Biological_Counts_Plot.png"
    plt.savefig(output_path)
    plt.show()
    plt.close()

    print(f" Saved estimated biological cell count plot to: {output_path}")

if __name__ == "__main__":


    "This is the way you plot your data. Remove the string and put your own path"

    SCAN_RESULTS_DIRECTORY_IL1RA = r"C:\Users\Gerard\Desktop\Parameter scan\IL_1RA_Results_k2_secretion"
    NUMBER_OF_ITERATIONS_IL1RA = 5  # You ran from 0 to 4

    SCAN_RESULTS_DIRECTORY_IL1RA_MU_FIRST = r"C:\Users\Gerard\Desktop\Parameter scan\Il_1RA_mu_results"
    NUMBER_OF_ITERATIONS_IL1RA_MU_FIRST = 6

    SCAN_RESULTS_DIRECTORY_IL1RA_MU_SECOND = r"C:\Users\Gerard\Desktop\Parameter scan\Il_1RA_mu_specific_Results"
    NUMBER_OF_ITERATIONS_IL1RA_MU_SECOND = 6

    SCAN_RESULTS_DIRECTORY_IL1RA_MU_FINAL = r"C:\Users\Gerard\Desktop\Parameter scan\Result_Final_Mu"
    NUMBER_OF_ITERATIONS_IL1RA_MU_FINAL = 7

    # this was to test when you switched the TNF-a or something else to CCl2 because of the plot you had
    SCAN_RESULTS_DIRECTORY = r"C:\Users\Gerard\Desktop\Parameter scan\CCL2_Results"

    SCAN_RESULTS_DIRECTORY_IL1RA_DIFF = r"C:\Users\Gerard\Desktop\Parameter scan\IL1RA_DIFF_result_Logarithmic"
    NUMBER_OF_ITERATIONS_IL1RA_DIFF = 5

    SCAN_RESULTS_DIRECTORY_START_M2 = r"C:\Users\Gerard\Desktop\Parameter scan\Starting_amount_results"

    SCAN_RESULTS = r"C:\Users\Gerard\Desktop\Parameter scan\3DParam_Scan_Roshan"

    SCAN_RESULTS_FIBRO = r"C:\Users\Gerard\Desktop\Code_8_full_finalized_model\Temp_Test_Fib_Result"

    SCAN_RESULTS_FIBRO_TEMP = r"C:\Users\Gerard\Desktop\Code_8_full_finalized_model\Temp_Temperature_Test_Fibro_deletion_RESULT"
    NUMBER_OF_ITERATIONS_FIBRO_TEMP = 3

    SCAN_RESULTS_DIRECTORY_Number_TWO = r"C:\Users\Gerard\Desktop\Parameter scan\3DParam_Scan_Roshan2"




    "This is used for multiple plots for different simulations (e.g. Appendix A)"
    # for i in range(NUMBER_OF_ITERATIONS_IL1RA):
    #
    #
    #     # Find and plot cytokine data for IL-1RA
    #     cytokine_file = find_scan_iteration_file(SCAN_RESULTS_DIRECTORY_IL1RA, i, "simulation_data.csv")
    #     if cytokine_file:
    #         create_master_cytokine_plot(cytokine_file)
    #
    #     # Find and plot cell count data for IL-1RA
    #     cell_count_file = find_scan_iteration_file(SCAN_RESULTS_DIRECTORY_IL1RA, i, "cell_counts.csv")
    #     if cell_count_file:
    #         create_cell_count_plot(cell_count_file)

    # for i in range(NUMBER_OF_ITERATIONS_IL1RA_MU_FIRST):
    #     cytokine_file = find_scan_iteration_file(SCAN_RESULTS_DIRECTORY_IL1RA_MU_FIRST, i, "simulation_data.csv")
    #     if cytokine_file:
    #         # create_master_cytokine_plot(cytokine_file)
    #         create_master_cytokine_plot_v_2(cytokine_file)
    #         create_master_cytokine_plot_with_SD(cytokine_file)
    #     else:
    #         continue
    #
    # for i in range(NUMBER_OF_ITERATIONS_IL1RA_MU_SECOND):
    #     # Find and plot cell count data for IL-1RA
    #     cytokine_file = find_scan_iteration_file(SCAN_RESULTS_DIRECTORY_IL1RA_MU_SECOND, i, "simulation_data.csv")
    #     if cytokine_file:
    #         create_master_cytokine_plot_v_2(cytokine_file)
    #         create_master_cytokine_plot_with_SD(cytokine_file)
    #     else:
    #         continue
    #




    # Call the new function to plot only IL-1RA
    # plot_aggregated_il1ra(SCAN_RESULTS_DIRECTORY)
    #
    # # 2. Create the aggregated cytokine plot with SD
    # plot_aggregated_cytokines(SCAN_RESULTS_DIRECTORY)
    # plot_aggregated_cell_counts(SCAN_RESULTS_DIRECTORY)


