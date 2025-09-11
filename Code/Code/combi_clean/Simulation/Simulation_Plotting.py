import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import os


def find_latest_run_file(workspace_path_str, project_name, file_to_find):
    """
    Finds the most recent simulation output folder and returns the path to a specific file inside it.
    """
    workspace_path = Path(workspace_path_str)

    if not workspace_path.is_dir():
        print(f" Workspace directory not found at '{workspace_path}'")
        return None


    project_folders = [d for d in workspace_path.iterdir() if d.is_dir() and d.name.startswith(project_name)]

    if not project_folders:
        print(f" No simulation output folders found for project '{project_name}' in '{workspace_path}'")
        return None


    latest_folder = max(project_folders, key=os.path.getmtime)


    data_file_path = latest_folder / file_to_find

    if not data_file_path.is_file():
        print(f"Found latest run folder, but '{file_to_find}' is missing inside.")
        return None

    print(f"Found latest data file: {data_file_path}")
    return data_file_path



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

    print(f" Successfully saved master cytokine plot to: {output_plot_path} ")





def create_individual_cytokine_plots(file_path_str):
    """
    Loads cytokine summary data from a CSV and creates a separate, individual
    plot for each cytokine found in the data file.

    Args:
        file_path_str (str): The full path to the simulation_data.csv file.
    """
    file_path = Path(file_path_str)


    if not file_path.is_file():
        print(f"ERROR: Data file not found at: {file_path}")
        return

    print(f"\nCreating individual cytokine plots from: {file_path}")
    df = pd.read_csv(file_path)


    cytokines_to_plot = [
        'ccl2', 'damps', 'pamps', 'tgf_b', 'pdgf', 'fgf',
        'tnf_a', 'il1a', 'il1b', 'il6', 'il8', 'il10', 'il1ra'
    ]


    color_map = {
        'ccl2': '#1f77b4', 'damps': '#ff7f0e', 'pamps': '#2ca02c',
        'tgf_b': '#d62728', 'pdgf': '#9467bd', 'fgf': '#8c564b',
        'tnf_a': '#e377c2', 'il1a': '#7f7f7f', 'il1b': '#bcbd22',
        'il6': '#17becf', 'il8': '#aec7e8', 'il10': '#ffbb78',
        'il1ra': '#98df8a'
    }


    output_dir = file_path.parent / "individual_cytokine_plots"
    os.makedirs(output_dir, exist_ok=True)


    for cytokine in cytokines_to_plot:
        mean_col = f'{cytokine}_mean'
        std_col = f'{cytokine}_std'


        if mean_col not in df.columns or std_col not in df.columns:
            print(f"WARNING: Columns for '{cytokine}' not found in data file. Skipping.")
            continue


        sns.set_theme(style="whitegrid")
        plt.figure(figsize=(12, 7))

        color = color_map.get(cytokine, 'black')


        plt.fill_between(
            df['mcs'],
            df[mean_col] - df[std_col],
            df[mean_col] + df[std_col],
            alpha=0.3,
            color=color,
            label='Std. Dev.'
        )


        plt.plot(df['mcs'], df[mean_col],
                 label=f'Mean {cytokine.upper()}',
                 color=color,
                 linewidth=3)


        plt.title(f'Concentration of {cytokine.upper()} Over Time', fontsize=18)
        plt.xlabel('Monte Carlo Steps (MCS)', fontsize=12)
        plt.ylabel('Mean Concentration', fontsize=12)

        plt.yscale('log')

        plt.legend()
        plt.tight_layout()


        output_plot_path = output_dir / f'{cytokine}_plot.png'
        plt.savefig(output_plot_path)
        plt.show()
        plt.close()

        print(f" Successfully saved plot to: {output_plot_path}")


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

    print(f"Successfully saved cell count plot to: {output_plot_path}")





if __name__ == "__main__":


    if len(sys.argv) > 1:

        #incorrect because it uses the same plot  function for different things
        input_file = Path(sys.argv[1])
        if 'movement' in input_file.name:
            create_master_cytokine_plot(input_file)
        elif 'simulation' in input_file.name:
            create_master_cytokine_plot(input_file)
        else:
            print(f"Unsure how to plot file: {input_file.name}")
    else:


        CC3D_WORKSPACE_PATH = "E:/CompuCell3D-py3-64bit/Miniconda3/Lib/site-packages/cc3d/player5/"
        PROJECT_NAME = "endothelial_cc3d"

        latest_movement_file = find_latest_run_file(CC3D_WORKSPACE_PATH, PROJECT_NAME, "movement_data.csv")




        # Find and plot cytokine data (BOTH graphs)
        latest_cytokine_file = find_latest_run_file(CC3D_WORKSPACE_PATH, PROJECT_NAME, "simulation_data.csv")
        if latest_cytokine_file:
            create_master_cytokine_plot(latest_cytokine_file)

            ############turn here ON INDIVIDUAL PLOTSS!!#####
            # create_individual_cytokine_plots(latest_cytokine_file)

        latest_cell_count_file = find_latest_run_file(CC3D_WORKSPACE_PATH, PROJECT_NAME, "cell_counts.csv")
        if latest_cell_count_file:
            create_cell_count_plot(latest_cell_count_file)
            # create_daily_cell_violin_plots(latest_cell_count_file)




