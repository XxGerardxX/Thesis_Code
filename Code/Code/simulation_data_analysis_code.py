import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns
import os
import matplotlib.ticker as mticker
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FixedLocator, FixedFormatter, MaxNLocator
from scipy.stats import zscore


sns.set(rc={'figure.figsize':(11.7,8.27)})

def cell_count(file, exp, show, outputdir=None):
    # Load data from the CSV file
    data = pd.read_csv(file, sep=",", usecols=["mcsteps", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])

    # Cell count labels and colors for the first plot
    cell_count_labels = ["Endothelial", "Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    colors = ["blue", "brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]

    # Cell count labels and colors for the second plot
    scoped_cols = ["1", "4", "5", "6", "7", "8", "9"]
    scoped_colors = ["blue", "violet", "red", "pink", "yellow", "orange", "darkblue"]
    scoped_labels = ["Endothelial", "Fibroblasts", "Activated Neutrophils",
                     "Necrotic Neutrophils", "Resting Monocytes", "Macrophages I", "Macrophages II"]

    # Create subplots with adjusted layout
    fig, axs = plt.subplots(ncols=1, nrows=2)
    plt.subplots_adjust(left=0.07, bottom=0.075, right=0.800, top=0.96)

    # Plot the data for the first plot
    for i, col in enumerate(cols):
        sns.lineplot(x="mcsteps", y=data[col], data=data, color=colors[i], label=cell_count_labels[i], ax=axs[0])

    # Plot the data for the second plot
    for i, col in enumerate(scoped_cols):
        sns.lineplot(x="mcsteps", y=data[col], data=data, color=scoped_colors[i], label=scoped_labels[i], ax=axs[1])

    # Set x-axis tick locations and labels for both plots
    hour_ticks = np.arange(-20, int(max(data['mcsteps']) / 10000) + 40, 20)
    ticks_loc = axs[0].get_xticks().tolist()
    xlabels = [r"{}".format(int(i)) for i in hour_ticks]
    
    axs[0].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))  # Ensure right tick locations
    axs[1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))  # Ensure right tick locations
    # set the labels on the ticks
    axs[0].set_xticklabels(xlabels, fontsize='medium')
    axs[1].set_xticklabels(xlabels, fontsize='medium')
    # Set y-axis limit for the second plot
    axs[1].set_ylim([0, max(data[cols[-2]]) / 2])

    # Set axis labels and legends
    axs[0].set_ylabel("Cell count", fontsize='medium')
    axs[1].set_ylabel("Cell count", fontsize='medium')
    axs[0].set_xlabel("Hours of simulation", fontsize='medium')
    axs[1].set_xlabel("Hours of simulation", fontsize='medium')
    axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save the figure
    fig.savefig("{}cell_count_data_{}.png".format(outputdir, exp), dpi=300)

    # Show or close the plot based on the 'show' parameter
    if show == "on":
        plt.show()
    else:
        plt.close()

    return



def cyt_concentrations(file, exp, show, outputdir=None):
    # Load data from the CSV file
    data = pd.read_csv(file, sep=",", usecols=["meanconcen","il8mean","il1mean",
                                               "il6mean","il10mean","tnfmean",
                                               "tgfmean","il8std","il1std","il6std",
                                               "il10std","tnfstd","tgfstd"])

    # Cytokine labels and colors for the first plot
    cytokine_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]
    cols = ["il8mean","il1mean", "il6mean","il10mean","tnfmean", "tgfmean"]
    colors = ["blue", "brown", "violet", "red", "yellow", "orange"]

    # Cytokine labels and colors for the second plot
    scoped_cols = ["il8mean","il1mean", "il6mean","il10mean","tnfmean"]
    scoped_colors = ["blue", "brown", "violet", "red", "yellow"]
    scoped_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$"]

    # Create subplots with adjusted layout
    fig, axs = plt.subplots(ncols=1, nrows=2)
    plt.subplots_adjust(left=0.07, bottom=0.075, right=0.900, top=0.96)

    # Plot the data for the first plot
    for i, col in enumerate(cols):
        sns.lineplot(x="meanconcen", y=data[col], data=data, color=colors[i], label=cytokine_labels[i], ax=axs[0])

    # Plot the data for the second plot
    for i, col in enumerate(scoped_cols):
        sns.lineplot(x="meanconcen", y=data[col], data=data, color=scoped_colors[i], label=scoped_labels[i], ax=axs[1])
    # Get x ticks locations   
    ticks_loc = axs[0].get_xticks().tolist()

    # Set x-axis tick locations and labels for both plots
    hour_ticks = np.arange(-20, (int(max(data['meanconcen']) / 10000)) + 40, 20)
    

    # Set major tick fixed locators
    axs[0].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
    axs[1].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))

    # Set x labels for the ticks
    axs[0].set_xticklabels(hour_ticks)
    axs[1].set_xticklabels(hour_ticks)

    # Set y-axis scale for the second plot to logarithmic
    axs[1].set_yscale('log')

    # Set axis labels and legends
    axs[0].set_ylabel("Concentration", fontsize='medium')
    axs[1].set_ylabel("Concentration", fontsize='medium')
    axs[0].set_xlabel("Hours of simulation", fontsize='medium')
    axs[1].set_xlabel("Hours of simulation", fontsize='medium')
    axs[0].legend(loc='center left', bbox_to_anchor=(1, 0.5))
    axs[1].legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # Save the figure
    fig.savefig("{}cytokine_concentration_{}.png".format(outputdir, exp), dpi=300)

    # Show or close the plot based on the 'show' parameter
    if show == "on":
        plt.show()
    else:
        plt.close()

    return


def percent_stacked(df, n_exp, calc, show, outputdir=None):
    # Define colors and cell count labels
    colors = ["blue", "brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    cell_count_labels = ["Endothelial", "Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]

    # Calculate the totals for each experiment
    totals = df[["blueBars", "brownBars", "cyanBars", "violetBars", "redBars", "pinkBars", "yellowBars",
                 "orangeBars", "darkblueBars", "greenBars"]].sum(axis=1)

    # Calculate the percentage for each cell type in each experiment
    percent_data = df[["blueBars", "brownBars", "cyanBars", "violetBars", "redBars", "pinkBars", "yellowBars",
                       "orangeBars", "darkblueBars", "greenBars"]].values / totals.values[:, np.newaxis] * 100

    # Create the stacked bar plot
    fig, ax = plt.subplots()
    barWidth = 0.85
    d = np.arange(0, n_exp)
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]

    # Loop through each cell type and create the stacked bars
    bottom = np.zeros(n_exp)
    for i, label in enumerate(cell_count_labels):
        ax.bar(d, percent_data[:, i], bottom=bottom, color=colors[i], edgecolor='white', width=barWidth, label=label)
        bottom += percent_data[:, i]

    # Set y-axis tick locations and limits
    ax.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(95))
    ax.set_ylim(0, 100)

    # Customize x-axis and labels
    ax.set_xticks(d)
    ax.set_xticklabels(names)
    ax.set_xlabel("Experiment", fontsize='medium')
    ax.set_ylabel(r"% of cells in total", fontsize='medium')

    # Add a legend
    ax.legend(title="Types of cells", fontsize='small', fancybox=True, loc='lower left', 
        bbox_to_anchor=(0., 1.02, 1., .102), ncol=4, mode="expand", borderaxespad=0.)
    # .legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',
    #                   ncols=2, mode="expand", borderaxespad=0.)
    # Save or show the plot based on the 'show' parameter
    plt.savefig(outputdir + calc + "_" + "cell_n_comparison.png", dpi=300)
    if show == "on":
        plt.show()
    else:
        plt.close()

    return



def cytokine_barplot(file, n_exp, exp, show, outputdir):
    all_bars = []
    all_stds = []
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]

    # Calculate the width of the bars and the number of bar groups
    barWidth = 0.8 / (len(names) + 1)
    error_bar_width = 2  # Half the width of the bars for error bars
    num_bars = len(names)
    colors = cm.rainbow(np.linspace(0, 1, num_bars))

    for i in range(n_exp):
        path = "scan_iteration_{}/{}mean_concentration.txt".format(i, file)
        data = pd.read_csv(path, sep=",", usecols=["meanconcen", "il8mean", "il1mean",
                                                   "il6mean", "il10mean", "tnfmean", "tgfmean",
                                                   "il8std", "il1std", "il6std", "il10std", "tnfstd", "tgfstd"])
        IL8m = np.mean(data['il8mean'])
        IL1m = np.mean(data['il1mean'])
        IL6m = np.mean(data['il6mean'])
        IL10m = np.mean(data['il10mean'])
        tnfm = np.mean(data['tnfmean'])
        tgfm = np.mean(data['tgfmean'])

        IL8std = np.mean(data['il8std'])
        IL1std = np.mean(data['il1std'])
        IL6std = np.mean(data['il6std'])
        IL10std = np.mean(data['il10std'])
        tnfstd = np.mean(data['tnfstd'])
        tgfstd = np.mean(data['tgfstd'])

        bars = [IL8m, IL1m, IL6m, IL10m, tnfm, tgfm]
        all_bars.append(bars)
        stds = [IL8std, IL1std, IL6std, IL10std, tnfstd, tgfstd]
        all_stds.append(stds)

        # Create positions for each bar group
        bar_pos = np.arange(len(bars)) + (i * barWidth)

        labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]

        plt.bar(bar_pos, bars, width=barWidth, color=colors[i], edgecolor='black', label=names[i])
        plt.errorbar(bar_pos, bars, yerr=stds, fmt='None', elinewidth=1, capsize=error_bar_width, ecolor='black')

    plt.xticks([r + (barWidth * (num_bars - 1) / 2) for r in range(len(bars))], labels)
    plt.ylabel('Mean Concentration', fontsize='medium')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale("log")
    # Save the figure
    plt.savefig("{}cytokine_data_{}.png".format(outputdir, exp), dpi=300)
    if show == "on":
        plt.show()
    else:
        plt.close()
    return

# day_barplot('combi_code/datafiles/', 7)


def cell_count_barplot(file, n_exp, exp, show, outputdir):
    all_bars = []
    all_stds = []
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]
    cell_count_labels = ["Endothelial", "Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    # Calculate the width of the bars and the number of bar groups
    barWidth = 0.7 / (len(names))
    error_bar_width = barWidth / 2  # Half the width of the bars for error bars
    num_bars = len(names)
    colors = plt.cm.rainbow(np.linspace(0, 1, num_bars))

    for i in range(n_exp):
        path = "scan_iteration_{}/{}cellcount.txt".format(i, file)
        data = pd.read_csv(path, sep=",")
        means = data.mean(axis=0)
        stds = data.std(axis=0)

        bars = [means[str(i)] for i in range(1, 11)]
        all_bars.append(bars)
        all_stds.append(stds)

        # Create positions for each bar group
        bar_pos = np.arange(len(bars)) + (i * barWidth)

        plt.bar(bar_pos, bars, width=barWidth, color=colors[i], edgecolor='black', label=names[i])

    # Extract error values for each experiment separately and plot error bars
    for i in range(n_exp):
        stds_exp = [all_stds[i][str(j)] for j in range(1, 11)]
        bar_pos_exp = np.arange(len(stds_exp)) + (i * barWidth)
        plt.errorbar(bar_pos_exp, all_bars[i], yerr=stds_exp, fmt='None', elinewidth=1, capsize=error_bar_width, ecolor='black')

    plt.xticks([r + (barWidth * (num_bars - 1) / 2) for r in range(len(bars))], cell_count_labels, rotation=30)
    plt.ylabel('Cell Count', fontsize='large')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.yscale('log')
    # Save the figure
    plt.savefig("{}cell_count_data_{}.png".format(outputdir, exp), dpi=300)
    if show == "on":
        plt.show()
    else:
        plt.close()
    return


# cell_count_barplot('combi_code/datafiles/', 7)



def stream_plot_cytokines(file, n_exp, exp, show, outputdir):
    full_data = []
    initial_conc = []
    epx = np.arange(0, n_exp, 1)
    for i in epx:
        path = "scan_iteration_{}/{}mean_concentration.txt".format(i, file)
        data = pd.read_csv(path, sep=",", usecols=["meanconcen", "il8mean", "il1mean",
                                                   "il6mean", "il10mean", "tnfmean", "tgfmean",
                                                   "il8std", "il1std", "il6std", "il10std", "tnfstd", "tgfstd"])
        full_data.append(data)

        # Get the initial concentration value for each cytokine
        initial_conc.append(full_data[i].iloc[0, 1:6])
        # print(initial_conc.head())

    fig, axs = plt.subplots(2, 3, sharex=True, sharey='row', figsize=(14, 7.5))
    plt.subplots_adjust(wspace=0.05, hspace=0.3)  # Adjust the spacing between subplots

    colors = cm.rainbow(np.linspace(0, 1, n_exp))  # Adjust the color map based on the number of experiments

    cytokine_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]

    for i in range(6):  # Loop through cytokines
        axs.flat[i].set_yscale('log')
        # min_ylim = max(initial_conc[i]/100 , 1e-15)
        axs.flat[i].set_ylim([10**(-13), 10**(-7)])
        axs.flat[i].tick_params(axis='x')
        for j in range(n_exp):  # Loop through experiments
            x = full_data[0]['meanconcen']
            y = full_data[j].iloc[:, i + 1]  # Select cytokine data

            # Normalize the data to start from the same initial value
            y_norm = y 

            try:
                if j == 6:
                    axs.flat[i].plot(x, y_norm, color=colors[j], lw=2.4, zorder=5, label="E{}".format(epx[j]+1))
                    # axs.flat[i].scatter(x, y_norm, fc="w", ec=colors[j], s=60, lw=2.4, zorder=5, label="E{} dot".format(epx[j]))
                else:
                    axs.flat[i].plot(x, y_norm, color=colors[j], lw=1.5, label="E{}".format(epx[j]+1))
            except:
                pass

        # Get x ticks locations   
        ticks_loc = axs.flat[i].get_xticks().tolist()

        # Set x-axis tick locations and labels for both plots
        hour_ticks = np.arange(-20, (int(max(data['meanconcen']) / 10000)) + 40, 20)

        # Set major tick fixed locators
        axs.flat[i].xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        
        # Set x labels for the ticks
        axs.flat[i].set_xticklabels(hour_ticks)


        # Add subtitle for each subplot with cytokine name
        axs.flat[i].set_title(cytokine_labels[i])
            
            

    # Set common x and y labels
    fig.text(0.5, 0.04, 'Time(h)', ha='center', fontsize='large')
    fig.text(0.04, 0.5, 'Concentration', va='center', rotation='vertical', fontsize='large')

    plt.legend(loc='center left', bbox_to_anchor=(1, 1.2))
    # Save the figure
    fig.savefig("{}cytokine_data_{}.png".format(outputdir, exp), dpi=300)
    if show == "on":
        plt.show()
    else:
        plt.close()

    return




# stream_plot_cytokines('combi_code/datafiles/', n_exp = 8)



def cell_count_boxplot(file, n_exp, exp, show, outputdir):

    
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]
    colors = ["brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    cell_count_labels = ["Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Create a figure with subplots, adjusting size to fit four plots in one row
    fig, axs = plt.subplots(n_exp, num_intervals, figsize=(16, 4 * n_exp), sharex='col')

    for exp_idx in range(n_exp):
        all_data = pd.DataFrame()  # Store the combined data for each experiment
        for interval_idx, (start_time, end_time) in enumerate(time_intervals):
            path = "scan_iteration_{}/{}cellcount.txt".format(exp_idx, file)
            data = pd.read_csv(path, sep=",")

            # Convert time intervals from hours to steps
            start_steps = start_time * 10000
            end_steps = end_time * 10000

            # Filter data based on time intervals
            interval_data = data.query('mcsteps >= @start_steps and mcsteps <= @end_steps')

            if interval_data.empty:
                # Skip creating the box plot for empty data
                continue

            # Melt the data from wide format to long format
            interval_data = pd.melt(interval_data, id_vars=['mcsteps'], value_vars=cols, var_name='Cell Type', value_name='Cell Count')

            all_data = pd.concat([all_data, interval_data], ignore_index=True)

            # Plot the boxplot and store handles and labels
            boxplot = sns.boxplot(x='Cell Type', y='Cell Count', data=interval_data, ax=axs[exp_idx, interval_idx], palette=colors)
            axs[exp_idx, interval_idx].set_title(f"{names[exp_idx]} - {start_time}-{end_time} hours", fontsize='large')
            axs[exp_idx, interval_idx].set_xlabel('Cell Type', fontsize='large')
            axs[exp_idx, interval_idx].set_ylabel('Cell Count', fontsize='large')

            # Set y-axis locator and formatter
            if interval_idx == 0:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(200))
            else:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(50))

        # Set y-label only for the first subplot in each row
        axs[exp_idx, 0].set_ylabel(f"{names[exp_idx]} - Cell Count")

    # Adjust x-axis ticks and labels for all subplots
    for interval_idx in range(num_intervals):
        axs[n_exp - 1, interval_idx].set_xticks(range(len(cell_count_labels)))
        axs[n_exp - 1, interval_idx].set_xticklabels(cell_count_labels, rotation=45, ha='right')

    # Create a common legend for all the subplots
    fig.legend(names, title="Types of cells", fontsize='small', fancybox=True, loc='lower left', 
        bbox_to_anchor=(0., 1.02, 1., .102), ncol=4, mode="expand", borderaxespad=0.)

    # Adjust spacing between subplots
    plt.tight_layout()
    # Save the figure
    fig.savefig("{}cell_count_data_{}.png".format(outputdir, exp), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()


def cell_count_boxplot(file, n_exp, exp, show, outputdir):
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]
    colors = ["brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    cell_count_labels = ["Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Create a figure with subplots, adjusting size to fit four plots in one row
    fig, axs = plt.subplots(n_exp, num_intervals, figsize=(16, 4 * n_exp), sharex='col')

    # Define a common y-range from 0 to 1200
    common_yrange = (0, 1200)

    for exp_idx in range(n_exp):
        all_data = pd.DataFrame()  # Store the combined data for each experiment
        for interval_idx, (start_time, end_time) in enumerate(time_intervals):
            path = "scan_iteration_{}/{}cellcount.txt".format(exp_idx, file)
            data = pd.read_csv(path, sep=",")

            # Convert time intervals from hours to steps
            start_steps = start_time * 10000
            end_steps = end_time * 10000

            # Filter data based on time intervals
            interval_data = data.query('mcsteps >= @start_steps and mcsteps <= @end_steps')

            if interval_data.empty:
                # Skip creating the box plot for empty data
                continue

            # Melt the data from wide format to long format
            interval_data = pd.melt(interval_data, id_vars=['mcsteps'], value_vars=cols, var_name='Cell Type', value_name='Cell Count')

            all_data = pd.concat([all_data, interval_data], ignore_index=True)

            # Normalize and scale the data to the common y-range
            interval_data['Cell Count'] = np.interp(interval_data['Cell Count'], (interval_data['Cell Count'].min(), interval_data['Cell Count'].max()), common_yrange)

            # Plot the boxplot and store handles and labels
            boxplot = sns.boxplot(x='Cell Type', y='Cell Count', data=interval_data, ax=axs[exp_idx, interval_idx], palette=colors)
            axs[exp_idx, interval_idx].set_title(f"{names[exp_idx]} - {start_time}-{end_time} hours", fontsize='large')
            axs[exp_idx, interval_idx].set_xlabel('Cell Type', fontsize='large')
            axs[exp_idx, interval_idx].set_ylabel('Normalized Cell Count', fontsize='large')

            # Set y-axis locator and formatter
            if interval_idx == 0:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(200))
            else:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(200))

        # Set y-label only for the first subplot in each row
        axs[exp_idx, 0].set_ylabel(f"{names[exp_idx]} - Normalized Cell Count")

    # Adjust x-axis ticks and labels for all subplots
    for interval_idx in range(num_intervals):
        axs[n_exp - 1, interval_idx].set_xticks(range(len(cell_count_labels)))
        axs[n_exp - 1, interval_idx].set_xticklabels(cell_count_labels, rotation=45, ha='right')

    # Create a common legend for all the subplots
    fig.legend(names, title="Types of cells", fontsize='small', fancybox=True, loc='lower left',
        bbox_to_anchor=(0., 1.02, 1., .102), ncol=4, mode="expand", borderaxespad=0.)

    # Adjust spacing between subplots
    plt.tight_layout()
    # Save the figure
    fig.savefig("{}cell_count_data_{}.png".format(outputdir, exp), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()

def cell_count_boxplot(file, n_exp, exp, show, outputdir):
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]
    colors = ["brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    cell_count_labels = ["Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Create a figure with subplots, adjusting size to fit four plots in one row
    fig, axs = plt.subplots(n_exp, num_intervals, figsize=(16, 4 * n_exp), sharex='col')

    # Initialize lists to store boxplot data for normalization
    all_box_data = []

    for exp_idx in range(n_exp):
        all_data = pd.DataFrame()  # Store the combined data for each experiment
        for interval_idx, (start_time, end_time) in enumerate(time_intervals):
            path = "scan_iteration_{}/{}cellcount.txt".format(exp_idx, file)
            data = pd.read_csv(path, sep=",")

            # Convert time intervals from hours to steps
            start_steps = start_time * 10000
            end_steps = end_time * 10000

            # Filter data based on time intervals
            interval_data = data.query('mcsteps >= @start_steps and mcsteps <= @end_steps')

            if interval_data.empty:
                # Skip creating the box plot for empty data
                continue

            # Melt the data from wide format to long format
            interval_data = pd.melt(interval_data, id_vars=['mcsteps'], value_vars=cols, var_name='Cell Type', value_name='Cell Count')

            # Normalize the cell count data using z-score
            interval_data['Cell Count'] = zscore(interval_data['Cell Count'])

            all_data = pd.concat([all_data, interval_data], ignore_index=True)
            all_box_data.append(interval_data['Cell Count'].values)

            # Plot the boxplot and store handles and labels
            boxplot = sns.boxplot(x='Cell Type', y='Cell Count', data=interval_data, ax=axs[exp_idx, interval_idx], palette=colors)
            axs[exp_idx, interval_idx].set_title(f"{names[exp_idx]} - {start_time}-{end_time} hours", fontsize='large')
            axs[exp_idx, interval_idx].set_xlabel('', fontsize='large')
            axs[exp_idx, interval_idx].set_ylabel('', fontsize='large')

            # Set y-axis locator and formatter
            if interval_idx == 0:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(1))  # Adjust as needed
            else:
                axs[exp_idx, interval_idx].yaxis.set_major_locator(MultipleLocator(1))  # Adjust as needed

        # Set y-label only for the first subplot in each row
        axs[exp_idx, 0].set_ylabel(f"{names[exp_idx]} - Z-score normalized cell count")
        
        # Set x-label only for the last row of subplots
        if exp_idx == n_exp - 1:
            axs[exp_idx, interval_idx].set_xlabel("Cell Type")

    # Normalize y-axis limits for all subplots based on the maximum and minimum values
    min_value = -1.5#min([min(box_data) for box_data in all_box_data])
    max_value = 5.1#max([max(box_data) for box_data in all_box_data])

    for exp_idx in range(n_exp):
        for interval_idx in range(num_intervals):
            axs[exp_idx, interval_idx].set_ylim(min_value, max_value)

    # Create a common legend for all the subplots
    fig.legend(names, title="Types of cells", fontsize='small', fancybox=True, loc='lower left', 
        bbox_to_anchor=(0., 1.02, 1., .102), ncol=4, mode="expand", borderaxespad=0.)

    # Adjust spacing between subplots
    plt.tight_layout()
    # Save the figure
    fig.savefig("{}cell_count_data_{}.png".format(outputdir, exp), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()

# output_folder = "Output/"
# cell_count_boxplot('combi_code/datafiles/', n_exp = 8, exp = None, show = "on", outputdir = output_folder)

# file, n_exp, exp, show, outputdir


def cytokine_boxplot(file, n_exp, exp):
    names = ['E{}'.format(i) for i in range(1, n_exp + 1)]
    cytokine_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]
    cols = ["il8mean", "il1mean", "il6mean", "il10mean", "tnfmean", "tgfmean"]
    colors = ["blue", "brown", "violet", "red", "yellow", "orange"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)


    # Create a figure with subplots, adjusting size to fit four plots in one row
    fig, axs = plt.subplots(n_exp, num_intervals, figsize=(16, 4 * n_exp), sharey='row', sharex='col')

    for exp_idx in range(n_exp):
        path = "scan_iteration_{}/{}mean_concentration.txt".format(exp_idx, file)
        # Load data from the file into a DataFrame
        data = pd.read_csv(path, sep=",", usecols=["meanconcen", "il8mean", "il1mean",
                                                   "il6mean", "il10mean", "tnfmean",
                                                   "tgfmean", "il8std", "il1std", "il6std",
                                                   "il10std", "tnfstd", "tgfstd"])

        # Extract the time from the 'meanconcen' column (assuming it's in hours)
        time = data['meanconcen']
        for interval_idx, (start_time, end_time) in enumerate(time_intervals):
            # Convert time intervals from hours to corresponding indices in the DataFrame
            start_idx = int(start_time * 10000)
            end_idx = int(end_time * 10000)

            # Filter data based on time intervals
            interval_data = data.iloc[start_idx:end_idx]

            # Plot the boxplot and store handles and labels
            boxplot = sns.boxplot(data=interval_data[cols], ax=axs[exp_idx, interval_idx], palette=colors)
            axs[exp_idx, interval_idx].set_title(f"{names[exp_idx]} - {start_time}-{end_time} hours",fontsize='large')
            axs[exp_idx, interval_idx].set_xlabel('')
            axs[exp_idx, interval_idx].set_ylabel('Concentration', fontsize='large')

            # Set y-axis to be logarithmic
            axs[exp_idx, interval_idx].set_yscale('log')

        # Set y-label only for the first subplot in each row
        axs[exp_idx, 0].set_ylabel(f"{names[exp_idx]} - Concentration", fontsize='large')

    # Adjust x-axis ticks and labels for all subplots
    for interval_idx in range(num_intervals):
        axs[n_exp - 1, interval_idx].set_xticks(range(len(cytokine_labels)))
        axs[n_exp - 1, interval_idx].set_xticklabels(cytokine_labels, rotation=45, ha='right')

    # Create a common legend for all the subplots
    fig.legend(cytokine_labels, loc='center left', bbox_to_anchor=(1, 0.5))

    # Adjust spacing between subplots
    plt.tight_layout()
    # Save the figure
    fig.savefig("{}cytokine_data_{}.png".format(outputdir, exp), dpi=300)
    if show == "on":
        plt.show()
    else:
        plt.close()

# cytokine_boxplot('combi_code/datafiles/', n_exp = 2)

def cell_count_mean_diff_boxplot(file, n_exp, exp, show, outputdir):
    cell_count_labels = ["Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]
    colors = ["brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Create a figure with subplots, 2 rows and 2 columns
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    for interval_idx, (start_time, end_time) in enumerate(time_intervals):
        mean_diff_data = []  # Initialize a list to hold the mean difference data for box plots
        std_diff_data = []   # Initialize a list to hold the standard deviation data for error bars

        for cell_idx, cell_type in enumerate(cell_count_labels):
            data_for_cell_type_first = np.array([])  # Data for the first experiment
            data_for_cell_type_last = np.array([])   # Data for the last experiment

            for exp_idx in range(n_exp):
                path = "scan_iteration_{}/{}cellcount.txt".format(exp_idx, file)
                data = pd.read_csv(path, sep=",")  # Load the data

                # Convert time intervals from hours to steps
                start_steps = start_time * 10000
                end_steps = end_time * 10000

                # Filter data based on time intervals and cell type
                interval_data = data.query('mcsteps >= @start_steps and mcsteps <= @end_steps')

                if interval_data.empty:
                    # Skip empty data
                    continue

                # Get the cell count data for the current time interval, experiment, and cell type
                cell_count_data = interval_data[cols[cell_idx]]  # Access the correct column

                if exp_idx == 3: # experiment 4 is control
                    data_for_cell_type_first = np.concatenate((data_for_cell_type_first, cell_count_data))
                elif exp_idx == n_exp - 1:
                    data_for_cell_type_last = np.concatenate((data_for_cell_type_last, cell_count_data))

            # Calculate the mean difference between the first and last experiments
            mean_diff = np.mean(data_for_cell_type_last) - np.mean(data_for_cell_type_first)
            mean_diff_data.append(mean_diff)

            # Calculate the standard deviation for the error bars
            std_diff = np.std(data_for_cell_type_last - data_for_cell_type_first)
            std_diff_data.append(std_diff)

        # Calculate subplot position
        row = interval_idx // 2
        col = interval_idx % 2

        # Create the box plot with error bars for the mean difference data
        box = axs[row, col].bar(cell_count_labels, mean_diff_data, color=colors, yerr=std_diff_data, capsize=5) 
        axs[row, col].set_title(f"{start_time}-{end_time} hours", fontsize='large')
        axs[row, col].set_xticks(range(len(cell_count_labels)))  # Set tick positions
        axs[row, col].set_xticklabels(cell_count_labels, rotation=45, ha='right')

        if row == 0:
            axs[row, col].set_xticklabels([])
            axs[row, col].set_xlabel('')

        if row == 1:
            axs[row, col].set_xlabel('Cell Type', fontsize='large')
            axs[row, col].tick_params(axis='x', rotation=45)

        if col == 0:
            axs[row, col].set_ylabel('Mean Difference in Cell Count', fontsize='large')

        if interval_idx == len(time_intervals) - 1:
            axs[row, col].set_xlabel('Cell Type', fontsize='large')
        # else:
        #     axs[row, col].set_xlabel('Cell Type')

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure
    fig.savefig("{}cell_count_mean_difference.png".format(outputdir), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()




def cell_count_mean_diff_violin(file, n_exp, exp, show, outputdir):
    cell_count_labels = ["Resting Neutrophils", "Monocytes", "Fibroblasts",
                         "Activated Neutrophils", "Necrotic Neutrophils", "Resting Monocytes",
                         "Macrophages I", "Macrophages II", "Myofibroblasts"]
    cols = ["2", "3", "4", "5", "6", "7", "8", "9", "10"]
    colors = ["brown", "cyan", "violet", "red", "pink", "yellow", "orange", "darkblue", "green"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    custom_colors = colors[:len(cell_count_labels)]

    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    ymin, ymax = float('inf'), float('-inf')

    for interval_idx, (start_time, end_time) in enumerate(time_intervals):
        data_for_interval = []  

        for cell_idx, cell_type in enumerate(cell_count_labels):
            data_for_cell_type = []  

            for exp_idx in range(n_exp):
                path = f"scan_iteration_{exp_idx}/{file}cellcount.txt"
                data = pd.read_csv(path, sep=",")

                start_steps = start_time * 10000
                end_steps = end_time * 10000

                interval_data = data.query('mcsteps >= @start_steps and mcsteps <= @end_steps')

                if interval_data.empty:
                    continue

                cell_count_data = interval_data[cols[cell_idx]]

                if exp_idx == 7:
                    if exp_idx == 3:
                        if not data_for_cell_type:
                            data_for_cell_type = [0] * len(cell_count_data)
                        diff_data = np.array(cell_count_data) - np.array(data_for_cell_type)
                        data_for_cell_type = list(diff_data)
                    else:
                        data_for_cell_type += list(cell_count_data)
                elif exp_idx == 3:
                    data_for_cell_type += list(cell_count_data)

            z_scores = zscore(data_for_cell_type)
            data_for_interval.append(list(z_scores))

        row = interval_idx // 2
        col = interval_idx % 2

        positions = np.arange(1, len(cell_count_labels) + 1)
        parts = axs[row, col].violinplot(data_for_interval, positions=positions, showmeans=True, showmedians=False, vert=True)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(custom_colors[i])
            pc.set_edgecolor('black')
            pc.set_alpha(1)
            pc.set_linewidth(2)

        parts['cmeans'].set_color('black')
        parts['cmeans'].set_linewidth(2)

        axs[row, col].set_title(f"{start_time}-{end_time} hours", fontsize='large')
        axs[row, col].set_xticks(positions)
        axs[row, col].set_xticklabels([])

        if col == 0:
            axs[row, col].set_ylabel('Z-scored normalized cell count difference', fontsize='large')
        else:
            axs[row, col].set_ylabel('')

        if row == 1:
            axs[row, col].set_xlabel('Cell Type', fontsize='large')
            axs[row, col].set_xticklabels(cell_count_labels, rotation=45, ha='right')
        else:
            axs[row, col].set_xlabel('')

        for partname in ('cbars', 'cmins', 'cmaxes'):
            vp = parts[partname]
            vp.set_color('black')
            vp.set_linewidth(2)

    for ax in axs.flat:
        ax.set_ylim(-3, 3)  # Adjust the y-limits as needed

    plt.tight_layout()
    fig.savefig(f"{outputdir}cell_count_mean_difference_violin.png", dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()

cell_count_mean_diff_violin('combi_code/datafiles/', 8, exp = "diff_boxplot", show="on", outputdir = "Output/")


def cytokine_mean_diff_boxplot(file, n_exp, exp, outputdir, show):
    cytokine_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]
    cols = ["il8mean", "il1mean", "il6mean", "il10mean", "tnfmean", "tgfmean"]
    colors = ["blue", "brown", "violet", "red", "yellow", "orange"]
    # cytokine_labels = ["IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]
    # cols = ["il1mean", "il6mean", "il10mean", "tnfmean", "tgfmean"]
    # colors = ["brown", "violet", "red", "yellow", "orange"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Initialize variables to store overall min and max values
    overall_min = float('inf')
    overall_max = float('-inf')

    # Create a figure with subplots, 2 rows and 2 columns
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    for interval_idx, (start_time, end_time) in enumerate(time_intervals):
        mean_diff_data = []  # Initialize a list to hold the mean difference data for box plots
        std_diff_data = []   # Initialize a list to hold the standard deviation data for error bars

        for cytokine_idx, cytokine_type in enumerate(cytokine_labels):
            data_for_cytokine_first = np.array([])  # Data for the first experiment
            data_for_cytokine_last = np.array([])   # Data for the last experiment

            for exp_idx in range(n_exp):
                path = "scan_iteration_{}/{}mean_concentration.txt".format(exp_idx, file)
                data = pd.read_csv(path, sep=",", usecols=["meanconcen"] + cols)  # Load the data

                # Convert time intervals from hours to steps
                start_steps = start_time * 10000 
                end_steps = end_time * 10000 

                # Filter data based on time intervals and cytokine type
                interval_data = data.query('meanconcen >= @start_steps and meanconcen <= @end_steps')

                if interval_data.empty:
                    # Skip empty data
                    continue

                # Get the cytokine concentration data for the current time interval, experiment, and cytokine type
                cytokine_concentration_data = interval_data[cols[cytokine_idx]]  # Access the correct column

                if exp_idx == 3: # experiment 4 is control
                    data_for_cytokine_first = np.concatenate((data_for_cytokine_first, cytokine_concentration_data))
                elif exp_idx == n_exp - 1:
                    data_for_cytokine_last = np.concatenate((data_for_cytokine_last, cytokine_concentration_data))

            # Calculate the mean difference between the first and last experiments
            mean_diff = np.mean(data_for_cytokine_last) - np.mean(data_for_cytokine_first)
            mean_diff_data.append(mean_diff)

            # Calculate the standard deviation for the error bars
            std_diff = np.std(data_for_cytokine_last - data_for_cytokine_first)
            std_diff_data.append(std_diff)


        # Calculate subplot position
        row = interval_idx // 2
        col = interval_idx % 2

        # Create the box plot with error bars for the mean difference data
        box = axs[row, col].bar(cytokine_labels, mean_diff_data, color=colors, yerr=std_diff_data, capsize=5)
        axs[row, col].set_title(f"{start_time}-{end_time} hours", fontsize='large')
        axs[row, col].set_xticks(range(len(cytokine_labels)))  # Set tick positions
        axs[row, col].set_xticklabels(cytokine_labels, rotation=45, ha='right')

        if row == 0:
            axs[row, col].set_xticklabels([])
            axs[row, col].set_xlabel('')

        if row == 1:
            axs[row, col].set_xlabel('Cytokine Type', fontsize='large')
            axs[row, col].tick_params(axis='x', rotation=45)

        if col == 0:
            axs[row, col].set_ylabel('Mean Difference in Concentration', fontsize='large')

        if interval_idx == len(time_intervals) - 1:
            axs[row, col].set_ylabel('')
        # else:
            # axs[row, col].set_xlabel('Cytokine Type')

        # Set y-axis to logarithmic scale
        # axs[row, col].set_yscale('log')

        # Set y-axis limits to be the same for all subplots
        axs[row, col].set_ylim(-1.6*10**-10, 8.0*10**-9)

    # Adjust spacing between subplots
    plt.tight_layout()

    # Save the figure
    fig.savefig("{}cytokine_mean_difference.png".format(outputdir), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()


def cytokine_mean_diff_violin(file, n_exp, exp, outputdir, show):
    cytokine_labels = ["IL-8", "IL-1", "IL-6", "IL-10", r"TNF-$\alpha$", r"TGF-$\beta$1"]
    cols = ["il8mean", "il1mean", "il6mean", "il10mean", "tnfmean", "tgfmean"]
    colors = ["blue", "brown", "violet", "red", "yellow", "orange"]
    time_intervals = [(0, 24), (24, 48), (48, 72), (72, 96)]
    num_intervals = len(time_intervals)

    # Define a custom color scheme using the colors you specified
    custom_colors = colors[:len(cytokine_labels)]

    fig, axs = plt.subplots(2, 2, figsize=(12, 12))

    # Initialize variables for y-limits
    ymin, ymax = float('inf'), float('-inf')

    for interval_idx, (start_time, end_time) in enumerate(time_intervals):
        data_for_interval = []  

        for cytokine_idx, cytokine_type in enumerate(cytokine_labels):
            data_for_cytokine_8 = []  # Data for experiment 8
            data_for_cytokine_4 = []  # Data for experiment 4

            for exp_idx in range(n_exp):
                path = "scan_iteration_{}/{}mean_concentration.txt".format(exp_idx, file)
                data = pd.read_csv(path, sep=",", usecols=["meanconcen"] + cols)

                start_steps = start_time * 10000
                end_steps = end_time * 10000

                interval_data = data.query('meanconcen >= @start_steps and meanconcen <= @end_steps')

                if interval_data.empty:
                    continue

                cytokine_concentration_data = interval_data[cols[cytokine_idx]]

                if exp_idx == 7:  # Experiment 8
                    data_for_cytokine_8 += list(cytokine_concentration_data)
                elif exp_idx == 3:  # Experiment 4
                    data_for_cytokine_4 += list(cytokine_concentration_data)

            # Calculate the mean difference by subtracting the two lists
            data_for_cytokine_diff = [a - b for a, b in zip(data_for_cytokine_8, data_for_cytokine_4)]

            # Z-score normalize the difference data
            data_for_cytokine_zscore = (data_for_cytokine_diff - np.mean(data_for_cytokine_diff)) / np.std(data_for_cytokine_diff)
            
            data_for_interval.append(data_for_cytokine_zscore)

        row = interval_idx // 2
        col = interval_idx % 2

        # Calculate positions for violins
        positions = np.arange(1, len(cytokine_labels) + 1)

        # Plot individual violins with specified colors and outlines
        parts = axs[row, col].violinplot(data_for_interval, positions=positions, showmeans=True, showmedians=False, vert=True)

        for i, pc in enumerate(parts['bodies']):
            pc.set_facecolor(custom_colors[i])
            pc.set_edgecolor('black')  # Set the outline to black
            pc.set_alpha(1)
            pc.set_linewidth(2)  # Set the linewidth to 2

        # Set the color of the mean bar to black and adjust its properties
        parts['cmeans'].set_color('black')
        parts['cmeans'].set_linewidth(2)

        axs[row, col].set_title(f"{start_time}-{end_time} hours", fontsize='large')
        axs[row, col].set_xticks(positions)
        axs[row, col].set_xticklabels(cytokine_labels, rotation=45, ha='right')

        if row == 0:
            axs[row, col].set_xticklabels([])
            axs[row, col].set_xlabel('')

        if row == 1:
            axs[row, col].set_xlabel('Cytokine Type', fontsize='large')
            axs[row, col].tick_params(axis='x', rotation=45)

        if col == 1:
            axs[row, col].set_ylabel('')  # Remove y label for the second column subfigures
        else:
            axs[row, col].set_ylabel('Z-Score Normalized Concentration Difference', fontsize='large')

        # Update y-limits
        ymin = min(ymin, np.min(data_for_interval))
        ymax = max(ymax, np.max(data_for_interval))

        # Set the color of the error bars and the line connecting error bar tips to black
        for partname in ('cbars', 'cmins', 'cmaxes'):
            vp = parts[partname]
            vp.set_color('black')
            vp.set_linewidth(2)

    # Set y-limits for all subplots
    for ax in axs.flat:
        ax.set_ylim(-3, 3)

    plt.tight_layout()
    fig.savefig("{}cytokine_mean_difference_violin.png".format(outputdir), dpi=300)

    if show == "on":
        plt.show()
    else:
        plt.close()

# cytokine_mean_diff_violin('combi_code/datafiles/', 8, exp = "diff_boxplot", outputdir = "Output/", show="on", )


def run_experiments(file, experiment, show, n_exp, calc=None):
    """
    experiment types:
    cell count
    cytokine concentration
    bar plot
    none = all
    """
    output_folder = "Output/"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    full_data = []  # Initialize full_data as an empty list

    if experiment == "cell count":
        for i in range(n_exp):
            cell_count("scan_iteration_{}/{}cellcount.txt".format(i, file), "E{}".format(i + 1),
                       show, output_folder)
            print("Experiment {} for cell count rendered!".format(i + 1))
    elif experiment == "cytokine concentration":
        for i in range(n_exp):
            cyt_concentrations("scan_iteration_{}/{}mean_concentration.txt".format(i, file), "E{}".format(i + 1),
                               show, output_folder)
            print("Experiment {} for cytokine concentration rendered!".format(i + 1))
    elif experiment == "cytokine bar plot":
        cytokine_barplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        print("Experiment bar cytokine plot rendered!")
    elif experiment == "cell count bar plot":
        cell_count_barplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        print("Experiment cell count bar plot rendered!")
    elif experiment == "stream plot cytokine":
        stream_plot_cytokines(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        print("Experiment stream plot cytokine rendered!")
    elif experiment == "cell count boxplot":
        cell_count_boxplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        print("Experiment cell count boxplot rendered!")

    elif experiment == "bar plot":
        for i in range(n_exp):
            data = pd.read_csv("scan_iteration_{}/{}cellcount.txt".format(i, file),
                               sep=",", usecols=["mcsteps", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])
            full_data.append(data)


        if calc == "max":
            raw_data = {'blueBars': [full_data[i]['1'].max() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].max() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].max() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].max() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].max() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].max() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].max() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].max() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].max() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].max() for i in range(n_exp)]}
        elif calc == "min":
            raw_data = {'blueBars': [full_data[i]['1'].min() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].min() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].min() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].min() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].min() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].min() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].min() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].min() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].min() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].min() for i in range(n_exp)]}
        else:
            raw_data = {'blueBars': [full_data[i]['1'].mean() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].mean() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].mean() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].mean() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].mean() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].mean() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].mean() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].mean() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].mean() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].mean() for i in range(n_exp)]}

        df = pd.DataFrame(raw_data)

        percent_stacked(df, n_exp, calc, show, outputdir=output_folder)
        print("Experiment bar graph rendered!")

    elif experiment is None:
        show = "Off"
        for i in range(n_exp):
            cell_count("scan_iteration_{}/{}cellcount.txt".format(i, file), "E{}".format(i + 1), show, output_folder)
            cyt_concentrations("scan_iteration_{}/{}mean_concentration.txt".format(i, file), "E{}".format(i + 1), show,
                               output_folder)
            print("Experiment {} rendered!".format(i + 1))

            # Data
            data = pd.read_csv("scan_iteration_{}/{}cellcount.txt".format(i, file),
                               sep=",", usecols=["mcsteps", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"])
            full_data.append(data)

        if calc == "max":
            raw_data = {'blueBars': [full_data[i]['1'].max() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].max() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].max() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].max() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].max() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].max() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].max() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].max() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].max() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].max() for i in range(n_exp)]}
        elif calc == "min":
            raw_data = {'blueBars': [full_data[i]['1'].min() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].min() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].min() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].min() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].min() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].min() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].min() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].min() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].min() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].min() for i in range(n_exp)]}
        else:
            raw_data = {'blueBars': [full_data[i]['1'].mean() for i in range(n_exp)],
                        'brownBars': [full_data[i]['2'].mean() for i in range(n_exp)],
                        'cyanBars': [full_data[i]['3'].mean() for i in range(n_exp)],
                        'violetBars': [full_data[i]['4'].mean() for i in range(n_exp)],
                        'redBars': [full_data[i]['5'].mean() for i in range(n_exp)],
                        'pinkBars': [full_data[i]['6'].mean() for i in range(n_exp)],
                        'yellowBars': [full_data[i]['7'].mean() for i in range(n_exp)],
                        'orangeBars': [full_data[i]['8'].mean() for i in range(n_exp)],
                        'darkblueBars': [full_data[i]['9'].mean() for i in range(n_exp)],
                        'greenBars': [full_data[i]['10'].mean() for i in range(n_exp)]}

        df = pd.DataFrame(raw_data)

        percent_stacked(df, n_exp, calc, show, outputdir=output_folder)
        experiment = "cytokine bar plot"
        cytokine_barplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        experiment = "cell count bar plot"
        cell_count_barplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        experiment = "stream plot cytokine"
        stream_plot_cytokines(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        experiment = "cell count boxplot"
        cell_count_boxplot(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        experiment = "cell count mean diff boxplot"
        cell_count_mean_diff_violin(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        experiment = "cytokine mean diff boxplot"
        cytokine_mean_diff_violin(file=file, n_exp=n_exp, exp=experiment, show= show, outputdir=output_folder)
        print("All Experiments rendered!")
       
    
    return

# run_experiments(file = 'combi_code/datafiles/', experiment= None, show = "off", n_exp = 8, calc ="mean")

