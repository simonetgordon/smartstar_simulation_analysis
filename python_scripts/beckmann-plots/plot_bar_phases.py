import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import matplotlib.cm as cm
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import *

##########################################################################################################
#                                     Plot mass growth of BHs                                            #
#
# to run: python plot_bar_phases.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
# python -i plot_bar_phases.py data_files/data-1B.m16-4dx.csv bar-phases-1B.m16_spectral
# python -i plot_bar_phases.py data_files/data-2B.RSb08.csv bar-phases-2B.b08
# list data files in order of low res -> high res
##########################################################################################################


# Function to set up the plot environment
def setup_plot_env(fontsize, linewidth):
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth


# Function to format and create the subplots
def create_subplots(num_subplots, xlim, ylim, time_cutoff, fontsize, title):
    fig, ax = plt.subplots(num_subplots, 1, sharex=True)

    i = 0
    # for i in range(num_subplots):
    ax.set_xticks(np.arange(0.1, time_cutoff+0.1, 0.1))
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.02))
    ax.tick_params(axis="x", which='minor', length=2, direction="in")
    ax.tick_params(axis="x", which='major', labelsize=fontsize-1, width=1.5, length=3, direction="in")
    ax.tick_params(axis="y", which='major', labelsize=fontsize-1)
    ax.tick_params(axis="y", which='minor', labelsize=fontsize-2)
    ax.set_yscale('log')
    ax.set_xlim([0, xlim+0.01]) # for truncated view

    ax.set_title(title)
    #ax.set_ylim([10.5, ylim_mass+0.01]) # for truncated view 270msun: [240, ylim_mass+0.01]
    #ax.set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    ax.set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    ax.set_ylim(ylim) # [2e-9, 8e-4] for 270msun, [5e-9, 6e-2] for 10.8msun-no-SN, [2e-9, 8e-4] for 10.8msun
    ax.set_xlabel('Black Hole Age (Myr)', fontdict=None)

    return fig, ax

def plot_extra_line_mass_growth(j, data_file, label_extra, alpha=0.8):
        # Load the CSV file into a DataFrame
        df = pd.read_csv(data_file)

        # Extract the columns you're interested in
        age = df['age'].values/1e6
        bh_mass = df['BH mass'].values
        accrate = df['accrate'].values

        # 1) BH Mass
        axs[0].plot(age, bh_mass, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)

        # 2) Accretion Rates
        axs[1].plot(age, accrate, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)
        axs[1].plot(age, eddington_rate(bh_mass), color=c[j], linestyle='dashed', label=label_extra, alpha=alpha)

        return 0

def moving_min_max(data, window_size):
    """Calculate the moving min and max values."""
    min_vals = np.array([np.min(data[i:i+window_size]) for i in range(len(data) - window_size + 1)])
    max_vals = np.array([np.max(data[i:i+window_size]) for i in range(len(data) - window_size + 1)])
    return min_vals, max_vals

def moving_average(data, window_size):
    """Calculate the moving average for smoothing."""
    return np.convolve(data, np.ones(window_size) / window_size, mode='valid')

if __name__ == "__main__":
    # Set up plot parameters
    j = 0
    alpha = 0.8
    xlim = 1
    ylim_mass = 110 # 4000 for 270msun, 200 for 10.8msun baseline, 1200 for no-sn, 110 for 10.8msun
    time_cutoff = 1
    smooth_simulations = 1 # number of simulations to smooth
    window = 4 # window size for smoothing
    extra_line = False # true for 10.8msun, false for 270msun
    n = 4 # half number of simulations to plot

    # Text format
    linewidth = 1.5
    fontsize = 14

    # Set up plot environment
    setup_plot_env(fontsize, linewidth)

    # Set up figure and subplots
    num_subplots = 1
    ylim = [8e-5, 2e-1]
    l = tidy_data_labels(bhl_object_labels)
    title = r"Bar Phases in {}".format(l[0]) 
    fig, ax = create_subplots(num_subplots, xlim, ylim, time_cutoff, fontsize, title)

    # Line colours
    c_s1 = extract_colors('viridis', n, portion="middle", start=0.35, end=0.92) # start=0.33, end=0.92 for 10.8msun-no-sn
    c_s2 = extract_colors('magma', n, portion="middle", start=0.3, end=0.85) # start=0.3, end=0.85 for 10.8msun-no-sn
    c = np.concatenate((c_s1, c_s2))[1]

    # Set BHL properties parameters and resample data    
    times = [BHL.ages/1e6 for BHL in bhl_object_list]
    min_time = min([min(t) for t in times])
    min_final_time = min([t[-1] for t in times])
    num_steps = int(len(times[-1])/5)
    common_time = np.linspace(min_time, min_final_time, num_steps)

    # Resample mass and accretion rates data
    accrates_og = [BHL.accrates for BHL in bhl_object_list][0]
    # mass = resample_data([BHL.mass for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    # accretion_rates = resample_data([BHL.accrates for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)

    window = 200
    accrate_smooth = moving_average(accrates_og, window)
    mass_smooth = moving_average([BHL.mass for BHL in bhl_object_list][0], window)
    common_time = np.linspace(min_time, min_final_time, len(accrate_smooth))
    accrate_min, accrate_max = moving_min_max(accrates_og, 1)
    min_max_time = np.linspace(min_time, min_final_time, len(accrate_min))

    # Plot BH Mass and Accretion Rates
    #axs[0].plot(common_time, mass[i], color=c, linestyle='solid', label=l[i], alpha=alpha)
    i = 0
    ax.plot(common_time, accrate_smooth, color='royalblue', linestyle='solid', label=l[i], alpha=1)
    ax.fill_between(min_max_time, accrate_min, accrate_max, color='cornflowerblue', alpha=0.3) # min-max region
    #ax.plot(common_time, eddington_rate(mass_smooth), color=c, linestyle='dashed', label=l[i], alpha=alpha)

    # Plot phases
    y_bins = np.logspace(np.log10(ylim[0]), np.log10(ylim[1]), 6)
    
    # Extract N evenly spaced colors
    cmap = plt.get_cmap('Pastel1')
    N = 5
    colors = cmap(np.linspace(0, 0.5, N))
    buffer = 1e-5

    # 4e9+ disc forming 0.18-0.21 for 1B.m16-4dx
    ax.fill_between([0.68, 1.1], y_bins[0], y_bins[1]-buffer, color=cmap([0.8]), alpha=0.8, label="Disc") # 'lightgreen'
    #ax.fill_between([0.84, 1.1], y_bins[0], y_bins[1]-buffer, color=cmap([0.8]), alpha=0.8)
    i += 1
    # spiral arms 0.33 - 0.38 for 1B.m16-4dx, 
    ax.fill_between([0.86, 1.1], y_bins[1], y_bins[2]-buffer*2, color=cmap([0.4]), alpha=0.8, label="Spiral Arms") # color='cornflowerblue'
    #ax.fill_between([0.93, 1.1], y_bins[1], y_bins[2]-buffer*2, color=cmap([0.4]), alpha=0.8)
    i += 1
    # bar 0.44-0.59 for 1B.m16-4dx
    ax.fill_between([0.29, 0.51], y_bins[2], y_bins[3]-buffer*20, color=colors[i], alpha=0.8, label="Bar") # color='purple'
    ax.fill_between([0.64, 0.87], y_bins[2], y_bins[3]-buffer*20, color=colors[i], alpha=0.8)
    ax.fill_between([0.9, 0.97], y_bins[2], y_bins[3]-buffer*20, color=colors[i], alpha=0.8)
    ax.fill_between([0.99, 1.1], y_bins[2], y_bins[3]-buffer*20, color=colors[i], alpha=0.8)
    i += 2
    # fragmentation 0.54-1 for 1B.m16-4dx, 0.95-1.1 for 2B.b08
    ax.fill_between([0.95, 1.1], y_bins[3], y_bins[4]-buffer*100, color=colors[i], alpha=0.8, label="Fragmentation") # color='gold'

    # Ring 0.94-0.99
    ax.fill_between([0.75, 0.78], y_bins[4], y_bins[5]-buffer*100, color=cmap([0.0]), alpha=0.8, label="Ring") # color='orange'
    ax.fill_between([0.94, 0.96], y_bins[4], y_bins[5]-buffer*100, color=cmap([0.0]), alpha=0.8)



    # Legend
    accrate_line = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth)]
    ax.legend(fontsize=fontsize-4, ncol=2, loc="upper left", handlelength=0.7) # "lower right" for no-sn
    #axs[1].legend(accrate_line, [r"$\rm \dot{M}_{Edd}$"], loc="lower right", fontsize=fontsize-2.2, ncol=1)

    # Save plot
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(5.5, 4.0)
    #plot_name = 'mass_growth-1S+2S-no-SN' + '.pdf'
    plot_name = sys.argv[-1] + '.pdf'
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + plot_name)
