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
# python -i plot_bar_phases.py data_files/data-1B.m16-4dx.csv bar-phases-1B.m16
# python -i plot_bar_phases.py data_files/data-2B.RSb08.csv bar-phases-2B.b08
# python -i plot_bar_phases.py data_files/data-1S.m04-no-SN.csv bar-phases-1S.m04-no-SN
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
    ylim = [2e-4, 1.2e-1]
    l = tidy_data_labels(bhl_object_labels)
    title = r"Instability Phases in {}".format(l[0]) 
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
    i = 0
    ax.plot(common_time, accrate_smooth, color='royalblue', linestyle='solid', label=l[0], alpha=1)

    # Plot phases
    y_bins = np.logspace(np.log10(ylim[0]), np.log10(ylim[1]), 7)
    
    # Extract N evenly spaced colors
    cmap = plt.get_cmap('Pastel1')
    N = 6
    colors = cmap(np.linspace(0, 0.5, N))
    buffer = 1e-5

    j = 1
    #4e9+ disc forming
    ax.fill_between([0.18, 0.52], y_bins[j], y_bins[j+1]-buffer*5, color=cmap([0.8]), alpha=0.8, label="Disc") # 
    ax.fill_between([0.68, 1.1], y_bins[j], y_bins[j+1]-buffer*5, color=cmap([0.8]), alpha=0.8)
    i += 1
    # spiral arms
    # ax.fill_between([0.21, 0.6], y_bins[j+1], y_bins[j+2]-buffer*10, color=cmap([0.4]), alpha=0.8, label="Spiral Arms") # color='cornflowerblue'
    # ax.fill_between([0.86, 1.1], y_bins[j+1], y_bins[j+2]-buffer*10, color=cmap([0.4]), alpha=0.8)
    ax.fill_between([0.32, 0.57], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8, label="Spiral Arms")
    ax.fill_between([0.61, 0.63], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8) # dual bars amongst higher modes
    ax.fill_between([0.67, 0.7], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8) # ring
    ax.fill_between([0.71, 0.8], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8) # m2 mode elevation moving outwards
    ax.fill_between([0.84, 0.87], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8) # strong ring-like feature as elevated m2 mode 0,03 pc
    ax.fill_between([0.94, 1.02], y_bins[j+1], y_bins[j+2]-buffer*14, color=cmap([0.4]), alpha=0.8) # 0.55 - 0.68 pc ring like feature showing in elevated m2 over m1 for a growing portion of the disk as multiple rings move out.
    i += 2
    # bar 
    # ax.fill_between([0.29, 0.51], y_bins[j+2], y_bins[j+3]-buffer*50, color=colors[i], alpha=0.8)
    # ax.fill_between([0.64, 0.87], y_bins[j+2], y_bins[j+3]-buffer*50, color=colors[i], alpha=0.8) # color='purple'
    # ax.fill_between([0.9, 0.97], y_bins[j+2], y_bins[j+3]-buffer*50, color=colors[i], alpha=0.8)
    # ax.fill_between([0.99, 1.1], y_bins[j+2], y_bins[j+3]-buffer*50, color=colors[i], alpha=0.8, label="Bar")
    # 0.88-94 strength, 0.018 pc - the bar radius isn't consistently at 0.018 pc, but there is a more or less constant peak in m2 at 0.018 pc with radius_bar being marked at or beyond r > radius_bar. 
    # Also verified visually by the focusing of the density along a plane of consistent orientation. 
    ax.fill_between([0.42, 0.55], y_bins[j+2], y_bins[j+3]-buffer*40, color=colors[i], alpha=0.8, label="Bar")
    ax.fill_between([0.66, 0.69], y_bins[j+2], y_bins[j+3]-buffer*40, color=colors[i], alpha=0.8) # 0.90 strength, 0.01 pc
    ax.fill_between([0.77, 0.8], y_bins[j+2], y_bins[j+3]-buffer*40, color=colors[i], alpha=0.8) # 0.90 strength, 0.01 pc
    i += 1
    i += 1
    # clumps
    #ax.fill_between([0.95, 1.1], y_bins[j+3], y_bins[j+4]-buffer*100, color='gold', alpha=0.2, label="Fragmentation") # color='gold'
    ax.fill_between([0.53, 1.1], y_bins[j+3], y_bins[j+4]-buffer*10, color='gold', alpha=0.2, label="Fragmentation") # color='gold'

    # Plot line
    ax.fill_between(min_max_time, accrate_min, accrate_max, color='cornflowerblue', alpha=0.3) # min-max region

    # Legend
    accrate_line = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth)]
    ax.legend(fontsize=fontsize-2, ncol=1, loc="upper left", handlelength=0.7) # "lower right" for no-sn

    # Save plot
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(5.5, 4.0)
    #plot_name = 'mass_growth-1S+2S-no-SN' + '.pdf'
    plot_name = sys.argv[-1] + '.pdf'
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + plot_name)
