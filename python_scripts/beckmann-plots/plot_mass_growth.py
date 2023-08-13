import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
import matplotlib.cm as cm
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import *

##########################################################################################################
#                                     Plot mass growth of BHs                                            #
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
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
def create_subplots(num_subplots, xlim, time_cutoff, fontsize, title):
    fig, axs = plt.subplots(num_subplots, 1, sharex=True)

    for i in range(num_subplots):
        axs[i].set_xticks(np.arange(0.1, time_cutoff+0.1, 0.1))
        axs[i].minorticks_on()
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(0.02))
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize-1, width=1.5, length=3, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize-1)
        axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)
        axs[i].set_yscale('log')
        axs[i].set_xlim([0, xlim+0.01]) # for truncated view

    axs[0].set_title(title)
    axs[0].set_ylim([240, ylim_mass+0.01]) # for truncated view
    axs[0].set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[-1].set_xlabel('Black Hole Age (Myr)', fontdict=None)

    return fig, axs


if __name__ == "__main__":
    # Set up plot parameters
    j = 0
    title = "Baseline Growth"
    alpha = 1
    xlim = 1
    ylim_mass = 4000
    time_cutoff = 1
    smooth_simulations = 4 # number of simulations to smooth
    window = 20

    # Text format
    linewidth = 1.5
    fontsize = 14

    # Set up plot environment
    setup_plot_env(fontsize, linewidth)

    # Set up figure and subplots
    num_subplots = 2
    fig, axs = create_subplots(num_subplots, xlim, time_cutoff, fontsize, title)

    # Line colours
    n = 2
    c_s1 = extract_colors('viridis', n, portion="middle")
    c_s2 = extract_colors('magma', n, portion="middle")
    c = np.concatenate((c_s1, c_s2))

    # Set BHL properties parameters and resample data
    l = tidy_data_labels(bhl_object_labels)
    times = [BHL.ages/1e6 for BHL in bhl_object_list]
    min_time = min([min(t) for t in times])
    min_final_time = min([t[-1] for t in times])
    num_steps = int(len(times[-1])/5)
    common_time = np.linspace(min_time, min_final_time, num_steps)

    # Resample mass and accretion rates data
    mass = resample_data([BHL.mass for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    accretion_rates = resample_data([BHL.accrates for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)

    # Plot BH Mass and Accretion Rates
    for i in range(len(mass)):
        axs[0].plot(common_time, mass[i], color=c[i], linestyle='solid', label=l[i], alpha=alpha)
        axs[1].plot(common_time, accretion_rates[i], color=c[i], linestyle='solid', label=l[i], alpha=alpha)
        axs[1].plot(common_time, eddington_rate(mass[i]), color=c[i], linestyle='dashed', label=l[i], alpha=alpha)
        j += 1

    accrate_line = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth)]

    # Include legends and save the plot
    axs[0].legend(fontsize=fontsize-2.2, ncol=2, loc="lower right")
    axs[1].legend(accrate_line, [r"$\rm \dot{M}_{Edd}$"], loc="lower right", fontsize=fontsize-2.2, ncol=1)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.7, 4.7)
    plot_name = 'mass_growth-baselines_resolution' + '.pdf'
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + plot_name)
