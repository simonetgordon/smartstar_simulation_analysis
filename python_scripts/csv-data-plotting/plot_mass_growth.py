import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import *

##########################################################################################################
#                                     Plot mass growth of BHs                                            #
#
# to run: python plot_mass_growth.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
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
    axs[0].set_ylim([10.5, ylim_mass+0.01]) # for truncated view 270msun: [240, ylim_mass+0.01]
    axs[0].set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[1].set_ylim([2e-9, 8e-4]) # [2e-9, 8e-4] for 270msun, [5e-9, 6e-2] for 10.8msun-no-SN, [2e-9, 8e-4] for 10.8msun
    axs[-1].set_xlabel('Black Hole Age (Myr)', fontdict=None)

    return fig, axs


def plot_extra_line_mass_growth(j, data_file, label_extra, alpha=0.8, n=10):
        # Load the CSV file into a DataFrame
        df = pd.read_csv(data_file)

        # Extract the columns you're interested in
        age = df['age'].values/1e6
        bh_mass = df['BH mass'].values
        accrate = df['accrate'].values

        # try adaptively smoothing the data
        age_av = adaptive_moving_average(age, window_size=n)
        bh_mass_av = adaptive_moving_average(bh_mass, window_size=n)
        accrate_av = adaptive_moving_average(accrate, window_size=n)

        # 1) BH Mass
        axs[0].plot(age_av, bh_mass_av, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)

        # 2) Accretion Rates
        axs[1].plot(age_av, accrate_av, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)
        axs[1].plot(age_av, eddington_rate(bh_mass_av), color=c[j], linestyle='dashed', label=label_extra, alpha=alpha)

        return 0

if __name__ == "__main__":
    # Set up plot parameters
    j = 0
    title = r"$\rm 10.8 \, M_\odot$ Baseline Growth" #"No-SN Growth"
    alpha = 0.8
    xlim = 1
    ylim_mass = 110 # 4000 for 270msun, 200 for 10.8msun baseline, 1200 for no-sn, 110 for 10.8msun
    time_cutoff = 1
    smooth_simulations = 8 # number of simulations to smooth
    window = 4 # window size for smoothing
    extra_line = True # true for 10.8msun, false for 270msun
    n = 4 # half number of simulations to plot

    # Text format
    linewidth = 1.5
    fontsize = 14

    # Set up plot environment
    setup_plot_env(fontsize, linewidth)

    # Set up figure and subplots
    num_subplots = 2
    fig, axs = create_subplots(num_subplots, xlim, time_cutoff, fontsize, title)

    # Line colours
    c_s1 = extract_colors('viridis', n, portion="middle", start=0.35, end=0.92) # start=0.33, end=0.92 for 10.8msun-no-sn
    c_s2 = extract_colors('magma', n, portion="middle", start=0.3, end=0.85) # start=0.3, end=0.85 for 10.8msun-no-sn
    c = np.concatenate((c_s1, c_s2))

    # Plot data
    data_files = [
        "data_files/data-1S.RSb01.csv", "data_files/data-1S.RSm01.csv",
        "data_files/data-2S.RSb01.csv", "data_files/data-2S.RSm01-2.csv",
        "data_files/data-1S.b01-no-SN.csv", "data_files/data-1S.m01-no-SN.csv",
        "data_files/data-2S.b01-no-SN.csv", "data_files/data-2S.m01-no-SN.csv"
    ]

    j = 0
    alpha2 = 0.8
    for data_file in data_files:
        label = extract_label(data_file)
        n = 20 if 'b' in label else 10
        plot_extra_line_mass_growth(j, data_file=data_file, label_extra=label, alpha=alpha2, n=n)
        j += 1

    accrate_line = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth)]

    # Include legends and save the plot
    axs[0].legend(fontsize=fontsize-4, ncol=2, loc="lower right", handlelength=0.7) # "lower right" for no-sn
    axs[1].legend(accrate_line, [r"$\rm \dot{M}_{Edd}$"], loc="lower right", fontsize=fontsize-2.2, ncol=1)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.7, 4.7)
    #plot_name = 'mass_growth-1S+2S-no-SN' + '.pdf'
    plot_name = sys.argv[-1] + '.pdf'
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + plot_name)
