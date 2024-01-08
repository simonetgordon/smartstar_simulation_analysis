import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.cm as cm
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d
import re
#import seaborn as sns

##########################################################################################################
#                                     Plot BHL Variables vs time
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
# list data files in order of low res -> high res
##########################################################################################################


def tidy_data_labels(labels: str or list):
    # for lists of labels
    if len(labels) < 50:
        data_labels = [i.replace("-2", "") for i in labels]
        data_labels = [i.replace("RS", "") for i in data_labels]
        data_labels = [i.replace("-gap", "") for i in data_labels]
        data_labels = [i.replace("-4dx", "") for i in data_labels]
        
    # for single label
    else:
        data_labels = labels.replace("-2", "")
        data_labels = data_labels.replace("RS", "")
        data_labels = data_labels.replace("-gap", "")
        data_labels = data_labels.replace("-4dx", "")
    return data_labels


def eddington_rate(mparticle_msun: float):
    # eddington rate constants in cgs units
    PI = 3.1415927
    GravConst = 6.6740831e-8
    mh = 1.67262171e-24
    clight = 2.99792458e10
    sigma_thompson = 6.65E-25
    eta_disk = 0.1 # optically thick and e_radiative = 0.11
    SolarMass = 1.9891e33
    yr_s = 3.1556952E7
    mparticle_g = mparticle_msun*SolarMass
    mdot_edd = 4.0 * PI * GravConst * mparticle_g * mh / (clight * eta_disk * sigma_thompson)
    return mdot_edd*yr_s/SolarMass


def identify_cell_widths(x: str):
    if x == "s1-270msun-":
        dx = [0.04919164, 1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04]
    elif x == "beck-comp-":
        dx = [1.229791e-02, 3.074475e-03]
    elif x == "s1-40msun-":
        dx = [0.09839468, 2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
    elif x == "s2-270msun-":
        dx = [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
    elif x == "s2-40msun-":
        dx = [8.4e-03, 2.1e-03, 5.2e-04, 1.3e-04]
    elif x == "s2-40msun-2-":
        dx = [1.3e-03, 5.2e-04, 1.3e-04]
    elif x == "s1+s2-270msun-":
        dx = [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04, 8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
    elif x == "s1+s2-40msun-":
        #dx = [2.459867e-02, 1.229940e-02, 3.074829e-03, 3.074829e-03, 7.692833e-04, 7.692833e-04, 8.4e-03, 2.1e-03, 5.2e-04, 1.3e-04, 1.3e-04]
        dx = [1.229940e-02, 3.074829e-03, 3.074829e-03, 7.692833e-04, 7.692833e-04, 8.4e-03, 5.2e-04, 1.3e-04, 1.3e-04] # for 1Sm/2Sm
        #dx = [1.229940e-02, 3.074829e-03, 3.074829e-03, 7.692833e-04, 8.4e-03, 5.2e-04, 1.3e-04, 1.3e-04, 7.692833e-04,] # for 1Sb/2Sb
    else:
        dx = None  # You may want to handle this case differently
    return dx


def extract_colors(cmap_name, n, portion=None, start=None, end=None):
    cmap = cm.get_cmap(cmap_name)

    if start is not None and end is not None:
        values = np.linspace(start, end, n)
    elif portion == "beginning":
        values = np.linspace(0, 0.3, n)
    elif portion == "middle":
        values = np.linspace(0.3, 0.95, n)
    elif portion == "end":
        values = np.linspace(0.7, 1, n)
    elif portion is None:
        values = np.linspace(0, 1, n)
    else:
        raise ValueError("Invalid portion specified. Use 'beginning', 'middle', 'end', or None.")

    colors = cmap(values)
    return colors


def adaptive_moving_average(data, window_size=5):
    """
    Compute an adaptive moving average on the logarithm of the data.
    
    :param data: The input data (list or array).
    :param window_size: The window size for the moving average.
    :return: An array of the moving average values in the original data space.
    """
    # Take the logarithm of data, excluding non-positive values
    log_data = np.log(data[data > 0])
    data_length = len(log_data)
    log_moving_avg = np.zeros(data_length)

    for i in range(data_length):
        start = max(i - window_size // 2, 0)
        end = min(i + window_size // 2 + 1, data_length)
        log_moving_avg[i] = np.mean(log_data[start:end])

    # Exponentiate to return to original data space
    moving_avg = np.exp(log_moving_avg)

    # Handle edge cases if original data had non-positive values
    moving_avg_full = np.full_like(data, np.nan)
    positive_indices = np.where(data > 0)[0]
    moving_avg_full[positive_indices] = moving_avg

    return moving_avg_full


def extract_label(file_path):
    # Extract the base file name without the path and extension
    file_name = file_path.split('/')[-1].split('.csv')[0]

    # Use regex to find the label pattern in the file name, excluding 'RS'
    match = re.search(r'(\dS)\.(?:RS)?([a-zA-Z0-9+]+(?:-no-SN)?)', file_name)
    if match:
        return match.group(1) + '.' + match.group(2)
    else:
        return 'Unknown'


def plot_lines(i, j, data_file, label_extra, alpha=0.8, n=10):

    # Load the CSV file into a DataFrame
    df = pd.read_csv(data_file)
    #group_labels = np.arange(len(df)) // window_size
    #df = df.groupby(group_labels).mean().reset_index(drop=True)

    # Extract the columns you're interested in
    age = df['age'].values/1e6
    bh_mass = df['BH mass'].values
    accrate = df['accrate'].values
    density = df['average density'].values
    vinfinity = df['average vinfinity'].values
    cinfinity = df['average cinfinity'].values
    hl_radius = df['HL radius'].values

    # try adaptively smoothing the data
    age_av = adaptive_moving_average(age, window_size=n)
    bh_mass_av = adaptive_moving_average(bh_mass, window_size=n)
    accrate_av = adaptive_moving_average(accrate, window_size=n)
    density_av = adaptive_moving_average(density, window_size=n)
    vinfinity_av = adaptive_moving_average(vinfinity, window_size=n)
    cinfinity_av = adaptive_moving_average(cinfinity, window_size=n)
    hl_radius_av = adaptive_moving_average(hl_radius, window_size=n)

    # 1) BH Mass
    axs[0].plot(age_av, bh_mass_av, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)

    # 2) Accretion Rates
    axs[1].plot(age_av, accrate_av, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)
    axs[1].plot(age_av, eddington_rate(bh_mass_av), color=c[j], linestyle='dashed', label=label_extra, alpha=alpha)

    # 3) Densities
    axs[2].plot(age_av, density_av, color=c[j], linestyle='solid', label=label_extra, alpha=alpha)

    # 4) Velocities
    axs[3].plot(age_av, vinfinity_av, color=c[j], linestyle='solid', label=label_extra+'-vinf', alpha=alpha)

    # 5) HL radius
    axs[4].plot(age_av, hl_radius_av/dx[i], color=c[j], linestyle='solid', label=label_extra, alpha=alpha)

    return 0


if __name__ == "__main__":

    ###################################### Parameters ######################################

    x = "s1+s2-40msun-"       # simulation type
    y = sys.argv[-1]        # naming plot
    xlim = 1                # Myrs
    time_cutoff = xlim      # Myrs
    i_start = 0             # start index
    include_beckmann_data = False
    alpha = 0.9        # transparency of lines
    num_subplots = 5        # number of subplots
    smooth_simulations = len(bhl_object_list)  # number of simulations to smooth (starting from last sim)
    window = 30            # window size to average over
    rtol=1e-6   
    atol=1e-4               # 1e-4 for 1S.m, 2e-4 for 1B.m
    title = None          # plot title
    extra_line = True      # plot extra line from data file

    ########################################################################################

    # text format
    linewidth = 1.5
    fontsize = 12
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth

    # set up figure
    fig = plt.figure()
    fig, axs = plt.subplots(num_subplots, 1, sharex=True)

    # identify cell widths for this simulation
    dx = identify_cell_widths(x)

    """""""""""""""""""""
    1) Plot Beckmann data
    """""""""""""""""""""

    if include_beckmann_data:
        # make array of Beckmann data
        beck_data_fp = '/cephfs/sgordon/smartstar_simulation_analysis/python_scripts/beckmann-data/D_126_tiny_Figure6'
        csv_files = ['mass.csv', 'accretion_rate.csv', 'number_density.csv', 'velocity.csv', 'radii.csv']
        beck_data_arr = []
        l_beck = "D_126_tiny" # label
        
        # Loop through each CSV file
        for k, file_name in enumerate(csv_files):
            # Create the full file path
            file_path = beck_data_fp + '/' + file_name
            
            # Load CSV file into DataFrame
            df = pd.read_csv(file_path)
            
            # Extract NumPy array and append to the list
            beck_data_arr.append(df.to_numpy())
            beck_data = df.to_numpy()

            axs[k].plot(beck_data[:,0], beck_data[:, 1], color="darkblue", linestyle='solid', label=l_beck, alpha=alpha)

        csv_files = ['mass.csv','accretion_rate.csv', 'number_density.csv']
        beck_data_fp = '/cephfs/sgordon/smartstar_simulation_analysis/python_scripts/beckmann-data/R_128_Figure11'
        beck_data_arr_2 = []
        l_beck = "R_128" # label

        # Loop through each CSV file
        for k, file_name in enumerate(csv_files):
            # Create the full file path
            file_path = beck_data_fp + '/' + file_name
            
            # Load CSV file into DataFrame
            df = pd.read_csv(file_path)
            
            # Extract NumPy array and append to the list
            beck_data_arr_2.append(df.to_numpy())
            beck_data = df.to_numpy()

            axs[k].plot(beck_data[:,0], beck_data[:, 1], color="royalblue", linestyle='solid', label=l_beck, alpha=alpha)


    """""""""""""""""
    2) Plot my data
    """""""""""""""""

    # set line colours
    c_s1 = extract_colors('magma', 5, portion="middle", start=0.3, end=0.8)
    c_s2 = extract_colors('viridis', 4, portion="middle", start=0.3, end=0.9)
    c = np.concatenate((c_s1, c_s2))

    # Plot data
    # 1S.b + 2S.b
    data_files = [
        'data_files/data-1S.RSbf4.csv', 'data_files/data-1S.RSb01.csv',
        'data_files/data-1S.b01-no-SN.csv', 'data_files/data-1S.RSb04.csv', 'data_files/data-1S.b04-no-SN.csv', 
        'data_files/data-2S.RSbf16.csv', 'data_files/data-2S.RSbf4.csv', 
        'data_files/data-2S.RSb01.csv',  'data_files/data-2S.b01-no-SN.csv'
    ]

    # 1S.m + 2S.m
    data_files = [
        'data_files/data-1S.mf4-no-derefine.csv', 'data_files/data-1S.RSm01.csv',
        'data_files/data-1S.m01-no-SN.csv', 'data_files/data-1S.RSm04.csv', 'data_files/data-1S.m04-no-SN.csv',
        'data_files/data-2S.RSmf16-2-gap.csv', 'data_files/data-2S.RSmf4-2.csv',
        'data_files/data-2S.m01-386+.csv', 'data_files/data-2S.m01-no-SN.csv'
    ]

    i = j = 0
    alpha2 = 0.8
    for data_file in data_files:
        label = extract_label(data_file)
        n = 20
        plot_lines(i, j, data_file=data_file, label_extra=label, alpha=alpha2, n=n)
        j += 1
        i += 1

    ############################### Format plots #################################

    for i in range(num_subplots):
        axs[i].set_xticks(np.arange(0.1, time_cutoff+0.1, 0.1))
        axs[i].minorticks_on()
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(0.02))
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1.5, length=3, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize-1)
        axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)
        axs[i].set_yscale('log')
        axs[i].set_xlim([0.01, xlim+0.01]) # for truncated view

    # ylabels
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[2].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
    axs[3].set_ylabel(r"$\rm \nu_\infty \, (km/s)$", fontdict=None)
    axs[4].set_ylabel(r"$\rm r_{HL}/dx $", fontdict=None)
    axs[-1].set_xlabel('BH Age (Myr)')

    dx_1s = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
    c1 = c3 = 'lightcoral'
    c2 = 'indianred'
    l1 = 'dashdot'
    l2 = 'dashdot'
    l3 = 'dotted'
    alpha_dx = 0.5

    # set ylims
    if x == "s1-40msun-":
        #axs[0].set_ylim([0, 80])
        #axs[1].set_ylim([5e-5, 2e-2])
        #axs[2].set_ylim([8e3, 3e8])
        axs[3].set_ylim([3e-2, 8e1])
    elif x == "s1-270msun-":
        axs[4].set_ylim([3e-2, 2e2])
        axs[5].set_ylim([-1.5e-2, 0.9e-2])
    elif x == "s2-270msun-":
        axs[4].set_ylim([2e-3, 9e1])
        #axs[5].set_ylim([-1.5e-2, 0.9e-2])
    elif x == "s1+s2-40msun-":
        axs[0].set_ylim([10.1, 3000])
        axs[1].set_ylim([7e-10, 4e-2])
        axs[2].set_ylim([18, 2e10])
        axs[3].set_ylim([1.1e-1, 30])
        axs[4].set_ylim([2e-3, 2e2])
        axs[0].axhline(y=60, color='grey', linestyle='dotted', linewidth=linewidth+2, alpha=alpha)
        axs[4].axhline(y=1, color='grey', linestyle='dotted', linewidth=linewidth+2,alpha=alpha)
    elif x == "s1+s2-270msun-":
        linewidth = 1.5
        axs[0].set_ylim([250, 3100])
        axs[0].axhline(y=1000, color='grey', linestyle='dotted', linewidth=linewidth+2,alpha=1)
        axs[1].set_ylim([3e-4, 2e-2])
        axs[2].set_ylim([2e4, 6e9])
        axs[3].set_ylim([1.1, 6e1])
        axs[4].set_ylim([9e-2, 7e1])
        axs[4].axhline(y=1, color='grey', linestyle='dashdot', linewidth=linewidth+1,alpha=1)

    ############################### Legends ################################

    # Set cell widths to be plotted
    if x == "s1-270msun-":
        dx = ["dx = 0.012 pc", "dx = 3.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
    elif x == "beck-comp-":
        dx = ["dx = 1.2e-02 pc", "dx = 3.1e-03 pc"]
    elif x == "s1-40msun-":
        dx = ["dx = 2.5e-02 pc", "dx = 1.2e-02 pc", "dx = 3.1-03 pc", "dx = 7.7e-04 pc"]
    elif x == "s2-270msun-":
        dx = ["dx = 8.3e-03 pc", "dx = 2.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
    elif x == "s2-40msun-":
        dx = ["dx = 8.4e-03 pc", "dx = 2.1e-03 pc", "dx = 5.2e-04 pc", "dx = 1.3e-04 pc"]
    elif x == "s2-40msun-2":
        dx = ["dx = 1.3e-03 pc", "dx = 5.2e-04 pc", "dx = 1.3e-04  pc"]

    # Configure the legend lines: color, linestyle, linewidth
    dx_lines = [Line2D([0], [0], color=c[0], lw=linewidth, linestyle=l1),
                Line2D([0], [0], color=c[1], lw=linewidth, linestyle=l1),
                # Line2D([0], [0], color=c[2], lw=linewidth, linestyle=l1),
                # Line2D([0], [0], color=c[3], lw=linewidth, linestyle=l1),
                ]
    vel_lines = [Line2D([0], [0], color='grey', lw=linewidth),
                Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
    radius_lines = [Line2D([0], [0], color='grey', lw=linewidth),
                    #Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)
                    ]
    accrate_line = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth)]

    # Include legends
    axs[0].legend(fontsize=fontsize-2, ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.58), handlelength=1) # 1.3 for m01, 

    ##########################################################################

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.8, 8)
    plot_name = 'time-' + str(x) + str(y) + '.pdf'
    fig.savefig('plots/'+plot_name, bbox_inches='tight')
    print("created plots/"+plot_name)
