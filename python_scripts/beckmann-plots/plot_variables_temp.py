import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np
import pandas as pd
from matplotlib import rc
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d
import matplotlib.ticker as ticker
import yt

##########################################################################################################
#                               Plot BHL Variables (+ temperature) vs time
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
# list data files in order of low res -> high res
##########################################################################################################


def resample_data(accretion_rates, times, common_time, smooth_simulations=2, window_size=3):
    resampled_acc_rates = []
    for i in range(len(accretion_rates)):
        # Create an interpolator for each simulation
        interpolator = interp1d(times[i], accretion_rates[i], kind='linear', fill_value='extrapolate')

        # Resample the accretion rates onto the common time grid
        resampled_acc_rates.append(interpolator(common_time))

    # Smooth the specified number of simulations with a moving average
    if smooth_simulations > 0 and len(accretion_rates) >= smooth_simulations:
        smooth_indices = range(-1, -smooth_simulations-1, -1)  # Indices of the last 'smooth_simulations' simulations (reversed)
        for idx in smooth_indices:
            resampled_acc_rates[idx] = convolve1d(resampled_acc_rates[idx], np.ones(window_size)/window_size, mode='reflect')

    return resampled_acc_rates


def tidy_data_labels(labels):
    # for lists of labels
    if len(labels) < 50:
        data_labels = [i.replace("-2", "") for i in labels]
        data_labels = [i.replace("RS", "") for i in data_labels]
    # for single label
    else:
        data_labels = labels.replace("-2", "")
        data_labels = data_labels.replace("RS", "")
    return data_labels


def first_index(a, val, rtol=0.1, atol=10):
    return next(m for m, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


def interpolate_data(arr, N=1000):
    # interpolate arr over N evenly spaced points
    min_val = np.min(arr)
    max_val = np.max(arr)

    t_orig = np.linspace(min_val, max_val, len(arr))
    t_interp = np.linspace(min_val, max_val, N)
    f = interp1d(x=t_orig, y=arr)
    interp_arr = f(t_interp)
    return interp_arr


def movingaverage(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)


def eddington_rate(mparticle_msun):
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


if __name__ == "__main__":

    ###################################### Parameters ######################################

    x = "s1-40msun-"        # simulation type
    y = sys.argv[-1]        # naming plot
    xlim = 1             # Myrs
    time_cutoff = xlim      # Myrs
    i_start = 0             # start index
    include_beckmann_data = False
    alpha = 0.9             # transparency of lines
    num_subplots = 6        # number of subplots
    smooth_simulations = 4  # number of simulations to smooth (starting from last sim)
    window = 10              # window size to average over
    rtol=1e-6   
    atol=2e-4               # 1e-4 for 1S.m, 
    title = '1S.b Local Properties Comparison'

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
    for i in [4]:
        if x == "s1-270msun-":
            dx = [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04]
        elif x == "beck-comp-":
            dx = [1.229791e-02, 3.074475e-03]
        elif x == "s1-40msun-":
            dx = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
        elif x == "s2-270msun-":
            dx = [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
        elif x == "s2-40msun-":
            dx = [8.4e-03, 2.1e-03, 5.2e-04, 1.3e-04]
        elif x == "s2-40msun-2-":
            dx = [1.3e-03, 5.2e-04, 1.3e-04]

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
    c = ['blueviolet', 'turquoise', 'limegreen', 'darkgreen']

    # set BHL properties parameters
    l = tidy_data_labels(bhl_object_labels)
    j = 0

    # try new resampling method
    times = [BHL.ages/1e6 for BHL in bhl_object_list]

    # Find the minimum and maximum times among all simulations
    min_time = min([min(t) for t in times])
    min_final_time = min([t[-1] for t in times])

    # Define the common time grid
    num_steps = int(len(times[-1])/5)  # take time steps from last (most refined) simulation
    print("num points: ", num_steps)
    common_time = np.linspace(min_time, min_final_time, num_steps)

    # resample BHL attributes over the common_time - produces a list of length 4 (for each simulations), 
    # where list elements = attribute array
    mass = resample_data([BHL.mass for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    accretion_rates = resample_data([BHL.accrates for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    density = resample_data([BHL.average_density for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    avg_vinf = resample_data([BHL.average_vinfinity for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    avg_cinf = resample_data([BHL.average_cinfinity for BHL in bhl_object_list], times, common_time, smooth_simulations=smooth_simulations, window_size=window)
    hl_radius = resample_data([BHL.hl_radius for BHL in bhl_object_list], times, common_time)
    bondi_radius = resample_data([BHL.hl_radius for BHL in bhl_object_list], times, common_time)
    #jeans = resample_data([BHL.jeans for BHL in bhl_object_list], times, common_time)

    # Print the resampled data
    for i in range(len(accretion_rates)):
        print(f"Resampled data for Simulation {i+1}:")
        print("Accretion Rates:", accretion_rates[i])
        print("Time:", common_time)
        print()

    # find gas temperature

    for i in range(len(mass)):

        # 1) BH Mass
        axs[0].plot(common_time, mass[i], color=c[j], linestyle='solid', label=l[i], alpha=alpha)

        # 2) Accretion Rates
        axs[1].plot(common_time, accretion_rates[i], color=c[j], linestyle='solid', label=l[i], alpha=alpha)
        axs[1].plot(common_time, eddington_rate(mass[i]), color=c[j], linestyle='dashed', label=l[i], alpha=alpha)

        # 3) Densities
        axs[2].plot(common_time, density[i], color=c[j], linestyle='solid', label=l[i], alpha=alpha)

        # 3) Densities
        axs[2].plot(common_time, density[i], color=c[j], linestyle='solid', label=l[i], alpha=alpha)

        # 4) Velocities
        axs[3].plot(common_time, avg_vinf[i]/avg_cinf[i], color=c[j], linestyle='solid', label=l[i]+'-Mach', alpha=alpha)
        #axs[3].plot(common_time, avg_vinf, color=c[j], linestyle='solid', label=l[i]+'-vinf', alpha=alpha)

        # 5) HL radius
        axs[4].plot(common_time, hl_radius[i]/dx[i], color=c[j], linestyle='solid', label=l[i], alpha=alpha)
        #axs[4].plot(common_time, bondi_radius, color=c[j], linestyle='dotted', label=l[i], alpha=alpha)
        print("average radius resolution: ", hl_radius[i].mean()/dx[i])

        # 6) Scatter Residual
        residual = (accretion_rates[i] - accretion_rates[-1])
            
        print("Min residual: ", residual.min())
        if i == 2:
            alpha2 = 0.2
        elif i == 1:
            alpha2 = 0.4
        else: 
            alpha2 = 0.6
        axs[5].scatter(common_time, residual, s=3, color=c[j], linestyle='solid', marker='o', label=l[i], alpha=alpha2)
        #axs[5].axhline(y=0, color='grey', linestyle='dashed', label=l[i], alpha=alpha)
        yscale_residual = 'linear'
        #axs[5].set_ylim([1e-8, 1e-1])

        print("=============================")

        j += 1


    ############################### Format plots #################################

    # include title (might remove later)
    axs[0].set_title(title)

    for i in range(num_subplots):
        axs[i].set_xticks(np.arange(0, time_cutoff, 0.1))
        axs[i].minorticks_on()
        axs[i].xaxis.set_minor_locator(plt.MultipleLocator(0.02))
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1.5, length=3, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize-1)
        axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)
        axs[i].set_yscale('log')
        axs[i].set_xlim([0, xlim+0.01]) # for truncated view


    # format plots
    axs[0].set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[2].set_ylabel(r"$\rm n \, (H \, cm^{-3})$", fontdict=None)
    #axs[3].set_ylabel(r"$\rm \nu \, (km/s)$", fontdict=None)
    axs[3].set_ylabel(r"$\rm Mach$", fontdict=None)
    #axs[4].set_ylabel(r"$\rm r_{HL} \, (pc)$", fontdict=None)
    axs[4].set_ylabel(r"$\rm r_{HL}/dx $", fontdict=None)
    #axs[5].set_ylabel(r"$\rm r_{jeans} \, (pc)$", fontdict=None)
    axs[5].set_ylabel(r"$\rm \Delta \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[5].set_yscale(yscale_residual)
    axs[5].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    axs[5].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    # Move the "x 10^-4" label below the tick labels
    axs[5].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axs[0].set_yscale('log')
    #axs[5].yaxis.get_offset_text().set_y(0.01)
    #axs[5].yaxis.offsetText.set_position((0, -0.15))
    offset_text = axs[5].yaxis.get_offset_text()
    offset_text.set_verticalalignment('top')
    offset_text.set_x(-0.121)
    axs[-1].set_xlabel('BH Age (Myr)')
    #axs[0].set_title(str(x) + str(y), fontdict=None)

    dx_1s = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
    c1 = c3 = 'lightcoral'
    c2 = 'indianred'
    l1 = 'dashdot'
    l2 = 'dashdot'
    l3 = 'dotted'
    alpha_dx = 0.5
    #axs[i].axhline(y=dx[0], color=c[0], linestyle=l1, lw=linewidth,  label="dx = " + str(dx[0]) + "pc", alpha=alpha_dx)
    #axs[i].axhline(y=dx[1], color=c[1], linestyle=l1, lw=linewidth,  label="dx = " + str(dx[1]) + "pc", alpha=alpha_dx)
    # axs[i].axhline(y=dx[2], color=c[2], linestyle=l1, lw=linewidth, label="dx = " + str(dx[2]) + "pc", alpha=alpha_dx)
    # axs[i].axhline(y=dx[3], color=c[3], linestyle=l1, lw=linewidth, label="dx = " + str(dx[3]) + "pc", alpha=alpha_dx)
    #axs[i].axhline(y=0.00077, color=c[6], linestyle='solid', label="dx = 7.7e-04 pc")
#axs[5].set_xlabel(r"BH Age (Myr)", fontdict=None)

    #if xlim == 1:
        # axs[4].set_ylim([1.01e-5, 9e-1])
    #     # axs[0].set_ylim([0, 2800])
    #     # axs[1].set_ylim([2e-8, 9e-3])
    #     # axs[2].set_ylim([7e3,8e7])
    #     # axs[3].set_ylim([0.1, 11])
    #     # #axs[5].set_ylim([8e-4, 12])
    if x == "s1-40msun-":
        #axs[0].set_ylim([0, 80])
        #axs[1].set_ylim([5e-5, 2e-2])
        #axs[2].set_ylim([8e3, 3e8])
        axs[3].set_ylim([3e-2, 8e1])
    elif x == "s1-270msun-":

        axs[4].set_ylim([3e-2, 2e2])
        axs[5].set_ylim([-1.5e-2, 0.9e-2])


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
    axs[0].legend(fontsize=fontsize-1, ncol=2, loc="lower right")  # upper/lower
    axs[1].legend(accrate_line, [r"$\rm \dot{M}_{Edd}$"], loc="lower right", fontsize=fontsize-1, ncol=2)
    #axs[3].legend(vel_lines, [r"$\rm \nu_{\infty}$", r"\rm $c_{\infty}$"], loc="upper left", fontsize=fontsize-1, ncol=1)
    #axs[4].legend(radius_lines, [r"$\rm r_{HL}$", r"$\rm r_{Bondi}$"], fontsize=fontsize-1, ncol=1)
    #axs[4].legend(dx_lines, dx, fontsize=fontsize-2, ncol=1)

    ##########################################################################

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(5.4, 9)
    plot_name = 'time-' + str(x) + str(y) + '.pdf'
    fig.savefig('plots/'+plot_name, bbox_inches='tight')
    print("created plots/"+plot_name)
