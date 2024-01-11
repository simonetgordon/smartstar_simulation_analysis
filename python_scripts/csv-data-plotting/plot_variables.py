import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from matplotlib import rc
from matplotlib.lines import Line2D
from helper_functions import adaptive_moving_average, extract_colors, eddington_rate, configure_font, extract_simulation_name_from_csv

##########################################################################################################
#                                     Plot BHL Variables vs time
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
# list data files in order of low res -> high res
##########################################################################################################


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
    elif x == "s1-baseline-fb":
        dx = [1.229940e-02, 1.543e-03] # l14, l17
    else:
        dx = None  # You may want to handle this case differently
    return dx


def plot_lines(i, j, data_file, label, alpha=0.8, n=10):

    # Load the CSV file into a DataFrame
    df = pd.read_csv(data_file)

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
    axs[0].plot(age_av, bh_mass_av, color=c[j], linestyle='solid', label=label, alpha=alpha)

    # 2) Accretion Rates
    axs[1].plot(age_av, accrate_av, color=c[j], linestyle='solid', label=label, alpha=alpha)
    axs[1].plot(age_av, eddington_rate(bh_mass_av), color=c[j], linestyle='dashed', label=label, alpha=alpha)

    # 3) Densities
    axs[2].plot(age_av, density_av, color=c[j], linestyle='solid', label=label, alpha=alpha)

    # 4) Velocities
    axs[3].plot(age_av, vinfinity_av, color=c[j], linestyle='solid', label=label+'-vinf', alpha=alpha)

    # 5) HL radius
    axs[4].plot(age_av, hl_radius_av/dx[i], color=c[j], linestyle='solid', label=label, alpha=alpha)

    return 0

def plot_beckmann_data(axs, alpha=0.8):
    # Beckmann data 1 - D_126_tiny_Figure6
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
    
    # Beckmann data 1 - R_128_Figure11
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



if __name__ == "__main__":

    ###################################### Parameters ######################################

    x = "s1-baseline-fb"    # simulation type
    y = sys.argv[-1]        # naming plot
    xlim = 2                # Myrs
    time_cutoff = xlim      # Myrs
    include_beckmann_data = False
    alpha = 0.9             # transparency of lines
    num_subplots = 5        # number of subplots
    title = None          # plot title
    linewidth = 1.5
    fontsize = 12

    ########################################################################################

    # set up figure
    configure_font(fontsize=fontsize, linewidth=linewidth)
    fig = plt.figure()
    fig, axs = plt.subplots(num_subplots, 1, sharex=True)

    # identify cell widths for this simulation
    dx = identify_cell_widths(x)


    if include_beckmann_data:
        plot_beckmann_data(axs, alpha=alpha)

    """""""""""""""""
       Plot data
    """""""""""""""""

    # set line colours
    c_s1 = extract_colors('magma', 1, portion="middle", start=0.3, end=0.8)
    c_s2 = extract_colors('viridis', 1, portion="middle", start=0.3, end=0.9)
    c = np.concatenate((c_s1, c_s2))

    # Plot data
    # # 1S.b + 2S.b
    # data_files = [
    #     'data_files/data-1S.RSbf4.csv', 'data_files/data-1S.RSb01.csv',
    #     'data_files/data-1S.b01-no-SN.csv', 'data_files/data-1S.RSb04.csv', 'data_files/data-1S.b04-no-SN.csv', 
    #     'data_files/data-2S.RSbf16.csv', 'data_files/data-2S.RSbf4.csv', 
    #     'data_files/data-2S.RSb01.csv',  'data_files/data-2S.b01-no-SN.csv'
    # ]

    # # 1S.m + 2S.m
    # data_files = [
    #     'data_files/data-1S.mf4-no-derefine.csv', 'data_files/data-1S.RSm01.csv',
    #     'data_files/data-1S.m01-no-SN.csv', 'data_files/data-1S.RSm04.csv', 'data_files/data-1S.m04-no-SN.csv',
    #     'data_files/data-2S.RSmf16-2-gap.csv', 'data_files/data-2S.RSmf4-2.csv',
    #     'data_files/data-2S.m01-386+.csv', 'data_files/data-2S.m01-no-SN.csv'
    # ]
    data_files = ['data_files/data-1B.RSb01-2.csv', 
              #'data_files/data-2B.RSb01.csv',
              'data_files/data-1B.resim.th.b01.csv']

    i = j = 0
    alpha2 = 0.8
    for data_file in data_files:
        label = extract_simulation_name_from_csv(data_file)
        n = 30
        plot_lines(i, j, data_file=data_file, label=label, alpha=alpha2, n=n)
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
        axs[i].set_xlim([1.4, xlim+0.01]) # for truncated view

    # ylabels
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"$\rm M_{BH} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm \dot{M} \, (M_{\odot}/yr)$", fontdict=None)
    axs[2].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
    axs[3].set_ylabel(r"$\rm \nu_\infty \, (km/s)$", fontdict=None)
    axs[4].set_ylabel(r"$\rm r_{HL}/dx $", fontdict=None)
    axs[-1].set_xlabel('BH Age (Myr)')

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
        axs[0].set_ylim([900, 3100])
        axs[0].axhline(y=1000, color='grey', linestyle='dotted', linewidth=linewidth+2,alpha=1)
        axs[1].set_ylim([2e-9, 6e-3])
        axs[2].set_ylim([2e2, 2e8])
        axs[3].set_ylim([6, 400])
        axs[4].set_ylim([6e-3, 3])
        axs[4].axhline(y=1, color='grey', linestyle='dashdot', linewidth=linewidth+1,alpha=1)
    elif x == "s1-baseline-fb":
        linewidth = 1.5
        axs[0].set_ylim([900, 3100])
        axs[0].axhline(y=1000, color='grey', linestyle='dotted', linewidth=linewidth+2,alpha=1)
        axs[1].set_ylim([2e-9, 6e-3])
        axs[2].set_ylim([2e2, 2e8])
        axs[3].set_ylim([6, 400])
        axs[4].set_ylim([6e-3, 3])
        axs[4].axhline(y=1, color='grey', linestyle='dashdot', linewidth=linewidth+1,alpha=1)

    ############################### Legends ################################
    # Include legend
    axs[0].legend(fontsize=fontsize-2, ncol=3, loc='upper center', 
                  bbox_to_anchor=(0.5, 1.58), handlelength=1) # 1.3 for m01, 

    ##########################################################################

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.8, 8)
    plot_name = 'time-' + str(x) + str(y) + '.pdf'
    fig.savefig('plots/'+plot_name, bbox_inches='tight')
    print("created plots/"+plot_name)
