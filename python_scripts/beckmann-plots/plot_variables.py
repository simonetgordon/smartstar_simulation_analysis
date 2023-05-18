import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np
import seaborn as sns
from matplotlib import rc

##########################################################################################################
#                                     Plot BHL Variables vs time
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
##########################################################################################################

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
    return next(j for j, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


# function that does the interpolation
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

    x = "s1-40msun-" # plot number
    y = sys.argv[-1] # naming plot
    xlim = 0.225
    linewidth = 1.5
    fontsize = 12

    # text format
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth

    # set up figure
    fig = plt.figure()
    num_subplots = 6
    fig, axs = plt.subplots(num_subplots, 1, sharex=True)


    c1 = ['#d46a7e', 'lightblue', 'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']
    # purple/lavender, greeny-blue, cameo green, dusky-pink/darker pink, mustard yellow, classic blues, fuschia/burgundy
    c2 = ['#856798', '#0a888a', '#56ae57', '#ba6873', '#ffc512', '#436bad', '#9d0759', '#9d0216', '#7bc8f6', '#d0c101',
        '#c4387f','#7bb274', '#06b48b', '#6c3461']
    #c = sns.color_palette("Paired", len(bhl_object_list*8))
    c0 = ['#9d0216', '#E66101', '#436bad', '#A6611A']
    #c = sns.diverging_palette(140, 300, s=100, l=50, sep=5, center="light", n=len(bhl_object_list*2))
    #c = sns.color_palette("spectral", len(bhl_object_list*2), desat=1)
    #c = sns.blend_palette(['blueviolet', 'orange', 'springgreen', "#2CA02C", '#17BECF', "#e377c2", "#76549A"], len(bhl_object_list*2))
    c = ['blueviolet', 'turquoise', 'limegreen', 'darkgreen']

    l = tidy_data_labels(bhl_object_labels)
    j = 0
    alpha = 0.9
    time_cutoff = 0.42 # Myrs
    i_start = 1
    for i, BHL in enumerate(bhl_object_list):
        # convert ages from yrs to Myrs
        BHL.ages = np.array(BHL.ages) / 1e6

        # find index of age that matches end age of time limit
        i_age = first_index(BHL.ages[i_start:], time_cutoff, rtol=1e-5, atol=0.2) 
        
        window_size = 50
        age = movingaverage(BHL.ages[i_start:i_age], window_size)
        mass =  movingaverage(BHL.mass[i_start:i_age], window_size)
        accrate = movingaverage(BHL.accrates[i_start:i_age], window_size)
        density = movingaverage(BHL.average_density[i_start:i_age], window_size)
        avg_vinf = movingaverage(BHL.average_vinfinity[i_start:i_age], window_size)
        avg_cinf = movingaverage(BHL.average_cinfinity[i_start:i_age], window_size)
        hl_radius = movingaverage(BHL.hl_radius[i_start:i_age], window_size)
        bondi_radius = movingaverage(BHL.bondi_radius[i_start:i_age], window_size)
        jeans = movingaverage(BHL.jeans_length[i_start:i_age], window_size)

        # 1) BH Mass
        axs[0].plot(age, mass, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

        # 2) Accretion Rates
        axs[1].plot(age, accrate, color=c[j], linestyle='solid', label=l[i], alpha=alpha)
        axs[1].plot(age, eddington_rate(mass), color=c[j], linestyle='dotted', label=l[i], alpha=alpha)

        # 3) Densities
        axs[2].plot(age, density, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

        # 4) Velocities
        axs[3].plot(age, avg_vinf, color=c[j], linestyle='solid', label=l[i]+'-vinf', alpha=alpha)
        axs[3].plot(age, avg_cinf, color=c[j], linestyle='dotted', label=l[i]+'-cinf', alpha=alpha)

        # 5) HL radius
        axs[4].plot(age, hl_radius, color=c[j], linestyle='solid', label=l[i], alpha=alpha)
        axs[4].plot(age, bondi_radius, color=c[j], linestyle='dotted', label=l[i], alpha=alpha)

        # 6) Jeans length
        axs[5].plot(age, jeans, color=c[j], linestyle='solid', label=l[i]+'-jeans-length', alpha=alpha)

        j += 1

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
    axs[3].set_ylabel(r"$\rm \nu \, (km/s)$", fontdict=None)
    axs[4].set_ylabel(r"$\rm r \, (pc)$", fontdict=None)
    axs[5].set_ylabel(r"$\rm r_{jeans} \, (pc)$", fontdict=None)
    axs[3].set_yscale('linear')
    axs[0].set_yscale('log')
    #axs[0].set_title(str(x) + str(y), fontdict=None)
    for i in [4, 5]:
        if x == "s1-270msun-":
            dx = [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04]
        elif x == "s1-40msun-":
            dx = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
        elif x == "s2-270msun-":
            dx = [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
        elif x == "s2-40msun-":
            dx = [8.4e-03, 4.2e-03, 2.1e-03]
        elif x == "s2-40msun-2-":
            dx = [1.3e-03, 5.2e-04, 1.3e-04]

        dx_1s = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
        c1 = c3 = 'lightcoral'
        c2 = 'indianred'
        l1 = 'dashdot'
        l2 = 'dashdot'
        l3 = 'dotted'
        alpha_dx = 0.5
        axs[i].axhline(y=dx[0], color=c[0], linestyle=l1, lw=linewidth,  label="dx = " + str(dx[0]) + "pc", alpha=alpha_dx)
        axs[i].axhline(y=dx[1], color=c[1], linestyle=l1, lw=linewidth,  label="dx = " + str(dx[1]) + "pc", alpha=alpha_dx)
        axs[i].axhline(y=dx[2], color=c[2], linestyle=l1, lw=linewidth, label="dx = " + str(dx[2]) + "pc", alpha=alpha_dx)
        axs[i].axhline(y=dx[3], color=c[3], linestyle=l1, lw=linewidth, label="dx = " + str(dx[3]) + "pc", alpha=alpha_dx)
        #axs[i].axhline(y=0.00077, color=c[6], linestyle='solid', label="dx = 7.7e-04 pc")
    axs[5].set_xlabel(r"BH Age (Myr)", fontdict=None)

    if xlim == 1:
        axs[0].set_ylim([0, 2800])
        axs[1].set_ylim([2e-8, 9e-3])
        axs[2].set_ylim([7e3,8e7])
        axs[3].set_ylim([0.1, 11])
        axs[5].set_ylim([8e-4, 12])
    if xlim == 3:
        axs[0].set_ylim([0, 8000])
        #axs[1].set_ylim([5e-5, 2e-2])
        #axs[2].set_ylim([8e3, 3e8])
        #axs[3].set_ylim([0, 25])


    # legends
    if x == "s1-270msun-":
        dx = ["dx = 0.012 pc", "dx = 3.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
    elif x == "s1-40msun-":
        dx = ["dx = 2.5e-02 pc", "dx = 1.2e-02 pc", "dx = 3.1-03 pc", "dx = 7.7e-04 pc"]
    elif x == "s2-270msun-":
        dx = ["dx = 8.3e-03 pc", "dx = 2.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
    elif x == "s2-40msun-":
        dx = ["dx = 8.4e-03 pc", "dx = 4.2e-03 pc", "dx = 2.1e-03 pc"]
    elif x == "s2-40msun-2":
        dx = ["dx = 1.3e-03 pc", "dx = 5.2e-04 pc", "dx = 1.3e-04  pc"]

    dx_lines = [Line2D([0], [0], color=c[0], lw=linewidth, linestyle=l1),
                Line2D([0], [0], color=c[1], lw=linewidth, linestyle=l1),
                Line2D([0], [0], color=c[2], lw=linewidth, linestyle=l1),
                Line2D([0], [0], color=c[3], lw=linewidth, linestyle=l1),
                ]
    vel_lines = [Line2D([0], [0], color='grey', lw=linewidth),
                Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
    radius_lines = [Line2D([0], [0], color='grey', lw=linewidth),
                    Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
    accrate_line = [Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
    axs[0].legend(fontsize=fontsize-1, ncol=1)  # upper/lower
    axs[1].legend(accrate_line, [r"$\rm \dot{M}_{Edd}$"], loc="lower right", fontsize=fontsize-1, ncol=2)  # upper/lower
    axs[3].legend(vel_lines, [r"$\rm \nu_{\infty}$", r"\rm $c_{\infty}$"], loc="upper left", fontsize=fontsize-1, ncol=1)  # upper/lower
    axs[4].legend(radius_lines, [r"$\rm r_{HL}$", r"$\rm r_{Bondi}$"], fontsize=fontsize-1, ncol=2)  # upper/lower
    axs[5].legend(dx_lines, dx, fontsize=fontsize-2, ncol=1)

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(5.4, 9)
    plot_name = 'time-' + str(x) + str(y) + '.pdf'
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/", plot_name)
