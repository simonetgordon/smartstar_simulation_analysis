########################################       Plot IC Attributes     ####################################
#
# to run: python plot_ic_attributes.py radial_profile_ics_halo_1kpc.pdf
##########################################################################################################

import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from smartstar_find import ss_properties
from derived_fields import add_fields_ds
import numpy as np
import yt
import os
from matplotlib import rc


def set_tick_params(axs, n_subplots, fontsize, custom_ticks=[], xlim=[], vline=0, linewidth=2):
    for i in range(n_subplots):
        axs[i].set_xscale('log')  # Set the x-axis scale to logarithmic

        # Add a vertical line at x = vline with customized visual properties
        if vline > 0:
            axs[i].axvline(x=vline, color='grey', linestyle='dashed', lw=linewidth-0.5, alpha=0.8)

        # Set the x-axis limits if provided
        if len(xlim) > 0:
            axs[i].set_xlim(xlim)  
        
        axs[i].tick_params(bottom=True, left=True)  # Set tick parameters for the plot
        axs[i].minorticks_on()  # Enable minor ticks on the plot

        # Set custom tick locations on the x-axis if provided
        if len(custom_ticks) > 0:
            axs[i].set_xticks(custom_ticks)  

        axs[i].tick_params(axis="x", which='minor', length=1, direction="in")  # Set minor tick parameters on the x-axis
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")

        # Set major tick parameters on the y-axis with specified label size and visual properties
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize)  # Set major tick label size on the y-axis
        axs[i].tick_params(axis="y", which='minor', length=2)  # Set minor tick parameters on the y-axis


def critical_density(ds, n_crit=False):
    # Get the cosmological parameters from the dataset parameters
    hubble_constant = yt.YTQuantity(ds.parameters.get("CosmologyHubbleConstantNow")*100, "km/s/Mpc")
    matter_density = yt.YTQuantity(ds.parameters.get("CosmologyOmegaMatterNow"), "")
    cosmological_constant_density = yt.YTQuantity(ds.parameters.get("CosmologyOmegaLambdaNow"), "")

    # Convert the Hubble constant to CGS units (cm/s/Mpc)
    hubble_constant_cgs = (hubble_constant).to("cm/s/Mpc")

    # Conver Hubble constant to 1/s units
    hubble_constant_pers = hubble_constant_cgs / yt.YTQuantity(3.1e19*1e5, "cm/Mpc") 

    # Calculate the critical density
    critical_density = (3 * hubble_constant_pers**2) / (8 * np.pi * yt.physical_constants.G) * (matter_density + cosmological_constant_density)

    # Convert the critical density to the desired units (1/cm^3)
    critical_density = critical_density.to("g/cm**3")

    # Convert to hydrogen nuclei density if desired
    if n_crit:
        critical_density /= yt.physical_constants.mh # /cm**3

    return critical_density


def set_font(fontsize=12, linewidth=1):
    # font settings
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['lines.linewidth'] = linewidth


def generate_labels(dds, root_dir, sim, s=True, tform=0):
    # generate 2 lists: labels, datasets (DS)
    labels = []
    DS = []
    for i, dd in enumerate(dds):
        ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
        add_fields_ds(ds)
        j = i + 1
        time_precision = 4 # no. characters in time str
        if s: # include seed
            label = "s" + str(j) + "_" + str(float(ds.current_time.to('Myr') - tform))[:time_precision] + "Myr"
        else: # include sim 
            label = str(tidy_data_labels(sim)[i]) + "_" + str(float(ds.current_time.to('Myr') - tform))[:time_precision] + "Myr"
        DS.append(ds)
        labels.append(label)

    return labels, DS


def generate_DS(dds, root_dir, sim):
    # make list of complete ds filepaths
    DS = []
    for i, dd in enumerate(dds):
        ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
        add_fields_ds(ds)
        DS.append(ds)
    return DS


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


if __name__ == "__main__":

    ################################## Parameters ##################################
    rvir_pc_s1 = 79.34 # at 2167.60571289 comoving. This is at z = 26.338 at l3
    rvir_pc_s2 = 187.0971 # 3600.0002 comoving. This is at z = 19.115 at l2
    M200 = 2.7491550e+05 # at z = 26.338 at l3 (same as Mvir)
    M200_s2 = 1292000.0 # at z = 19.115 at l2
    n_subplots = 5
    y = sys.argv[-1] # naming plot
    fontsize = 10 # for projection annotations
    linewidth = 2
    r_lim_kpc = 10 # kpc
    alpha = 0.9

    root_dir = ["/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                "/disk14/sgordon/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/"]
    sim = ["1B.RSb16", "2B.RSb16"]
    dds = ["DD0128/DD0128", "DD0198/DD0198"]

    # line colours
    c2 = ['cyan', 'salmon', 'salmon', 
        'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']
    #c = ['#415BFF', '#FF1F71']
    c = [(0.1294, 0.5686, 0.5490), (0.7176, 0.2157, 0.4745)] # teal, magenta
    ################################################################################

    # generate line labels
    labels, DS = generate_labels(dds, root_dir, sim)

    # set font format
    set_font(fontsize, linewidth)

    # set up plot
    fig = plt.figure()
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)

    # create profiles
    for i, ds in enumerate(DS):
        ss_pos = ss_properties(ds)[0]
        sp = ds.sphere(ss_pos, (r_lim_kpc, "kpc"))

        # create weight_field=None profile for cumulative sum rp
        rp = yt.create_profile(
            sp,
            "radius",
            [("gas", "mass"), ("gas", "temperature"),
            ("deposit", "all_mass"),
            #("gas", "theta_vel_dynamical_timescale"), 
            ],
            units={"radius": "pc", ("gas", "mass"): "Msun", ("deposit", "all_mass"):"Msun"},
            logs={"radius": True, ("gas", "mass"): True, ("gas", "temperature"): True, ("deposit", "all_mass"):True},
            weight_field=None,
            accumulation=True # -True to reverse direction
        )

        # create weight_field=cell_mass profile for non-cumulative rp
        rp2 = yt.create_profile(
            sp,
            "radius",
            [("gas", "temperature"), ("gas", "mass"), ("gas", "dynamical_time"), ("gas", "H2_p0_fraction"),
            ("gas", "number_density"), ("gas", "blackhole_freefall_timescale"), ("gas", "density"),
            ('enzo', 'Dark_Matter_Density')],
            units={"radius": "pc", ("gas", "mass"): "Msun", ('enzo', 'Dark_Matter_Density'): "g/cm**3"},
            logs={"radius": True, ("gas", "temperature"): True}
        )

         # create weight_field=density profile for non-cumulative rp
        rp3 = yt.create_profile(
            sp,
            "radius",
            [("gas", "H2_p0_fraction"), ("gas", "H_p0_fraction"), ("gas", "density"), ("gas", "radial_velocity")],
            units={"radius": "pc", "radial_velocity": "km/s"},
            logs={"radius": True},
            weight_field=('gas', 'density')
        )

        # total mass
        axs[0].loglog(rp.x.value, rp[("gas", "mass")].value + rp[("deposit", "all_mass")].value,
                color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # temperature
        axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # total density
        mh = 1.6735575e-24*yt.units.g # g
        axs[2].loglog(rp2.x[rp2.used], (rp2[("gas", "density")][rp2.used] + rp2[("enzo", "Dark_Matter_Density")][rp2.used])/mh,
                      color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # dynamical timescale
        axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "blackhole_freefall_timescale")][rp2.used].to("Myr"),
                    color=c[i], linestyle='solid', label=r"BH ff time $\propto M$", alpha=alpha)

        # h2 fraction
        # axs[4].loglog(rp3.x[rp3.used], rp3[("gas", "H2_p0_fraction")][rp3.used],
        #             color=c[i], linestyle='solid', label=labels[i], alpha=alpha)
        
        # radial velocity
        axs[4].plot(rp3.x[rp3.used], rp3[("gas", "radial_velocity")][rp3.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)
    # set ticks
    xticks = np.logspace(-2, 4, 7)
    set_tick_params(axs, n_subplots, fontsize, custom_ticks=xticks, xlim=[2e-3, r_lim_kpc*1e3])

    # make lines for legend
    lines = [Line2D([0], [0], color=c[0], linestyle='solid', lw=linewidth),
                Line2D([0], [0], color=c[1], linestyle='solid', lw=linewidth),
                Line2D([0], [0], color='grey', linestyle='dashdot', lw=linewidth),
                Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
    labels = ['Halo1_124.7Myr', 'Halo2_195.5Myr', r'$R_{200}$', r'$M_{200}$']
    # set axis labels
    axs[-1].set_xlabel(r"$\rm Radius \, (pc)$", fontdict=None)
    axs[4].set_ylim([-7, 1])
    axs[4].set_ylabel(r"$\nu_r$ (km/s)", fontdict=None)
    axs[4].set_yscale('linear')
    axs[2].set_ylabel(r"$\rm n \, (1/cm^{3})$", fontdict=None)
    for i in range(n_subplots):
        axs[i].axvline(x=rvir_pc_s1, color=c[0], linestyle='dashdot', lw=linewidth-0.5, alpha=0.7, label=r"$R_{200}$")
        axs[i].axvline(x=rvir_pc_s2, color=c[1], linestyle='dashdot', lw=linewidth-0.5, alpha=0.7, label=r"$R_{200}$")
    axs[3].set_ylabel(r"$t_{\rm ff} \, \rm (Myr)$", fontdict=None)
    axs[1].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
    axs[0].set_ylabel(r"$\rm M_{encl} \, (M_{\odot})$", fontdict=None)
    axs[0].axhline(y=M200, color=c[0], linestyle='dotted', lw=linewidth-0.5, alpha=0.8, label=r"$M_{200}$")
    axs[0].axhline(y=M200_s2, color=c[1], linestyle='dotted', lw=linewidth-0.5, alpha=0.8, label=r"$M_{200}$")
    axs[0].legend(lines, labels, loc="upper left", fontsize=fontsize-1, ncol=1)  # upper/lower
    for i in range(n_subplots):
        axs[i].set_xlim([2e-3, 1e3])
        axs[i].grid(which='major', linestyle='solid', linewidth=0.5, color='grey', alpha=0.5)
        axs[i].tick_params(axis="x", which='major', length=3, direction="in")
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
    #axs[0].set_title("Halo properties at time of BH formation", fontsize=fontsize, fontdict=None)
    #axs[3].axhline(y=360*critical_density(ds), color='grey', linestyle='dotted', lw=linewidth, alpha=1, label=r"$n_{200}$")

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(3.6, 6.2)
    plot_name = 'radial_profile_ics_halo_{}kpc.pdf'.format(r_lim_kpc)
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + str(plot_name))
