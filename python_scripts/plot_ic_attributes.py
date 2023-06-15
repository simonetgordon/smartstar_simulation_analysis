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


def generate_labels(dds, root_dir, sim, s=True):
    # generate 2 lists: labels, datasets (DS)
    labels = []
    DS = []
    for i, dd in enumerate(dds):
        ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
        add_fields_ds(ds)
        j = i + 1
        if s: # include seed
            label = "s" + str(j) + "_" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"
        else: # include sim 
            sim = tidy_data_labels(sim)
            label = sim + "_" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"
        DS.append(ds)
        labels.append(label)

    return labels, DS


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
    rvir_pc = 2167.60571289 # at z = 26.338 at l3
    M200 = 2.7491550e+05 # at z = 26.338 at l3 (same as Mvir)
    n_subplots = 4
    y = sys.argv[-1] # naming plot
    fontsize = 12 # for projection annotations
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
    c = ['#415BFF', '#FF1F71']
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
            [("gas", "temperature"), ("gas", "mass"), ("gas", "dynamical_time"), ("gas", "baryon_overdensity_sg"),
            ("gas", "number_density"), ("gas", "blackhole_freefall_timescale"), ("gas", "theta_vel_dynamical_timescale"), ("gas", "density"),
            ('enzo', 'Dark_Matter_Density')],
            units={"radius": "pc", ("gas", "mass"): "Msun", ('enzo', 'Dark_Matter_Density'): "g/cm**3"},
            logs={"radius": True, ("gas", "temperature"): True}
        )

        # rp3 = yt.create_profile(
        #     sp,
        #     "radius",
        #     [("all", "particle_mass") ],
        #      units={"radius": "pc", ("all", "particle_mass"): "Msun"},
        #     logs={"radius": True},
        #     weight_field=None
        # )

        # plot data
        axs[0].loglog(rp.x.value, rp[("gas", "mass")].value + rp[("deposit", "all_mass")].value,
                color=c[i], linestyle='solid', label=labels[i] + "_total", alpha=alpha)
        axs[0].loglog(rp.x.value, rp[("deposit", "all_mass")].value,
                    color=c2[i], linestyle='solid', label= labels[i] + "_DM", alpha=alpha)
        r200 = rp[("deposit", "all_mass")].value.cumsum().max() + rp[("gas", "mass")].value.cumsum().max()
        print("r200 mass: ", r200)

        axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "blackhole_freefall_timescale")][rp2.used].to("Myr"),
                    color=c[i], linestyle='solid', label=r"BH ff time $\propto M$", alpha=alpha)
        # axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "theta_vel_dynamical_timescale")][rp2.used].to("yr"),
        #               color=c[j], linestyle='solid', label=r"$2\pi r / v_\theta $", alpha=alpha)

        # Omega_b = 0.0449
        # axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "baryon_overdensity_sg")][rp2.used],
        #               color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "density")][rp2.used] + rp2[("enzo", "Dark_Matter_Density")][rp2.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)
        axs[3].loglog(rp2.x[rp2.used], rp2[("enzo", "Dark_Matter_Density")][rp2.used],
                    color=c2[i], linestyle='solid', label=labels[i], alpha=alpha)

        # t_dyn is larger when using non-mass-weighted densities
        # axs[2].loglog(rp.x[rp.used], rp[("gas", "dynamical_time")][rp.used].to("yr"),
        #               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)


    # set ticks
    xticks = np.logspace(-2, 5, 8)
    axs[3].set_xticks(xticks)
    for i in range(n_subplots):
        axs[i].tick_params(bottom=True, left=True)
        axs[i].minorticks_on()
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
        axs[i].tick_params(axis="y", which='minor', length=2)
        #axs[i].grid(color='grey', linestyle='solid', linewidth=0.5, alpha=0.7)

    # make lines for legend
    r_lines = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth),
                Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]

    # set axis labels
    axs[n_subplots-1].set_xlabel(r"$\rm Radius \, (pc)$", fontdict=None)
    axs[n_subplots-1].set_xlim(2e-3, r_lim_kpc*1e3)
    axs[3].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
    for i in range(n_subplots):
        axs[i].axvline(x=rvir_pc, color='grey', linestyle='dashed', lw=linewidth, alpha=1, label=r"$r_{200/vir}$")
    axs[2].set_ylabel(r"$\rm t_{ff} \, (Myr)$", fontdict=None)
    #axs[2].set_ylim(200, rp2[("gas", "blackhole_freefall_timescale")][rp2.used].to("Myr").d.max()+100)
    axs[1].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
    axs[0].set_ylabel(r"$\rm M_{encl} \, (M_{\odot})$", fontdict=None)
    axs[0].axhline(y=M200, color='grey', linestyle='dotted', lw=linewidth, alpha=1, label=r"$M_{200/vir}$")
    axs[0].legend(loc="upper left", fontsize=fontsize-4, ncol=1)  # upper/lower
    axs[0].set_title("Gas properties at time of BH formation", fontdict=None)
    axs[3].axhline(y=360*critical_density(ds), color='grey', linestyle='dotted', lw=linewidth, alpha=1, label=r"$n_{200}$")

    # save plot as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.6, 6.2)
    plot_name = 'radial_profile_ics_halo_{}kpc.pdf'.format(r_lim_kpc)
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + str(plot_name))
