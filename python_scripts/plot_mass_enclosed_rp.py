##################################### Enclosed Mass Radial Profile  ######################################
#
# to run: python -i plot_mass_enclosed_rp.py radial_profile_mass_enclosed.pdf
##########################################################################################################

from plot_ic_attributes import *
import matplotlib.cm as cm
import matplotlib.ticker as ticker

def set_tick_params(n_subplots, fontsize, custom_ticks=[], xlim=[]):
    for i in range(n_subplots):
        axs[i].set_xscale('log')
        if len(xlim) > 0:
            axs[i].set_xlim(xlim)
        axs[i].tick_params(bottom=True, left=True)
        axs[i].minorticks_on()
        if len(custom_ticks) > 0:
            axs[i].set_xticks(custom_ticks)
        axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
        axs[i].tick_params(axis="y", which='minor', length=2)

if __name__ == "__main__":

    ############################################# Parameters #############################################
    rvir_pc = 2167.60571289 # at z = 26.338 at l3
    M200 = 2.7491550e+05    # at z = 26.338 at l3 (same as Mvir)
    n_subplots = 4          # number of subplots
    y = sys.argv[-1]        # naming plot
    fontsize = 12           # for projection annotations
    linewidth = 2
    r_lim_kpc = 0.1         # xlim in kpc
    alpha = 0.9             # line opacity
    custom_ticks=[]
    xlim = [0.00101, 90]

    # input data
    #dds = ["DD0128/DD0128", "DD0136/DD0136", "DD0148/DD0148", "DD0159/DD0159", "DD0161/DD0161", "DD0178/DD0178"] 1B.b16
    dds = ["DD0128/DD0128", "DD0129/DD0129", "DD0130/DD0130", "DD0132/DD0132", "DD0134/DD0134", "DD0138/DD0138"]
    root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"]*len(dds)
    sim = ["1B.RSb01-2"]*len(dds)

    # Get the RGB values for the Plasma colormap ('cool', 'brg', 'rainbow', 'jet', 'turbo' also good)
    c = cm.get_cmap('plasma', len(dds)+1).colors
    ######################################################################################################

    # generate line labels
    DS = generate_DS(dds, root_dir, sim)
    ss_creation_time = DS[0].r['SmartStar', 'creation_time'].to('Myr')
    labels, DS = generate_labels(dds, root_dir, sim, s=False, tform=ss_creation_time)

    # set font format
    set_font(fontsize, linewidth)

    # set up plot
    fig = plt.figure()
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)

    # create radial profiles
    for i, ds in enumerate(DS):
        ss_pos = ss_properties(ds)[0]
        sp = ds.sphere(ss_pos, (r_lim_kpc, "kpc"))
        bulk_vel = sp.quantities.bulk_velocity()
        sp.set_field_parameter("bulk_velocity", bulk_vel)

        # create weight_field=None profile for cumulative sum rp
        rp = yt.create_profile(
            sp,
            "radius",
            [("gas", "mass"), ("deposit", "all_mass"),
            ],
            units={"radius": "pc", ("gas", "mass"): "Msun", ("deposit", "all_mass"): "Msun"},
            logs={"radius": True, ("gas", "mass"): True, ("deposit", "all_mass"):True},
            weight_field=None,
            accumulation=True # -True to reverse direction
        )

        # create weight_field=cell_mass profile for non-cumulative rp
        rp2 = yt.create_profile(
            sp,
            "radius",
            [("gas", "temperature"), ("gas", "mass"),
            ("gas", "number_density"), ("gas", "density"), ('enzo', 'Dark_Matter_Density'),
            ("gas", "radial_velocity")],
            units={"radius": "pc", ("gas", "mass"): "Msun", ('enzo', 'Dark_Matter_Density'): "g/cm**3",
                    ("gas", "radial_velocity"): "km/s"},
            logs={"radius": True, ("gas", "temperature"): True, ("gas", "number_density"): True},
            weight_field=("gas", "density")
        )

        # cumulative mass enclosed
        # axs[0].loglog(rp.x[rp.used], rp[("gas", "mass")][rp.used] + rp[("deposit", "all_mass")][rp.used],
        #             color=c[i], linestyle='solid', label="total_" + labels[i], alpha=alpha)
        axs[0].loglog(rp.x[rp.used], rp[("gas", "mass")][rp.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # number density
        axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "number_density")][rp2.used], 
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # temperature
        axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

        # radial velocity
        axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "radial_velocity")][rp2.used],
                    color=c[i], linestyle='solid', label=labels[i], alpha=alpha)



    # plot annotations
    set_tick_params(n_subplots, fontsize, custom_ticks=custom_ticks, xlim=xlim)
    axs[0].legend(loc="lower right", fontsize=fontsize-4, ncol=1) 
    axs[0].set_ylabel(r"$\rm M_{encl} \, (M_{\odot})$", fontdict=None)
    axs[1].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
    axs[2].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
    axs[3].set_ylabel(r"$\rm v_r \, (km/s)$", fontdict=None)
    axs[3].set_yscale("linear")
    axs[-1].set_xlabel(r"$\rm Radius \, (pc)$", fontdict=None)

    # save as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.6, 6.2)
    plot_name = 'radial_profile_mass_enclosed_gas_{}_{}kpc.pdf'.format(sim[0], r_lim_kpc)
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + str(plot_name))
