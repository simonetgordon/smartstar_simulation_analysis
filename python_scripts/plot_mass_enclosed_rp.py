##################################### Enclosed Mass Radial Profile  ######################################
#
# to run: python -i plot_mass_enclosed_rp.py radial_profile_mass_enclosed.pdf
##########################################################################################################

from plot_ic_attributes import *
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from derived_fields import _mass_flux

if __name__ == "__main__":

    ############################################# Parameters #############################################
    rvir_pc = 79.34 # at 2167.60571289 comoving. This is at z = 26.338 at l3
    r_disc_pc = 0.15 # at 1 Myr
    M200 = 2.7491550e+05    # at z = 26.338 at l3 (same as Mvir)
    n_subplots = 6          # number of subplots
    y = sys.argv[-1]        # naming plots
    fontsize = 12           # for projection annotations
    linewidth = 2
    r_lim_kpc = 0.1         # xlim in kpc
    alpha = 0.9             # line opacity
    custom_ticks=[]
    xlim = [0.0002, 80]     # 0.003 1B.b01, 0.0002 1B.b16
    plot_title = 'Evolution of 1B.b16 disc properties over 1 Myr'

    # input data
    #dds = ["DD0128/DD0128", "DD0136/DD0136", "DD0148/DD0148", "DD0159/DD0159", "DD0161/DD0161", "DD0170/DD0170"]  # 1B.b16
    #dds = ["DD0128/DD0128", "DD0129/DD0129", "DD0130/DD0130", "DD0132/DD0132", "DD0134/DD0134", "DD0138/DD0138"] # 1B.b01
    dds = ["DD0198/DD0198","DD0218/DD0218", "DD0220/DD0220", "DD0234/DD0234", "DD0335/DD0335", "DD0534/DD0534"] # 2B.b08
    #root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/"]*len(dds)
    root_dir = ["/ceph/cephfs/sgordon/disk14/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/"]*3 + \
                ["/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/"]*3
    sim = ["2B.RSb08"]*len(dds)

    # Get the RGB values for the Plasma colormap ('coolwarm', 'brg', 'rainbow', 'jet', 'turbo' also good)
    c = cm.get_cmap('plasma', len(dds)).colors
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

        # add mass flux field
        ds.add_field(
            ("gas", "mass_flux"),
            function=_mass_flux,
            sampling_type="local",
            units="g/(cm**2 * s)"
            )

        # create weight_field=None profile for cumulative sum rp
        rp = yt.create_profile(
            sp,
            "radius",
            [("gas", "mass"), ("deposit", "all_mass"), ("gas", "mass_flux")],
            units={"radius": "pc", ("gas", "mass"): "Msun", ("deposit", "all_mass"): "Msun", 
                    ("gas", "mass_flux"): "Msun/yr"},
            logs={"radius": True, ("gas", "mass"): True, ("deposit", "all_mass"):True,},
            weight_field=None,
            accumulation=True # -True to reverse direction
        )

        # create weight_field=cell_mass profile for non-cumulative rp
        rp2 = yt.create_profile(
            sp,
            "radius",
            [("gas", "temperature"), ("gas", "mass"),
            ("gas", "number_density"), ("gas", "density"), ('enzo', 'Dark_Matter_Density'),
            ("gas", "radial_velocity"), ("gas", "velocity_cylindrical_theta"), ("gas", "sound_speed"),
            ("gas", "mass_flux")],
            units={"radius": "pc", ("gas", "mass"): "Msun", ('enzo', 'Dark_Matter_Density'): "g/cm**3",
                    ("gas", "radial_velocity"): "km/s", ("gas", "mass_flux"): "Msun/yr"},
            logs={"radius": True, ("gas", "temperature"): True, ("gas", "number_density"): True},
            weight_field=("gas", "density")
        )

        # cumulative mass enclosed
        # axs[0].loglog(rp.x[rp.used], rp[("gas", "mass")][rp.used] + rp[("deposit", "all_mass")][rp.used],
        #             color=c[i], linestyle='solid', label="total_" + labels[i], alpha=alpha)
        axs[0].loglog(rp.x[rp.used], rp[("gas", "mass")][rp.used],
                    color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # # number density
        axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "number_density")][rp2.used], 
                    color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # temperature
        axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                    color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # radial velocity
        axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "radial_velocity")][rp2.used],
                    color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # circular `mach number'
        axs[4].loglog(rp2.x[rp2.used], rp2[("gas", "velocity_cylindrical_theta")][rp2.used] / 
                    rp2[("gas", "sound_speed")][rp2.used],
                    color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        #  # shell mass 
        # axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "mass")][rp2.used],
        #         color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)
        
        #  how long it would take the mass to move to r = 0 if it were to to keep its current radial velocity 
        # vel_r = np.array(rp2[("gas", "radial_velocity")].to("pc/yr")[rp2.used])
        # vel_r[vel_r >= 0] = 0 # set postive velocities to 0        
        # axs[3].loglog(rp2.x[rp2.used], rp2.x[rp2.used]/np.abs(vel_r),
        #         color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # # mass flux
        # axs[4].loglog(rp2.x[rp2.used], rp2[("gas", "mass_flux")][rp2.used],
        #         color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

        # accretion time
        axs[5].loglog(rp2.x[rp2.used], rp[("gas", "mass")][rp.used] / 
                rp2[("gas", "mass_flux")][rp2.used],
                color=c[i], linestyle='solid', label= labels[i][7:], alpha=alpha)

    # plot annotations
    set_tick_params(axs, n_subplots, fontsize, custom_ticks=custom_ticks, xlim=xlim, vline=r_disc_pc)
    #axs[0].set_title(plot_title, fontsize=fontsize-2)
    axs[0].legend(loc="lower right", fontsize=fontsize-4, ncol=1)

    axs[0].set_ylabel(r"$\rm M_{encl} \, (M_{\odot})$", fontdict=None)
    axs[0].set_ylim([5e-7, 1e5])
    axs[1].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
    axs[1].set_ylim([2e-1, 2e10])
    axs[2].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
    axs[2].set_ylim([2e1, 1.1e5])
    #axs[2].set_ylabel(r"$\rm M(r) \, M_\odot$", fontdict=None)
    axs[3].set_ylabel(r"$\rm v_r \, (km/s)$", fontdict=None)
    axs[3].set_yscale("linear")
    axs[3].set_ylim([-80, 20])
    #axs[3].set_ylabel(r"$\rm r/v_r \, (yrs)$", fontdict=None)
    axs[4].set_ylabel(r"$\rm \nu_\theta/c_s$", fontdict=None)
    #axs[4].set_ylabel(r"$\rm Mass \, Flux \, (M_\odot/yr)$", fontdict=None)
    axs[4].set_yscale("linear")
    axs[4].set_ylim([-25, 19])

    axs[5].set_ylabel(r"$\rm M/\dot{M} \, (yrs)$", fontdict=None)
    axs[5].set_ylim([1, 5e8])
    #axs[5].set_ylabel(r"$\rm Accretion \, Time \, (yrs)$", fontdict=None)
    axs[-1].set_xlabel(r"$\rm Radius \, (pc)$", fontdict=None)

    # save as pdf
    fig = plt.gcf()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(4.4, 8.2)
    plot_name = 'radial_profile_mass_enclosed_gas_{}_{}kpc_x6.pdf'.format(sim[0], r_lim_kpc)
    fig.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + str(plot_name))
