"""
Plots 6x1 radial profile of disc attributes (excluding surface density)
python -i plot_disc_attributes.py "1Bb01-1Sb01-x3" 
"""
import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds
from yt.utilities.math_utils import ortho_find
import matplotlib as mpl
import pandas as pd
from matplotlib import rc
from find_disc_attributes import _make_disk_L
from plot_multi_projections import tidy_data_labels
from plot_cooling_rates import extract_colors


def orbital_velocity(ds, disk):
    G = 6.67e-8 * (yt.units.cm ** 3)/(yt.units.g*yt.units.s**2) # cgs
    return np.sqrt(G * ds.r['SmartStar', 'particle_mass'].to('g') / disk['index', 'radius'].to('cm'))


def radial_profile(field, disk, n_bins, cell_width_pc):
    bins = np.logspace(np.log10(cell_width_pc), np.log10(disk["index", "radius"].to('pc').max()), n_bins+1)
    counts_r, r_bin_edges = np.histogram(disk["index", "radius"].to('pc'), bins=bins)
    y, radius = np.histogram(disk["index", "radius"].to('pc'), weights=field, bins=bins)
    y = np.nan_to_num(y)
    y = y/counts_r
    y = np.nan_to_num(y)
    y_no_zeros = y[y > 0]
    r_no_zeros = radius[:-1][y > 0]
    return r_no_zeros, y_no_zeros


if __name__ == "__main__":

    ###################################### Parameters ######################################

    x = "s1-270msun-"        # simulation type
    y = sys.argv[-1]        # naming plot
    xlim = 1             # Myrs
    time_cutoff = xlim      # Myrs
    i_start = 0             # start index
    include_beckmann_data = False
    alpha = 0.9             # transparency of lines
    num_subplots = 7        # number of subplots
    smooth_simulations = 4  # number of simulations to smooth (starting from last sim)
    window = 10              # window size to average over
    rtol=1e-6   
    atol=2e-4               # 1e-4 for 1S.m, 
    title = '1B.b Local Properties Comparison'

    # font
    fontsize = 12 # for projection annotations
    linewidth = 2

    # datasets
    input = sys.argv[-1]
    # root_dir =["/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/", 
    #            "/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/",
    #            "/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/" ]
    # sim = ["1B.RSb01-2", "1S.RSb01", "1B.RSm01"]
    # dds = ["DD0138/DD0138", "DD0140/DD0140", "DD0138/DD0138"]
    # at t = 0.5 kyr in 1S group
    # root_dir =["/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/",
    #            "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/",
    #            "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/", 
    #            "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/"]
    # sim = ["1S.RSb04", "1S.RSm04", "1S.RSbf4", "1S.RSmf4"]
    ## at t = 1 Myr in 1B.b group
    root_dir = ["/cephfs/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/", 
                "/cephfs/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                "/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/",
                "/cephfs/sgordon/disk14/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/",
                "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/",
                "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/"]
               
    sim = ["1B.RSm01", "1B.RSm04-2", "1B.m16-4dx", \
           "2B.RSm01", "2B.RSm04", "2B.m08-4dx"]
    # dds = ["DD0138/DD0138", "DD0138/DD0138", "DD0166/DD0166", \ # 1B.b, 2B.b
    #        "DD0208/DD0208", "DD0208/DD0208", "DD0368/DD0368"]
    # 1B.m, 2B.m
    dds = ["DD0138/DD0138", "DD0138/DD0138", "DD0202/DD0202", \
           "DD0208/DD0208", "DD0370/DD0370", "DD0291/DD0291"]

    ########################################################################################

    # make and tidy data labels
    labels = []
    DS = []
    for i, dd in enumerate(dds):
        ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
        add_fields_ds(ds)
        j = []
        formation_time = 124.76 if sim[i].startswith('1') else 195.59
        label = str(sim[i]) + "-" + str(float(ds.current_time.to('Myr')) - formation_time)[:3] + "Myr"
        DS.append(ds)
        labels.append(label)

    labels = tidy_data_labels(labels)

    # font settings
    mpl.rcParams['pgf.texsystem'] = 'pdflatex'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    #rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    mpl.rcParams['text.usetex'] = True
    plt.rcParams['lines.linewidth'] = linewidth

    # create figure
    n_subplots = num_subplots
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)


    """""""""""""""""""""
    1) Plot Beckmann data
    """""""""""""""""""""

    if include_beckmann_data:
        # make array of Beckmann data
        beck_data_fp = '/cephfs/sgordon/smartstar_simulation_analysis/python_scripts/beckmann-data/R_128_Figure14/'
        csv_files = ['omega.csv', 'number_density.csv', 'temperature.csv', 'height.csv', 'theta_c_ratio.csv', 'radial_velocity.csv']
        beck_data_arr = []
        l_beck = "R_128-1.6Myr" # label
        alpha = 0.9      # transparency of lines

        # Loop through each CSV file
        for k, file_name in enumerate(csv_files):
            # Create the full file path
            file_path = beck_data_fp + file_name
            
            # Load CSV file into DataFrame
            df = pd.read_csv(file_path)
            
            # Extract NumPy array and append to the list
            beck_data_arr.append(df.to_numpy())
            beck_data = df.to_numpy()

            if k == 5:
                beck_data[:, 1] = beck_data[:, 1]/1e5

            axs[k].plot(beck_data[:,0], beck_data[:, 1], color="darkblue", linestyle='solid', label=l_beck, alpha=alpha)

            if k == 3:
                axs[6].plot(beck_data[:,0], beck_data[:, 1]/beck_data[:,0], color="darkblue", linestyle='solid', label=l_beck, alpha=alpha)


    """""""""""""""""
    2) Plot my data
    """""""""""""""""

    for k, ds in enumerate(DS):

        #cell_width_pc = [0.012, 0.003, 0.012]
        #cell_width_pc = [0.00078, 0.00078, 0.012, 0.012] # BHL vs MF
        cell_width_pc = [0.0124, 0.00309, 0.00077, 0.0083, 0.00207, 0.000518 ] # 1B.b, 2B.b res++

        # grab bh particle properties
        ss_pos, ss_mass, ss_age = ss_properties(ds)

        # make small disk data container to define L
        disc_r_pc = cell_width_pc[k]*10*yt.units.pc
        disc_h_pc = 0.02*yt.units.pc
        disk_small, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

        # make larger disc from this L for plotting
        disc_r_pc = 10*yt.units.pc
        disc_h_pc = 0.1*yt.units.pc
        disk = ds.disk(ss_pos, L, disc_r_pc, disc_h_pc)

        ##########################################################################################################
        #                                           Create Profiles
        ##########################################################################################################

        n_bins = 64
        profile = yt.create_profile(
            data_source=disk,
            bin_fields=[("index", "radius")],
            fields=[("gas", "velocity_cylindrical_theta"), ("gas", "number_density"), ('index', 'cylindrical_z'),
                    ('gas', 'temperature'), ('gas', 'sound_speed'), ('gas', 'radial_velocity'), ('gas', 'velocity_magnitude'),
                    ('gas', 'angular_frequency'), ('gas', 'velocity_spherical_theta'),
                    ('index', 'radius'), ('gas', 'keplerian_frequency_BH'), ("gas", "tangential_velocity"), ("index", "height"),
                    ('gas', 'velocity_spherical_phi'),
                    ("gas", "omega"), ("gas", "omega_k"), ('gas', 'orbital_velocity')],
            n_bins=n_bins,
            units=dict(
                    radius="pc", velocity_cylindrical_theta="km/s", sound_speed="km/s", velocity_spherical_theta="km/s",
                    cylindrical_z="pc", radial_velocity="km/s", keplerian_frequency_BH="1/s", tangential_velocity="km/s",
                    height="pc", orbital_velocity="km/s"),
            logs=dict(cylindrical_radius=False),
            weight_field=("gas", "cell_mass"),
        )

        profile2 = yt.create_profile(
        data_source=disk,
        bin_fields=[("index", "cylindrical_radius")],
        fields=[("gas", "H_nuclei_density")],
        n_bins=n_bins,
        units=dict(cylindrical_radius="pc"),
        logs=dict(cylindrical_radius=False),
        weight_field=("gas", "cell_mass"),
        )

        ##########################################################################################################
        #                                           Calculate disc height
        ##########################################################################################################

        # querying height (abs(cylindrical_z)) on a cut region of max small disc density/5
        cut = disk_small[('gas', 'number_density')].d.max()/100
        #cut = 1e7
        cr = ds.cut_region(disk, ["obj[('gas', 'number_density')] > {}".format(cut)])
        h_disc = cr[("index", "height")].to('pc')
        r_disc = cr[("index", "radius")].to('pc')
        r_h, h = radial_profile(h_disc, cr, n_bins, cell_width_pc[k])

        ##########################################################################################################
        #                                           Plot All Disc Attributes
        ##########################################################################################################


        # ignore invalid value error in plot_theta divide
        np.seterr(invalid='ignore')

        # define plots
        #c = ['blueviolet', 'turquoise', 'limegreen'] # for beck vs me rp
        #c = [ '#C51400', '#114993', '#E27396', '#68BFDE'] # for BHL vs MF rp red, blue, pink, light blue
        # colors
        n = int(len(DS)/2)
        c_s1 = extract_colors('viridis', n, portion="middle")
        c_s2 = extract_colors('magma', n, portion="middle", start=0.25, end=0.8)
        c = np.concatenate((c_s1, c_s2))

        plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "omega")][profile.used] /
                    profile[("gas", "omega_k")][profile.used], color=c[k], label=labels[k])
        plot_density = axs[1].plot(profile.x[profile.used], profile[("gas", "number_density")][profile.used], color=c[k], label=labels[k])
        plot_temp = axs[2].loglog(profile.x[profile.used], profile[("gas", "temperature")][profile.used], color=c[k], label=labels[k])
        plot_h = axs[3].loglog(r_h, h, color=c[k], label=labels[k])
        plot_theta = axs[4].plot(profile.x.value, np.abs(profile[("gas", "tangential_velocity")].value/profile[("gas", "sound_speed")].value), 
                                 color=c[k], label=labels[k])
        plot_vr = axs[5].plot(profile.x[profile.used], profile[("gas", "radial_velocity")][profile.used], color=c[k], label=labels[k])
        r, vorb = radial_profile(orbital_velocity(ds, disk).to('km/s'), disk, n_bins, cell_width_pc[k])
        #plot_vorb = axs[6].plot(r, vorb, color=c[k])
        plot_hratio = axs[6].plot(r_h, h/r_h, color=c[k], label=labels[k])


        ##########################################################################################################
        #                                               Format Plots
        ##########################################################################################################

        # format plots
        axs[-1].set_xlabel("Radius (pc)", fontdict=None)
        #axs[6].set_ylabel(r"$\rm \nu_{orbit} \, (km/s)$", fontdict=None)
        axs[6].set_ylabel(r"$\rm H/r $", fontdict=None)
        axs[6].axhline(y=1, color='grey', linestyle='dashed', alpha=1)
        axs[6].set_yscale('log')
        axs[5].set_ylabel(r"$\rm \nu_r \, (km/s)$", fontdict=None)
        axs[5].set_yscale('linear')
        axs[5].set_ylim([-100,6])
        axs[4].set_ylabel(r"$\rm \nu_{\theta}/c_s$", fontdict=None)
        axs[4].set_yscale('log')
        axs[3].set_ylabel(r"$\rm H \, (pc)$", fontdict=None)
        axs[2].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
        axs[2].set_ylim([50,2.1e5])
        axs[1].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
        axs[1].set_yscale('log')
        axs[0].set_ylabel(r"$\rm \omega / \omega_K $", fontdict=None)
        axs[0].set_yscale('linear')
        axs[0].axhline(y=1, color='grey', linestyle='dashed', alpha=1)
        axs[0].legend(loc="upper left", fontsize=fontsize-1, ncol = 2)
        #axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]))

        for i in range(n_subplots):
            axs[i].set_xlim([4e-4, 1e1])
            axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
            axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
            axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
            axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)


    # save plot as pdf
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(5.4, 10)
    plot_name = 'disc_rp_' + str(input) + '.pdf'
    plt.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/" + str(plot_name))
