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
from scipy import stats
from matplotlib import rc
from find_disc_attributes import _make_disk_L
from plot_multi_projections import tidy_data_labels

def orbital_velocity(ds, disk):
    #G = ds.parameters['GravitationalConstant']
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


# datasets
input = sys.argv[-1]
root_dir =["/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/", 
           "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/",
            "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/" ]
sim = ["1B.RSb01-2", "1S.RSb01", "1B.RSm01"]
dds = ["DD0144/DD0144", "DD0154/DD0154", "DD0144/DD0144"]
labels = []
DS = []
for i, dd in enumerate(dds):
    ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
    add_fields_ds(ds)
    j = []
    label = str(sim[i]) + "-" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"
    DS.append(ds)
    labels.append(label)

labels = tidy_data_labels(labels)

# naming plot
seed = int(root_dir[0][43:44])
if seed == 1:
    index = 82
elif seed == 2:
    index = 84


if __name__ == "__main__":

    # font settings
    fontsize = 12 # for projection annotations
    linewidth = 2
    mpl.rcParams['pgf.texsystem'] = 'pdflatex'
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    #rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    mpl.rcParams['text.usetex'] = True
    plt.rcParams['lines.linewidth'] = linewidth

    # create figure
    n_subplots = 7
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)
    

    for k, ds in enumerate(DS):

        cell_width_pc = [0.012, 0.003, 0.012]

        # grab bh particle properties
        ss_pos, ss_mass, ss_age = ss_properties(ds)

        # make disk data container and define L
        disc_r_pc = 10*yt.units.pc
        disc_h_pc = 1*yt.units.pc
        disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

        ##########################################################################################################
        #                                           Create Profiles
        ##########################################################################################################

        n_bins = 64
        profile = yt.create_profile(
            data_source=disk,
            bin_fields=[("index", "radius")],
            fields=[("gas", "velocity_cylindrical_theta"), ("gas", "H_nuclei_density"), ('index', 'cylindrical_z'),
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
        #                                           Create Height Plot
        ##########################################################################################################

        cr = ds.cut_region(disk, ["obj[('gas', 'number_density')] > 4e6"])
        h_disc = cr[("index", "height")].to('pc')
        r_disc = cr[("index", "radius")].to('pc')
        r_h, h = radial_profile(h_disc, cr, n_bins, cell_width_pc[k])

        frb_resolution=int(disk.radius.to('pc')/cell_width_pc[k])
        disc_height_sum=cr.sum('dz',axis="z")
        disc_frb=disc_height_sum.to_frb(width=2*disk.radius.to('pc'), resolution=frb_resolution, center=disk.center)
        height_data=disc_frb['dz'].in_units('pc')

        # get radius
        bds = disc_frb.bounds
        shape = disc_frb.buff_size
        dx = (bds[1] - bds[0]) / shape[0]
        dy = (bds[3] - bds[2]) / shape[1]
        px, py = np.meshgrid(np.arange((bds[0] + dx / 2), (bds[1] + dx / 2), dx),
                            np.arange((bds[2] + dy / 2), (bds[3] + dy / 2), (dy)))
        pr = ds.arr(np.sqrt(px ** 2 + py ** 2), "code_length").to('pc')


        ##########################################################################################################
        #                                           Plot All Disc Attributes
        ##########################################################################################################


        # ignore invalid value error in plot_theta divide
        np.seterr(invalid='ignore')

        # define plots
        c = ['blueviolet', 'turquoise', 'limegreen']
        plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "omega")][profile.used] /
                    profile[("gas", "omega_k")][profile.used], color=c[k], label=labels[k])
        plot_density = axs[1].plot(profile.x[profile.used], profile[("gas", "H_nuclei_density")][profile.used], color=c[k])
        plot_temp = axs[2].loglog(profile.x[profile.used], profile[("gas", "temperature")][profile.used], color=c[k])
        plot_h = axs[3].loglog(r_h, h, color=c[k])
        plot_theta = axs[4].plot(profile.x.value, np.abs(profile[("gas", "tangential_velocity")].value), color=c[k])
        plot_vr = axs[5].plot(profile.x[profile.used], profile[("gas", "radial_velocity")][profile.used], color=c[k])
        r, vorb = radial_profile(orbital_velocity(ds, disk).to('km/s'), disk, n_bins, cell_width_pc[k])
        plot_vorb = axs[6].plot(r, vorb, color=c[k])


        # format plots
        axs[-1].set_xlabel("Radius (pc)", fontdict=None)
        axs[6].set_ylabel(r"$\rm \nu_{orbit} \, (km/s)$", fontdict=None)
        axs[5].set_ylabel(r"$\rm \nu_r \, (km/s)$", fontdict=None)
        axs[5].set_yscale('linear')
        axs[4].set_ylabel(r"$\rm \nu_{\theta} \, (km/s)$", fontdict=None)
        axs[4].set_ylim([-1,25])
        axs[3].set_ylabel(r"$\rm H \, (pc)$", fontdict=None)
        axs[2].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
        axs[1].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
        axs[1].set_yscale('log')
        axs[0].set_ylabel(r"$\rm \omega / \omega_K $", fontdict=None)
        axs[0].set_yscale('linear')
        axs[0].axhline(y=1, color='grey', linestyle='dashed', alpha=1)
        axs[0].legend(fontsize=fontsize-1)
        #axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]))

        for i in range(n_subplots):
            axs[i].set_xlim([1.5e-3, 1e1])
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
