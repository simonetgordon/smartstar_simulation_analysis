"""
Plots 6x1 radial profile of disc attributes (excluding surface density)
python -i plot_disc_attributes.py "1Bb01-1Sb01"
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


# datasets
input = sys.argv[-1]
root_dir =["/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/", 
           "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/" ]
sim = ["1B.RSb01-2", "1S.RSb01"]
dds = ["DD0144/DD0144", "DD0154/DD0154"]
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
    for k, ds in enumerate(DS):
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
                    ("gas", "omega"), ("gas", "omega_k")],
            n_bins=n_bins,
            units=dict(
                    radius="pc", velocity_cylindrical_theta="km/s", sound_speed="km/s", velocity_spherical_theta="km/s",
                    cylindrical_z="pc", radial_velocity="cm/s", keplerian_frequency_BH="1/s", tangential_velocity="km/s",
                    height="pc"),
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
        #                                           Plot All Disc Attributes
        ##########################################################################################################

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
        n_subplots = 6
        fig, axs = plt.subplots(n_subplots, 1, sharex=True)

        # ignore invalid value error in plot_theta divide
        np.seterr(invalid='ignore')

        # define plots
        l = ["1Bb01", "1Sb01"]
        plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "omega")][profile.used] /
                    profile[("gas", "omega_k")][profile.used])
        plot_density = axs[1].plot(profile.x[profile.used], profile[("gas", "H_nuclei_density")][profile.used])
        plot_temp = axs[2].loglog(profile.x[profile.used], profile[("gas", "temperature")][profile.used])
        plot_h = axs[3].loglog(profile.x[profile.used], profile[("index", "height")][profile.used])
        plot_theta = axs[4].plot(profile.x.value, np.abs(profile[("gas", "tangential_velocity")].value /
                                                        profile[("gas", "sound_speed")].value))
        plot_vr = axs[5].loglog(profile.x[profile.used], np.abs(profile[("gas", "radial_velocity")][profile.used]))


        # format plots
        axs[-1].set_xlabel("Radius (pc)", fontdict=None)
        axs[5].set_ylabel(r"$\rm \nu_r \, (km/s)$", fontdict=None)
        axs[5].set_yscale('linear')
        axs[4].set_ylabel(r"$\rm \nu_{\theta} \,/ c_s$", fontdict=None)
        axs[3].set_ylabel(r"$\rm H \, (pc)$", fontdict=None)
        axs[2].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
        axs[1].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
        axs[1].set_yscale('log')
        axs[0].set_ylabel(r"$\rm \omega / \omega_K $", fontdict=None)
        axs[0].set_yscale('linear')
        axs[0].axhline(y=1, color='grey', linestyle='dashed', alpha=1)
        axs[0].legend()
        #axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]))

        for i in range(n_subplots):
            axs[i].set_xlim([7e-3, 1e1])
            axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
            axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
            axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
            axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)


        # save plot as pdf
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.set_size_inches(6, 10)
        plot_name = 'disc_rp_' + str(input) + '.pdf'
        plt.savefig('plots/' + plot_name, bbox_inches='tight')
        print("created plots/" + str(plot_name))
