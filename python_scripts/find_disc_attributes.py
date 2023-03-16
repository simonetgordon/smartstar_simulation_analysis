"""
For analysing the nuclear disc feeding a black hole
python find_disc_attributes.py DD0148/DD0148
"""

import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds
from yt.utilities.math_utils import ortho_find

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))
add_fields_ds(ds)
beckmann_method = 1

# naming plot
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84


# make disc data container
def _make_disk_L(ds, center, width_pc, height_pc):
    width = width_pc*yt.units.pc
    height = height_pc*yt.units.pc
    sp = ds.sphere(center, width)
    L = sp.quantities.angular_momentum_vector()
    L /= np.sqrt((L ** 2).sum())
    disk = ds.disk(ss_pos, L, width, height)
    return disk, L


def radius(s_area):
    return np.sqrt(s_area/4*np.pi)

ss_pos, ss_mass, ss_age = ss_properties(ds)

if __name__ == "__main__":

    disc_r_pc = 10
    disc_h_pc = 1

    # make disk data container and define L
    disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    vecs = ortho_find(L)

    # for i, vec in enumerate(vecs):
    #     north = vecs[0] if i > 0 else vecs[1]
    #     p = yt.ProjectionPlot(ds, vec, ("gas", "density"), weight_field=("gas", "density"), north_vector=north, center=disk.center, width=2 * disk.radius, data_source=disk)
    #     p.set_axes_unit("pc")
    #     p.set_cmap(("gas", "density"), "turbo")
    #     p.save()
    #
    # for field in ("height", "cylindrical_radius"):
    #     p = yt.ProfilePlot(disk, ("index", field), ("gas", "density"), weight_field=("gas", "cell_mass"))
    #     p.set_unit(("index", field), "pc")
    #     p.save()
    #
    # p = yt.ProfilePlot(disk, ("index", "radius"), ("index", "height"), weight_field=("gas", "cell_mass"))
    # p.set_unit(("index", "height"), "pc")
    # p.set_unit(("index", "radius"), "pc")
    # p.save()

    # plot disc
    #disc_prj = yt.ProjectionPlot(ds, "z", ("gas", "H_nuclei_density"), center=ss_pos, width=(2, "pc"), data_source=disk)
    #slice = yt.OffAxisSlicePlot(ds, L, ("gas", "H_nuclei_density"), center=ss_pos, width=(2, "pc"), data_source=disk)
    #slice.save()
    #disc_prj.save()

    # Make surface density plot using FRB
    w = disc_r_pc*2 # pc
    dx = 1.229791e-02  # pc
    frb_resolution = int(w/dx)
    proj = ds.proj(("gas", "H_nuclei_density"), "x", center=ss_pos, ds=ds, data_source=disk)
    disc_frb = proj.to_frb((disc_r_pc*2, "pc"), resolution=(int(frb_resolution)), center=ss_pos)
    disc_frb_ds = disc_frb.export_dataset(fields=[("gas", "H_nuclei_density"), ("index", "radius"), ("index", "x")])

    # use frb to calculate radius
    y = disc_frb_ds.r[("gas", "H_nuclei_density")]
    #y = y[int(y.shape[0] / 2):]
    # x = np.linspace(disc_frb_ds.domain_width.in_units('pc') / frb_resolution, disc_frb_ds.domain_width.in_units('pc')/2,
    #                 y.shape[0])
    #x = np.logspace(-1.91009, 1, y.shape[0])
    # x = x[:, 0]

    xleft = np.abs(disc_frb_ds.domain_left_edge[0] - ss_pos).to('pc')
    xright = np.abs(disc_frb_ds.domain_right_edge[0] - ss_pos[0]).to('pc')
    yleft = np.abs(disc_frb_ds.domain_left_edge[1] - ss_pos[1]).to('pc')
    yright = np.abs(disc_frb_ds.domain_right_edge[1] - ss_pos[1]).to('pc')

    xwidth = (disc_frb_ds.domain_right_edge[0] - disc_frb_ds.domain_left_edge[0]).to('pc')
    xs = np.linspace(xleft, xright, y.shape[0]) # xwidth/frb_resolution
    ywidth = (disc_frb_ds.domain_right_edge[1] - disc_frb_ds.domain_left_edge[1]).to('pc')
    ys = np.linspace(-yleft, yright, y.shape[0])
    r = np.sqrt(xs**2 + ys**2)

    nbins = 520
    sigma, radius = np.histogram(r, weights=y, bins=nbins)
    plt.plot(radius[:nbins], sigma)
    plt.ylabel('$\Sigma \, (cm^{-2})$')
    plt.xlabel('$R \, (pc)$')
    plt.yscale('log')
    plt.xscale('log')
    plot_name = 'disc_sigma_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.pdf'
    plt.savefig('plots/' + plot_name, dpi=100)
    print("created plots/" + str(plot_name))

    create_profiles = 1
    n_bins = 120
    if create_profiles:
        profile = yt.create_profile(
            data_source=disk,
            bin_fields=[("index", "radius")],
            fields=[("gas", "velocity_cylindrical_theta"), ("gas", "H_nuclei_density"), ('index', 'cylindrical_z'),
                    ('gas', 'temperature'), ('gas', 'sound_speed'), ('gas', 'radial_velocity'),
                    ('gas', 'angular_frequency'), ('gas', 'velocity_spherical_theta'), ('gas', 'angular_frequency_keplerian'),
                    ('index', 'radius'), ('gas', 'keplerian_frequency_BH'), ("gas", "tangential_velocity"), ("index", "height"),
                    ("gas", "epicycle_frequency_k"), ("gas", "angular_frequency_square"), ('gas', 'velocity_spherical_phi')],
            n_bins=64,
            units=dict(
                       radius="pc", velocity_cylindrical_theta="km/s", sound_speed="km/s", velocity_spherical_theta="km/s",
                       cylindrical_z="pc", radial_velocity="km/s", keplerian_frequency_BH="1/s", tangential_velocity="km/s",
                       height="pc"),
            logs=dict(cylindrical_radius=False),
            weight_field=("gas", "cell_mass"),
        )

        profile2 = yt.create_profile(
            data_source=disk,
            bin_fields=[("index", "cylindrical_radius")],
            fields=[("gas", "H_nuclei_density")],
            n_bins=64,
            units=dict(cylindrical_radius="pc"),
            logs=dict(cylindrical_radius=False),
            weight_field=("gas", "cell_mass"),
        )

        # plot cylindrical Z vs density
        plt.plot(profile2.x[profile2.used], profile2[("gas", "H_nuclei_density")][profile2.used])
        plt.ylabel('n ($cm^{3}$)')
        plt.xlabel('Cylindrical Z (pc)')
        plt.yscale('log')
        plot_name = 'disc_z_dens_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.png'
        plt.savefig('plots/' + plot_name, dpi=100)
        print("created plots/" + str(plot_name))

        # make Toomre Q profile
        m_p = 1.67262192e-24 * yt.units.g # g
        G = 6.67e-8 * (yt.units.cm ** 3) / (yt.units.g * yt.units.s ** 2)  # cgs
        num = profile[("gas", "sound_speed")].to('cm/s') * profile[("gas", "epicycle_frequency_k")]
        sigma, radius = np.histogram(r, weights=y, bins=n_bins)
        denom = np.pi * G * sigma[:nbins] * m_p / yt.units.cm**2

        fig = plt.figure()
        fig, axs = plt.subplots(7, 1, sharex=True)
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["mathtext.default"] = "regular"
        linewidth = 1
        plt.rcParams['lines.linewidth'] = linewidth

        fontsize = 8
        font = {'family': 'serif',
                'weight': 'normal',
                'size': fontsize,
                }

        # ignore invalid value error in plot_theta divide
        np.seterr(invalid='ignore')

        # define plots
        plot_theta = axs[4].plot(profile.x.value, np.abs(profile[("gas", "tangential_velocity")].value /
                                                   profile[("gas", "sound_speed")].value))
        plot_density = axs[1].plot(profile.x[profile.used], profile[("gas", "H_nuclei_density")][profile.used])
        plot_temp = axs[2].loglog(profile.x[profile.used], profile[("gas", "temperature")][profile.used])
        plot_vr = axs[5].loglog(profile.x[profile.used], np.abs(profile[("gas", "radial_velocity")][profile.used]))
        plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "angular_frequency")][profile.used] /
                      profile[("gas", "keplerian_frequency_BH")][profile.used])
        plot_sigma = axs[6].loglog(radius[:n_bins], sigma)
        plot_h = axs[3].loglog(profile.x[profile.used], profile[("index", "height")][profile.used])

        axs[6].set_xlabel(r"$Radius \, (pc)$", fontdict=font)
        axs[4].set_ylabel(r"$\nu_{\theta} \,/ c_s$", fontdict=font)
        axs[1].set_ylabel(r"$n \, (cm^{-3})$", fontdict=font)
        axs[1].set_yscale('log')
        axs[3].set_ylabel(r"$H \, (pc)$", fontdict=font)
        axs[2].set_ylabel(r"$T \, (K)$", fontdict=font)
        axs[5].set_ylabel(r"$\nu_r \, (km/s)$", fontdict=font)
        axs[6].set_ylabel(r"$\Sigma \, (cm^{-2})$", fontdict=font)
        axs[6].set_ylim([8e24, 8e26])
        axs[0].set_ylabel(r"$\omega / \omega_K $", fontdict=font)
        axs[0].set_yscale('linear')
        axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]),
                         fontproperties=font)

        for i in range(7):
            axs[i].set_xlim([7e-3, 1e1])
            axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
            axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
            axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
            axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)

        # save plot as pdf
        fig = plt.gcf()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.set_size_inches(6, 10)
        #fig.suptitle("2 Myr " + str(root_dir[index:]), fontproperties=font)
        plot_name = 'disc_attributes_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.pdf'
        fig.savefig('plots/' + plot_name, dpi=100)
        print("created plots/" + str(plot_name))

    # Produce a 2D array of uniform pixel sizes of the disc height at the maximum resolution of simulation
    beckmann_method = 0
    if beckmann_method:
        # make disc object
        rho_disc = 0  # by density projection plot inspection
        print(ss_mass)
        ad = ds.all_data()
        disc = ad.cut_region(disk, ["obj['density'].in_units('amu/cm**3') > {0}".format(rho_disc)])
        disc_height_sum = disc.sum(('index', 'z'), axis=('index', 'y'))
        print("disc height sum: ", disc_height_sum)

        sphere_pc = (width * ds.length_unit.in_units("code_length")).in_units("pc")
        dx = 7.687095e-04 # pc
        print("sphere_pc = ", sphere_pc)
        frb_resolution = int(sphere_pc/dx)
        print("frb_res = ", frb_resolution)
        disc_frb = disc_height_sum.to_frb(width=(2*sphere_pc, 'pc'), resolution=frb_resolution, center=ss_pos)
        height_data = disc_frb['dy'].in_units('pc')
        print(height_data > 0)


    # Used a fixed resolution buffer to grid the height data onto something I could work with. Here “pc” is the total
    # size of my sphere in pc and dx is the smallest cell size in the simulation. Sink.pos is the center of my sphere
    # frb_resolution=int(pc/dx)
    # disc_frb=disc_height_sum.to_frb(width=(2*pc,'pc'),resolution=frb_resolution,center=sink.pos)
    # height_data=disc_frb['dz'].in_units('pc')
    # If I recall correctly this should give you a 2D array of uniform pixel sizes of the disc height at the maximum
    # resolution of your simulation. You could look at it with imshow if you wanted but to make the radial profiles
    # I simply binned this data in concentric radial bins centred on the array (making sure I rescaled the x-axis from
    # “pixel count” to “pc” using the frb_resolution computed above.