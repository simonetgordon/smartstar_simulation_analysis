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
import matplotlib as mpl
mpl.rcParams['agg.path.chunksize'] = 10000
from scipy import stats

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


def myExpFunc(x, a, b):
    return a * np.power(x, b)

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

    for i, vec in enumerate(vecs):
        north = vecs[0] if i > 0 else vecs[1]
        p = yt.ProjectionPlot(ds, vec, ("gas", "density"), weight_field=("gas", "density"), north_vector=north,
                              center=disk.center, width=2 * disk.radius, data_source=disk)
        p.set_axes_unit("pc")
        p.set_cmap(("gas", "density"), "turbo")
        p.save("disk_{}".format(i))

    for field in ("height", "cylindrical_radius"):
        p = yt.ProfilePlot(disk, ("index", field), ("gas", "density"), weight_field=("gas", "cell_mass"))
        p.set_unit(("index", field), "pc")
        p.save()

    p = yt.ProfilePlot(disk, ("index", "radius"), ("index", "height"), weight_field=("gas", "cell_mass"))
    p.set_unit(("index", "height"), "pc")
    p.set_unit(("index", "radius"), "pc")
    p.save()

    # plot 2d hist of sigma from off-axis projection
    p = yt.ProjectionPlot(ds, L, ("gas", "H_nuclei_density"), weight_field=None, north_vector=vecs[1],
                                 center=disk.center, width=2*disk.radius, data_source=disk)
    density_frb = p.frb[("gas", "number_density")]
    bds = p.frb.bounds
    shape = p.frb.buff_size
    dx = ((bds[1] - bds[0]) / shape[0])
    dy = ((bds[3] - bds[2]) / shape[1])

    #px, py = np.mgrid[(bds[0]+dx/2):(bds[1]+dx/2):dx, (bds[2]+dy/2):(bds[3]+dy/2):(dy)]
    px, py = np.meshgrid(np.arange((bds[0] + dx / 2), (bds[1] + dx / 2), dx), np.arange((bds[2] + dy / 2), (bds[3] + dy / 2), (dy)))
    pr = ds.arr(np.sqrt(px**2 + py**2), "pc")

    hist, xedges, yedges = np.histogram2d(pr.flatten(), density_frb.flatten(), bins=64)

    fig = plt.figure()
    H = hist.T
    X, Y = np.meshgrid(xedges, yedges)
    plt.pcolor(X, Y, H, vmin=1, vmax=3000)
    plt.colorbar()
    plt.xscale('log')
    plt.yscale('log')
    plot_name = "hist2D_proj.png"
    plt.ylabel('$\Sigma \, (cm^{-2})$')
    plt.xlabel('$R \, (pc)$')
    fig.savefig('plots/' + plot_name, dpi=100)
    plt.show()

    #p.set_cmap(("gas", "H_nuclei_density"), "turbo")
    #p.save("disk_sigma_L")

    # make density projection and use frb to extract 2D data array
    w = disc_r_pc*2*yt.units.pc # pc -  0.00076892 in code units
    dx = 1.229791e-02  # pc
    frb_resolution = int(w/dx)
    proj = ds.proj(("gas", "H_nuclei_density"), "x", center=ss_pos, ds=ds, data_source=disk)
    proj_dens = proj[("gas", "H_nuclei_density")]
    disc_frb = proj.to_frb((disc_r_pc*2, "pc"), resolution=(int(frb_resolution)), center=ss_pos)
    #yt.write_image(np.log10(disc_frb[("gas", "H_nuclei_density")]), "disc_frb_proj.png")
    disc_frb_ds = disc_frb.export_dataset(fields=[("gas", "H_nuclei_density"), ("index", "radius"), ("index", "x")])

    # find proj radius
    xs = np.abs(disk.center[1] - proj["px"])
    ys = np.abs(disk.center[2] - proj["py"])
    proj_rad = ds.arr(np.sqrt(xs**2 + ys**2), "pc")

    create_profiles = 1
    n_bins = 120
    if create_profiles:
        profile = yt.create_profile(
            data_source=disk,
            bin_fields=[("index", "radius")],
            fields=[("gas", "velocity_cylindrical_theta"), ("gas", "H_nuclei_density"), ('index', 'cylindrical_z'),
                    ('gas', 'temperature'), ('gas', 'sound_speed'), ('gas', 'radial_velocity'),
                    ('gas', 'angular_frequency'), ('gas', 'velocity_spherical_theta'),
                    ('index', 'radius'), ('gas', 'keplerian_frequency_BH'), ("gas", "tangential_velocity"), ("index", "height"),
                    ('gas', 'velocity_spherical_phi'),
                    ("gas", "omega"), ("gas", "omega_k")],
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
        plot_z_dens = 0
        if plot_z_dens:
            plt.plot(profile2.x[profile2.used], profile2[("gas", "H_nuclei_density")][profile2.used])
            plt.ylabel('n ($cm^{3}$)')
            plt.xlabel('Cylindrical Z (pc)')
            plt.yscale('log')
            plot_name = 'disc_z_dens_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.png'
            plt.savefig('plots/' + plot_name, dpi=100)
            print("created plots/" + str(plot_name))

        ##########################################################################################################
        #                                           Plot Sigma
        ##########################################################################################################

        # 1) divide sigma from weighted 1d histogram by r bin counts (works)
        # define sigma_frb and pr from ProjectionPlot p
        sigma_frb = p.frb[("gas", "number_density")]
        bds = p.frb.bounds
        shape = p.frb.buff_size
        dx = (bds[1] - bds[0]) / shape[0]
        dy = (bds[3] - bds[2]) / shape[1]
        px, py = np.meshgrid(np.arange((bds[0] + dx / 2), (bds[1] + dx / 2), dx),
                             np.arange((bds[2] + dy / 2), (bds[3] + dy / 2), (dy)))
        pr = ds.arr(np.sqrt(px ** 2 + py ** 2), "code_length").to('pc') # pr.min() = unyt_quantity(0.01767767, 'pc')

        counts_r, r_bin_edges = np.histogram(pr, bins=64)
        sigma, radius = np.histogram(pr, weights=sigma_frb, bins=64)
        sigma = sigma / counts_r

        fig = plt.figure()
        plot_name = "sigma_r.png"
        plt.loglog(radius[:-1], sigma)
        plt.ylabel('$\Sigma \, (cm^{-2})$')
        plt.xlabel('$R \, (pc)$')
        fig.savefig('plots/' + plot_name, dpi=100)
        plt.show()

        # 2) plot all r, y points (no zeros) with a fitted line
        # r = pr.flatten()
        # y = sigma_frb.flatten()
        # fig = plt.figure()
        # plt.loglog(r, y)
        # slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(r), np.log10(y))
        # sigma_fitted = myExpFunc(r, 10 ** intercept, slope)
        # plt.plot(r, sigma_fitted, 'green', label="({0:.2E}*x**{1:.3f})".format(10 ** intercept, slope))
        # plt.legend()
        # plot_name = "disc_sigma_with_fit.png"
        # fig.savefig('plots/' + plot_name, dpi=100)
        # print("created plots/" + str(plot_name))
        # plt.show()


        ##########################################################################################################
        #                                           Plot Toomre Q
        ##########################################################################################################

        # make Toomre Q profile
        def beckmann_fit(r):
            return ((10 ** (25.19)) * (r ** (-0.19)))

        m_p = 1.67262192e-24 * yt.units.g # g
        G = 6.67e-8 * (yt.units.cm ** 3) / (yt.units.g * yt.units.s ** 2)  # cgs
        num = profile[("gas", "sound_speed")].to('cm/s') * profile[("gas", "omega")].to('1/s')
        denom = np.pi * G * sigma * m_p
        #denom_fit_sigma = np.pi * G * myExpFunc(profile.x.value, 10 ** intercept, slope) * m_p / yt.units.cm**2
        denom_beck_fit = np.pi * G * beckmann_fit(profile.x.value) * m_p / yt.units.cm ** 2
        denom_sigma_all = np.pi * G * sigma * m_p / yt.units.cm ** 2

        ##########################################################################################################
        #                                           Plot All Disc Attributes
        ##########################################################################################################

        fig = plt.figure()
        n_subplots = 8
        fig, axs = plt.subplots(n_subplots, 1, sharex=True)
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
        plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "omega")][profile.used] /
                      profile[("gas", "omega_k")][profile.used])
        plot_sigma = axs[6].loglog(radius[:-1], sigma)
        #plot_sigma_fitted = axs[6].plot(r, sigma_fitted)
        #plot_sigma_weighted_hist = axs[6].loglog(radius[:-1], sigma)
        #plot_sigma_rp = axs[6].loglog(proj_rad, proj_dens)
        #plot_sigma_rp = axs[6].loglog(rad)
        plot_h = axs[3].loglog(profile.x[profile.used], profile[("index", "height")][profile.used])
        plot_toomreq = axs[7].loglog(profile.x.value, num/denom_sigma_all)
        #plot_toomreq_beck_fit = axs[7].loglog(profile.x.value, num / denom_beck_fit)
        #plot_toomreq_sigma_fit = axs[7].loglog(profile.x.value, num / denom_fit_sigma)

        axs[7].set_xlabel(r"$Radius \, (pc)$", fontdict=font)
        axs[7].axhline(y=1, color='grey', linestyle='dashed', lw=linewidth, alpha=1)
        axs[7].set_ylabel("Toomre Q", fontdict=font)
        axs[6].set_ylabel(r"$\Sigma \, (cm^{-2})$", fontdict=font)
        axs[6].set_ylim(2e19, 2e25)
        axs[5].set_ylabel(r"$\nu_r \, (km/s)$", fontdict=font)
        axs[4].set_ylabel(r"$\nu_{\theta} \,/ c_s$", fontdict=font)
        axs[3].set_ylabel(r"$H \, (pc)$", fontdict=font)
        axs[2].set_ylabel(r"$T \, (K)$", fontdict=font)
        axs[1].set_ylabel(r"$n \, (cm^{-3})$", fontdict=font)
        axs[1].set_yscale('log')
        axs[0].set_ylabel(r"$\omega / \omega_K $", fontdict=font)
        axs[0].set_yscale('linear')
        axs[0].axhline(y=1, color='grey', linestyle='dashed', lw=linewidth, alpha=1)
        axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]),
                         fontproperties=font)

        for i in range(n_subplots):
            axs[i].set_xlim([7e-3, 1e1])
            axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
            axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
            axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
            axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)

        # save plot as pdf
        fig = plt.gcf()
        fig.subplots_adjust(wspace=0, hspace=0)
        fig.set_size_inches(6, 10)
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