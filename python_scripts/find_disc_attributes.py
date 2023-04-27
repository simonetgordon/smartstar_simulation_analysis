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
from matplotlib import rc

# font settings
fontsize = 12 # for projection annotations
linewidth = 2
plt.rcParams['font.size'] = fontsize
plt.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)

plt.rcParams["mathtext.default"] = "regular"
plt.rcParams['lines.linewidth'] = linewidth

# make disc data container
def _make_disk_L(ds, center, width, height):
    sp = ds.sphere(center, width)
    L = sp.quantities.angular_momentum_vector()
    L /= np.sqrt((L ** 2).sum())
    disk = ds.disk(center, L, width, height)
    return disk, L


def radius(s_area):
    return np.sqrt(s_area/4*np.pi)


def myExpFunc(x, a, b):
    return a * np.power(x, b)

if __name__ == "__main__":

    # set by user
    root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm01"
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

    # grab bh particle properties
    ss_pos, ss_mass, ss_age = ss_properties(ds)

    # make disk data container and define L
    disc_r_pc = 10*yt.units.pc
    disc_h_pc = 1*yt.units.pc
    disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    vecs = ortho_find(L)

    # get surface density and (sigma) and radius arrays from off-axis density projection
    p = yt.ProjectionPlot(ds, L, ("gas", "H_nuclei_density"), weight_field=None, north_vector=vecs[1],
                          center=disk.center, width=2*disk.radius, data_source=disk)

    # define frb resolution by (width of domain)/(minimum cell width in simulation)
    cell_width_pc = 1.229791e-02*yt.units.pc # level 15
    frb_resolution = int(disc_r_pc*2/cell_width_pc)

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
    pr = ds.arr(np.sqrt(px ** 2 + py ** 2), "code_length").to('pc')  # pr.min() = unyt_quantity(0.01767767, 'pc')

    # make radial profile from histogram
    bins = np.logspace(np.log10(cell_width_pc), np.log10(disk.radius.to('pc')), 65)
    counts_r, r_bin_edges = np.histogram(pr, bins=bins)
    sigma, radius = np.histogram(pr, weights=sigma_frb, bins=bins)
    print("r_bin_edges == radius: ", (r_bin_edges == radius))
    sigma = np.nan_to_num(sigma)
    sigma = sigma / counts_r

    # make sigma vs r plot
    fig = plt.figure()
    plot_name = "sigma_r.png"
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(radius[:-1], sigma)
    plt.ylabel('$\Sigma \, (cm^{-2})$')
    plt.xlabel('$R \, (pc)$')
    plt.xlim(0.001, 10)
    fig.savefig('plots/' + plot_name, dpi=100)
    print("created plots/" + str(plot_name))
    plt.show()


    # 2) make 2D histogram of all points
    hist, xedges, yedges = np.histogram2d(pr.ravel(), sigma_frb.ravel(), bins=512)
    fig = plt.figure()
    H = hist.T
    X, Y = np.meshgrid(xedges, yedges)
    #plt.pcolor(X, Y, H, vmin=1, vmax=200)
    #plt.colorbar()
    plt.scatter(pr, sigma_frb)
    plt.xscale('log')
    plt.yscale('log')
    plot_name = "hist2D_all_lines_2.png"
    #plot_name = "hist2D_proj.png"
    plt.ylabel('$\Sigma \, (cm^{-2})$')
    plt.xlabel('$R \, (pc)$')
    plt.ylim(1e18, 1e25)
    fig.savefig('plots/' + plot_name, dpi=100)
    print("created plots/" + str(plot_name))
    plt.show()



    # 3) plot all r, y points (no zeros) with a fitted line
    r = pr.flatten()
    y = sigma_frb.flatten()
    y_no_zeros = y[y > 0]
    r_no_zeros = r[y > 0]
    fig = plt.figure()
    plt.scatter(r_no_zeros, y_no_zeros)
    plt.xscale('log')
    plt.yscale('log')

    # fit all data - green
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(r_no_zeros), np.log10(y_no_zeros))
    sigma_fitted = myExpFunc(r_no_zeros, 10 ** intercept, slope)
    plt.plot(r_no_zeros, sigma_fitted, 'green', label="({0:.2E}*x**{1:.3f})".format(10 ** intercept, slope))

    # fit for r < 0.5 pc - orange
    cut = 1
    sub_pc_r = np.log10(r_no_zeros[r_no_zeros < cut])
    sub_pc_sigma = np.log10(y_no_zeros[r_no_zeros < cut])
    slope, intercept, r_value, p_value, std_err = stats.linregress(sub_pc_r, sub_pc_sigma)
    sigma_fitted_sub_pc = myExpFunc(r_no_zeros[r_no_zeros < cut], 10 ** intercept, slope)
    plt.plot(r_no_zeros[r_no_zeros < cut], sigma_fitted_sub_pc, 'orange',
             label="({0:.2E}*x**{1:.3f})".format(10 ** intercept, slope))

    # fit for r < 0.2 pc - purple
    cut = 0.2
    sub_pc_r = np.log10(r_no_zeros[r_no_zeros < cut])
    sub_pc_sigma = np.log10(y_no_zeros[r_no_zeros < cut])
    slope, intercept, r_value, p_value, std_err = stats.linregress(sub_pc_r, sub_pc_sigma)
    sigma_fitted_sub_pc = myExpFunc(r_no_zeros[r_no_zeros < cut], 10 ** intercept, slope)
    plt.plot(r_no_zeros[r_no_zeros < cut], sigma_fitted_sub_pc, 'purple',
             label="({0:.2E}*x**{1:.3f})".format(10 ** intercept, slope))

    plt.legend()
    plot_name = "disc_sigma_with_fit.png"
    fig.savefig('plots/' + plot_name, dpi=100)
    print("created plots/" + str(plot_name))
    plt.show()

    ##########################################################################################################
    #                                           Create Profiles
    ##########################################################################################################

    create_profiles = 1
    n_bins = 120
    profile = yt.create_profile(
        data_source=disk,
        bin_fields=[("index", "radius")],
        fields=[("gas", "velocity_cylindrical_theta"), ("gas", "H_nuclei_density"), ('index', 'cylindrical_z'),
                ('gas', 'temperature'), ('gas', 'sound_speed'), ('gas', 'radial_velocity'), ('gas', 'velocity_magnitude'),
                ('gas', 'angular_frequency'), ('gas', 'velocity_spherical_theta'),
                ('index', 'radius'), ('gas', 'keplerian_frequency_BH'), ("gas", "tangential_velocity"), ("index", "height"),
                ('gas', 'velocity_spherical_phi'),
                ("gas", "omega"), ("gas", "omega_k")],
        n_bins=64,
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
        n_bins=64,
        units=dict(cylindrical_radius="pc"),
        logs=dict(cylindrical_radius=False),
        weight_field=("gas", "cell_mass"),
    )


    ##########################################################################################################
    #                                           Plot Toomre Q
    ##########################################################################################################

    # beckmann's fit for surface density
    def beckmann_fit_sigma(r):
        return ((10 ** (25.19)) * (r ** (-0.19)))

    def my_fit_purple(r):
        return ((1.07e23) * r **(-0.963))

    def my_fit_green(r):
        return ((5.44e21) * r ** (-1.809))

    def my_fit_orange(r):
        return ((7.51e21) * r **(-2.121))


    # make Toomre Q profile
    m_p = 1.67262192e-24 * yt.units.g  # g
    G = 6.67e-8 * (yt.units.cm ** 3) / (yt.units.g * yt.units.s ** 2)  # cgs
    num = profile[("gas", "sound_speed")].to('cm/s') * profile[("gas", "velocity_magnitude")].to('cm/s')
    #denom = np.pi * G * sigma * m_p
    # denom_fit_sigma = np.pi * G * myExpFunc(profile.x.value, 10 ** intercept, slope) * m_p / yt.units.cm**2
    denom_beck_fit = np.pi * G * beckmann_fit_sigma(profile.x.to('cm')).d * m_p * profile.x.to('cm') / yt.units.cm ** 2
    denom_fit_purple = np.pi * G * my_fit_purple(profile.x.to('pc')).d * m_p * profile.x.to('cm') / yt.units.cm ** 2
    denom_fit_green = np.pi * G * my_fit_green(profile.x.to('pc')).d * m_p * profile.x.to('cm') / yt.units.cm ** 2
    denom_fit_orange = np.pi * G * my_fit_orange(profile.x.to('pc')).d * m_p * profile.x.to('cm') / yt.units.cm ** 2
    denom_sigma_all = np.pi * G * sigma * m_p * profile.x.to('cm') / yt.units.cm ** 2

    ##########################################################################################################
    #                                 Plot cylindrical Z vs Density
    ##########################################################################################################

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
    #                                           Plot All Disc Attributes
    ##########################################################################################################

    n_subplots = 8
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)

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
    plot_toomreq_green = axs[7].loglog(profile.x.value, num / denom_fit_green, color='green')
    plot_toomreq_orange = axs[7].loglog(profile.x.value, num / denom_fit_orange, color='orange')
    plot_toomreq_purple = axs[7].loglog(profile.x.value, num / denom_fit_purple, color='purple')

    #plot_toomreq_beck_fit = axs[7].loglog(profile.x.value, num / denom_beck_fit)

    # format plots
    axs[7].set_xlabel("Radius (pc)", fontdict=None)
    axs[7].axhline(y=1, color='grey', linestyle='dashed', alpha=1)
    axs[7].set_ylabel("Toomre Q", fontdict=None)
    axs[6].set_ylabel(r"$\rm \Sigma \, (cm^{-2})$", fontdict=None)
    axs[6].set_ylim(2e19, 2e25)
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
    axs[0].set_title("BH Age = " + "{:.2f}".format(ss_age[0]/1e6) + " Myr" + ", " + str(root_dir[index:]),
                     fontproperties=None)

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

    # # Produce a 2D array of uniform pixel sizes of the disc height at the maximum resolution of simulation
    # beckmann_method = 0
    # if beckmann_method:
    #     # make disc object
    #     rho_disc = 0  # by density projection plot inspection
    #     print(ss_mass)
    #     ad = ds.all_data()
    #     ad = ds.all_data()
    #     disc = ad.cut_region(disk, ["obj['density'].in_units('amu/cm**3') > {0}".format(rho_disc)])
    #     disc_height_sum = disc.sum(('index', 'z'), axis=('index', 'y'))
    #     print("disc height sum: ", disc_height_sum)
    #
    #     sphere_pc = (width * ds.length_unit.in_units("code_length")).in_units("pc")
    #     dx = 7.687095e-04 # pc
    #     print("sphere_pc = ", sphere_pc)
    #     frb_resolution = int(sphere_pc/dx)
    #     print("frb_res = ", frb_resolution)
    #     disc_frb = disc_height_sum.to_frb(width=(2*sphere_pc, 'pc'), resolution=frb_resolution, center=ss_pos)
    #     height_data = disc_frb['dy'].in_units('pc')
    #     print(height_data > 0)


    # Used a fixed resolution buffer to grid the height data onto something I could work with. Here “pc” is the total
    # size of my sphere in pc and dx is the smallest cell size in the simulation. Sink.pos is the center of my sphere
    # frb_resolution=int(pc/dx)
    # disc_frb=disc_height_sum.to_frb(width=(2*pc,'pc'),resolution=frb_resolution,center=sink.pos)
    # height_data=disc_frb['dz'].in_units('pc')
    # If I recall correctly this should give you a 2D array of uniform pixel sizes of the disc height at the maximum
    # resolution of your simulation. You could look at it with imshow if you wanted but to make the radial profiles
    # I simply binned this data in concentric radial bins centred on the array (making sure I rescaled the x-axis from
    # “pixel count” to “pc” using the frb_resolution computed above.