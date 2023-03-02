"""
For analysing the nuclear disc feeding a black hole
python find_disc_attributes.py DD0130/DD0130
"""

import yt
import sys
import os
from pathlib import Path
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt


# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))
beckmann_method = 0

# naming sphere container directory
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84
sp_container_dir = root_dir[index:]
sp_container_path = Path("sphere_containers_{}/{}_sphere.h5".format(sp_container_dir, ds))

# center objects on BH
ss_pos, ss_mass, ss_age = ss_properties(ds)

# make sphere object
width = 0.0001
sp = ds.sphere(ss_pos, width)

# angular momentum vector
L = sp.quantities.angular_momentum_vector()
L_kms = L.in_units('km**2/s')

G = 6.67e-8 * yt.units.g # cgs
M = 270*1.989e+33 * yt.units.g # cgs

def _keplerian_frequency_BH(field, data):
    return np.sqrt(ds.units.newtons_constant*M / (data[("index", "radius")]**3))

# add epicyclic frequency kappa / angular frequency omega derived field
def _angular_frequency_keplerian(field, data):
    return np.abs(data["gas", "radial_velocity"]) / (2*np.pi*data["index", "radius"]) #cm/s / cm = 1/s

def _angular_frequency(field, data):
    omega = np.abs(data["gas", "velocity_cylindrical_theta"]) / (2 * np.pi * data["index", "radius"])  # cm/s / cm
    return omega

def _epicyclic_frequency_ratio(field, data):
    omega = np.abs(data["gas", "velocity_cylindrical_theta"])/(2*np.pi*data["index", "radius"]) #cm/s / cm
    return omega/data["gas", "angular_frequency_keplerian"]


ds.add_field(
    name=("gas", "angular_frequency_keplerian"),
    function=_angular_frequency_keplerian,
    sampling_type="local",
    units="1/s",
)

ds.add_field(
    name=("gas", "epicyclic_frequency_ratio"),
    function=_epicyclic_frequency_ratio,
    sampling_type="local",
    units="dimensionless",
)

ds.add_field(
    name=("gas", "angular_frequency"),
    function=_angular_frequency,
    sampling_type="local",
    units="1/s",
)

ds.add_field(
    name=("gas", "keplerian_frequency_BH"),
    function=_keplerian_frequency_BH,
    sampling_type="local",
    units="1/s",
)

# make disc out of disk function
disk = ds.disk(ss_pos, L, (10, "pc"), (0.01, "pc"))

disc_prj = yt.ProjectionPlot(ds, "z", ("gas", "H_nuclei_density"), center=ss_pos, width=(2, "pc"), data_source=disk)
slice = yt.OffAxisSlicePlot(ds, L, ("gas", "H_nuclei_density"), center=ss_pos, width=(2, "pc"), data_source=disk)
slice.save()
disc_prj.save("disc_prj")

profile = yt.create_profile(
    data_source=disk,
    bin_fields=[("index", "radius")],
    fields=[("gas", "velocity_cylindrical_theta"), ("gas", "H_nuclei_density"), ('index', 'cylindrical_z'),
            ('gas', 'temperature'), ('gas', 'sound_speed'), ('gas', 'radial_velocity'), ('gas', 'epicyclic_frequency_ratio'),
            ('gas', 'angular_frequency'), ('gas', 'velocity_spherical_theta'), ('gas', 'angular_frequency_keplerian'),
            ('index', 'radius'), ('gas', 'keplerian_frequency_BH')],
    n_bins=64,
    units=dict(#cylindrical_radius="pc",
               radius="pc", velocity_cylindrical_theta="km/s", sound_speed="km/s", velocity_spherical_theta="km/s",
               cylindrical_z="pc", radial_velocity="km/s", keplerian_frequency_BH="1/s"),
    logs=dict(cylindrical_radius=False),
    weight_field=("gas", "mass"),
    #extrema=dict(cylindrical_radius=(0, 3)),
)

profile2 = yt.create_profile(
    data_source=disk,
    bin_fields=[("index", "cylindrical_radius")],
    fields=[("gas", "H_nuclei_density")],
    n_bins=64,
    units=dict(cylindrical_radius="pc"),
    logs=dict(cylindrical_radius=False),
    weight_field=("gas", "mass"),
)

plt.plot(profile2.x[profile2.used], profile2[("gas", "H_nuclei_density")][profile2.used])
plt.yscale('log')
plt.savefig('z-dens')
print("save z-dens")

fig = plt.figure()
fig, axs = plt.subplots(6, 1, sharex=True)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"
linewidth = 1
plt.rcParams['lines.linewidth'] = linewidth

fontsize = 8
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': fontsize,
        }

# ignore invalid value error in plot_theta divide
np.seterr(invalid='ignore')

# define plots
plot_theta = axs[4].plot(profile.x.value, np.abs(profile[("gas", "velocity_cylindrical_theta")].value /
                                           profile[("gas", "sound_speed")].value))
plot_density = axs[1].plot(profile.x[profile.used], profile[("gas", "H_nuclei_density")][profile.used])
plot_z = axs[3].loglog(profile.x[profile.used], np.abs(profile[("index", "cylindrical_z")][profile.used]))
plot_temp = axs[2].loglog(profile.x[profile.used], profile[("gas", "temperature")][profile.used])
plot_vr = axs[5].loglog(profile.x[profile.used], np.abs(profile[("gas", "radial_velocity")][profile.used]))
plot_omega = axs[0].loglog(profile.x[profile.used], profile[("gas", "angular_frequency")][profile.used] /
              profile[("gas", "keplerian_frequency_BH")][profile.used])
# plot_omega = axs[5].loglog(profile.x[profile.used], profile[("gas", "epicyclic_frequency_ratio")][profile.used])
# axs[5].loglog(profile.x[profile.used], np.abs(profile[("gas", "velocity_spherical_theta")][profile.used]) /
#               (2*np.pi*profile[("index", "radius")][profile.used]))
# axs[5].loglog(profile.x[profile.used], profile[("gas", "angular_frequency_keplerian")][profile.used])
# axs[5].loglog(profile.x[profile.used], profile[("gas", "keplerian_frequency_BH")][profile.used])

axs[5].set_xlabel(r"$Radius \, (pc)$", fontdict=font)
axs[4].set_ylabel(r"$\nu_{\theta} \,/ c_s$", fontdict=font)
axs[1].set_ylabel(r"$n \, (cm^{-3})$", fontdict=font)
axs[1].set_yscale('log')
axs[3].set_ylabel(r"$Cylindrical \,Z \, (pc)$", fontdict=font)
axs[2].set_ylabel(r"$T \, (K)$", fontdict=font)
axs[5].set_ylabel(r"$\nu_r \, (km/s)$", fontdict=font)
axs[0].set_ylabel(r"$\omega / \omega_K $", fontdict=font)
axs[0].set_yscale('linear')

for i in range(6):
    axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
    axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(6, 8)
plot_name = 'disc_attributes_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created plots/" + str(plot_name))

# Produce a 2D array of uniform pixel sizes of the disc height at the maximum resolution of simulation
if beckmann_method:
    # make disc object
    rho_disc = 0  # by density projection plot inspection
    print(ss_mass)
    disc = ds.cut_region(disk, ["obj['density'].in_units('amu/cm**3') > {0}".format(rho_disc)])
    disc_height_sum = disc.sum('dy', axis="y")
    print(disc_height_sum)

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