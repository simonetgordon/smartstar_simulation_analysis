# Before supernova plots in one 2x2 figure
# Star position currently set on position at time of formation
import yt
import ytree
import os
import numpy as np
import matplotlib.pyplot as plt

# macros
i = 1
# load data
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/BF-0.0078125"

ds1 = yt.load(os.path.join(root_dir, "DD0126/DD0126")) # t = 124.55, Time of star formation.
label1 = "DD0126"
ds2 = yt.load(os.path.join(root_dir, "DD0129/DD0129")) # t = 125.55, 1 Myr after star particle inserted
label2 = "DD0129"
# ds3 = yt.load(os.path.join(root_dir, "DD0560/DD0560")) # t = 126.55, 2 Myr after star particle inserted
# ds4 = yt.load(os.path.join(root_dir, "DD0560/DD0560")) # 127.55. 3 Myr after star particle inserted
# ds5 = yt.load(os.path.join(root_dir, "DD0560/DD0560")) # 128.41,  3.87 Myr, last dd in stellar lifetime
# ds6 = yt.load(os.path.join(root_dir, "DD0560/DD0560")) # 128.53 0.1Myr after BH is formed.

# all ds NEEDS CHANGING EACH TIME
DD = [ds1, ds2]

# define dark_matter_mass field + add to all ds
def _dm_mass(field, data):
    return data[("gas", "dark_matter_density")].in_units("g/cm**3") * data[("index", "cell_volume")].in_units("cm**3")

for ds in DD:
    ds.add_field(("gas", "dm_mass"),function=_dm_mass, units = "g", sampling_type= "cell") # need to define with units here.


""" Make Sphere of Most Massive Halo """

# Load seed1 merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds1.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 3 * virial radius
star_pos0 = [0.49053118, 0.49467636, 0.50964283]
spheres = []
for i, ds in enumerate(DD):
    sphere = ds.sphere(star_pos0, 3*r)
    spheres.append(sphere)


""" Make Profile Plots weighted by cell mass """

all_fields = ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass",
              "cell_volume", "H_p1_fraction", "H_p0_fraction"] # 'dm_mass' raises error if put here.
all_units = {"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "g/cm**3"}

# Radial profile weighted by default = 'cell_mass'
radial_profiles = []
for i, sp in enumerate(spheres):
    radial_profile = yt.create_profile(
        sp,
        "radius",
        all_fields,
        units=all_units,
        logs={"radius": True} # Set true for more bins
    )
    radial_profiles.append(radial_profile)


# delete zeros
radial_profiles_0_temp = radial_profiles[0]["temperature"][radial_profiles[0].used]
radial_profiles_0_rad = radial_profiles[0].x[radial_profiles[0].used]

radial_profiles_0_h2 = radial_profiles[0]["H_p1_fraction"][radial_profiles[0].used] # changed from H2 TO hii
radial_profiles_0_h1 = radial_profiles[0]["H_p0_fraction"][radial_profiles[0].used] # for HI fraction

if len(spheres) > 1:
    # delete zeros for temperature plot
    radial_profiles_1_temp = radial_profiles[1]["temperature"][radial_profiles[1].used]
    radial_profiles_1_rad = radial_profiles[1].x[radial_profiles[1].used]

    # delete zeros for h2 fraction plot
    radial_profiles_1_h2 = radial_profiles[1]["H_p1_fraction"][radial_profiles[1].used]
    radial_profiles_1_h1 = radial_profiles[1]["H_p0_fraction"][radial_profiles[1].used] # for HI fraction

if len(spheres) > 2:
    # delete zeros for temperature plot
    radial_profiles_2_temp= radial_profiles[2]["temperature"][radial_profiles[2].used]
    radial_profiles_2_rad = radial_profiles[2].x[radial_profiles[2].used]

    # delete zeros for h2 fraction plot
    radial_profiles_2_h2 = radial_profiles[2]["H_p1_fraction"][radial_profiles[2].used]
    radial_profiles_2_h1 = radial_profiles[2]["H_p0_fraction"][radial_profiles[2].used] # for HI fraction

if len(spheres) > 3:
    # delete zeros for temperature plot
    radial_profiles_3_temp = radial_profiles[3]["temperature"][radial_profiles[3].used]
    radial_profiles_3_rad = radial_profiles[3].x[radial_profiles[3].used]

    # delete zeros for h2 fraction plot
    radial_profiles_3_h2 = radial_profiles[3]["H_p1_fraction"][radial_profiles[3].used]
    radial_profiles_3_h1 = radial_profiles[3]["H_p0_fraction"][radial_profiles[3].used] # for HI fraction

if len(spheres) > 4:
    # delete zeros for temperature plot
    radial_profiles_4_temp = radial_profiles[4]["temperature"][radial_profiles[4].used]
    radial_profiles_4_rad = radial_profiles[4].x[radial_profiles[4].used]

    # delete zeros for h2 fraction plot
    radial_profiles_4_h2 = radial_profiles[4]["H_p1_fraction"][radial_profiles[4].used]
    radial_profiles_4_h1 = radial_profiles[4]["H_p0_fraction"][radial_profiles[4].used] # for HI fraction

if len(spheres) > 5:
    # delete zeros for temperature plot
    radial_profiles_5_temp = radial_profiles[5]["temperature"][radial_profiles[5].used]
    radial_profiles_5_rad = radial_profiles[5].x[radial_profiles[5].used]

    # delete zeros for h2 fraction plot
    radial_profiles_5_h2 = radial_profiles[5]["H_p1_fraction"][radial_profiles[5].used]
    radial_profiles_5_h1 = radial_profiles[5]["H_p0_fraction"][radial_profiles[5].used] # for HI fraction


""" Set up fig and 1x2 subplots """

fig = plt.figure()
fig, axs = plt.subplots(1, 2)
plt.rcParams["font.family"] = "serif"
plt.rcParams['lines.linewidth'] = 2

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

font2 = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }

# 1) H2 molecule fraction vs. Radius
axs[0].loglog(radial_profiles_0_rad, radial_profiles_0_h2, color='b', linestyle='solid', label=label1)
axs[0].loglog(radial_profiles_1_rad, radial_profiles_1_h2, color='r', linestyle='solid', label=label2)
# axs[0].loglog(radial_profiles_2_rad, radial_profiles_2_h2, color = 'r', linestyle = 'dashed', label = '2Myr post-SF')
#
#
# axs[0].loglog(radial_profiles_3_rad, radial_profiles_3_h2, color = 'r', linestyle = 'dotted', label = '3Myr post-SF')
# axs[0].loglog(radial_profiles_4_rad, radial_profiles_4_h2, color = 'g', linestyle = 'solid', label = '3.8Myr lifetime of 40Msun')
# axs[0].loglog(radial_profiles_5_rad, radial_profiles_5_h2, color = 'g', linestyle = 'dashed', label = '0.1Myr post-BH')

#axs[0].legend(loc = "lower right", fontsize=12,  ncol=1)
axs[0].set_xlabel(r"r (pc)", fontdict=font)
axs[0].set_ylabel(r"Gas fraction in HII", fontdict=font)

axs[0].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[0].tick_params(axis="x", which='major', labelsize = 14, width=2, length=7, direction="in")
axs[0].tick_params(axis="y", which='major', labelsize = 14)


# 2) Temperature vs. Radius
axs[1].loglog(radial_profiles_0_rad, radial_profiles_0_temp, color='b', linestyle='solid', label=label1)
axs[1].loglog(radial_profiles_1_rad, radial_profiles_1_temp, color='r', linestyle='dashed', label=label2)
# axs[1].loglog(radial_profiles_2_rad, radial_profiles_2_temp, color = 'r', linestyle = 'dotted', label = '2Myr post-SF')
#
# axs[1].loglog(radial_profiles_3_rad, radial_profiles_3_temp, color = 'r', linestyle = 'solid', label = '3Myr post-SF')
# axs[1].loglog(radial_profiles_4_rad, radial_profiles_4_temp, color = 'g', linestyle = 'dashed', label = '3.8Myr lifetime of 40Msun')
# axs[1].loglog(radial_profiles_5_rad, radial_profiles_5_temp, color = 'g', linestyle = 'dotted', label = '0.1Myr post-BH')

axs[1].set_xlabel(r"r (pc)", fontdict=font)
axs[1].set_ylabel(r"T (K)", fontdict=font)

axs[1].legend(loc="upper left", fontsize=14,  ncol=1)

axs[1].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[1].tick_params(axis="x", which='major', labelsize = 14,  width=2, length=7, direction="in")
axs[1].yaxis.tick_right()
axs[1].yaxis.set_label_position("right")
axs[1].tick_params(axis="y", which='major', labelsize = 14)


# save 2 plots in in 1 figure
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(14.129921, 6.8661417)
plot_name = 'star-profile-plot-' + str(i) + '.pdf'
fig.savefig(plot_name, dpi=100)

