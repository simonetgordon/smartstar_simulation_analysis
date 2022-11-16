""""
Radial profiles in one 2x2 figure
Particle position set on position at time of formation and must be found prior
"""
import yt
import ytree
import os
import sys
import matplotlib.pyplot as plt

# macros
x = yt.load(sys.argv[1]) # plot number
mass_density = True # to create volume-weight radial profiles
star_pos0 = [0.49048811, 0.49467262, 0.50964459] # found with smartstar-find.py

# load data
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun-thermal-only/BF-0.0078125"

# all ds
DD = []
ds1 = yt.load(os.path.join(root_dir, "DD0125/DD0125")) # t = 124.760, before particle formation
label1 = "DD0125"
DD.append(ds1)
ds2 = yt.load(os.path.join(root_dir, "DD0126/DD0126")) # t = 124.7618, 0.0004323 Myr (400 yrs) after particle formation
label2 = "DD0126"
DD.append(ds2)


# define dark_matter_mass field + add to all ds
def _dm_mass(field, data):
    return data[("gas", "dark_matter_density")].in_units("g/cm**3") * data[("index", "cell_volume")].in_units("cm**3")


for ds in DD:
    ds.add_field(("gas", "dm_mass"), function=_dm_mass, units="g", sampling_type="cell") # need to define with units here.


""" Make Sphere of Most Massive Halo """

# Load seed1 merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds1.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 3 * virial radius
spheres = []
for i, ds in enumerate(DD):
    sphere = ds.sphere(star_pos0, 3*r)
    spheres.append(sphere)


""" Make Profile Plots weighted by cell mass """

all_fields = ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass",
              "cell_volume", "H_p1_fraction", "H_p0_fraction", "number_density"] # 'dm_mass' raises error if put here.
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

if mass_density:
    # Radial profile weighted by default = 'cell_volume'
    radial_profiles_vol = []
    for i, sp in enumerate(spheres):
        radial_profile = yt.create_profile(
            sp,
            "radius",
            all_fields,
            units=all_units,
            weight_field="cell_volume",
            accumulation=True,
            logs={"radius": True} # Set true for more bins
        )
        radial_profiles_vol.append(radial_profile)


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
fig, axs = plt.subplots(2, 2, sharex=True)
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
axs[0, 0].loglog(radial_profiles_0_rad, radial_profiles_0_h2, color='b', linestyle='solid', label=label1)
axs[0, 0].loglog(radial_profiles_1_rad, radial_profiles_1_h2, color='r', linestyle='solid', label=label2)

axs[0, 0].set_xlabel(r"r (pc)", fontdict=font)
axs[0, 0].set_ylabel(r"Gas fraction in HII", fontdict=font)

axs[0, 0].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[0, 0].tick_params(axis="x", which='major', labelsize = 14, width=2, length=7, direction="in")
axs[0, 0].tick_params(axis="y", which='major', labelsize = 14)


# 2) Temperature vs. Radius
axs[0, 1].loglog(radial_profiles_0_rad, radial_profiles_0_temp, color='b', linestyle='solid', label=label1)
axs[0, 1].loglog(radial_profiles_1_rad, radial_profiles_1_temp, color='r', linestyle='solid', label=label2)

axs[0, 1].set_xlabel(r"r (pc)", fontdict=font)
axs[0, 1].set_ylabel(r"T (K)", fontdict=font)

axs[0, 1].legend(loc="lower left", fontsize=14,  ncol=1) # upper/lower

axs[0, 1].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[0, 1].tick_params(axis="x", which='major', labelsize = 14,  width=2, length=7, direction="in")
axs[0, 1].yaxis.tick_right()
axs[0, 1].yaxis.set_label_position("right")
axs[0, 1].tick_params(axis="y", which='major', labelsize = 14)


# 3) Mass Density
if mass_density:
    axs[1, 0].loglog(radial_profiles_vol[0].x.value, radial_profiles_vol[0]["number_density"].value, color='b', linestyle='solid', label=label1)
    axs[1, 0].loglog(radial_profiles_vol[1].x.value, radial_profiles_vol[1]["number_density"].value, color='r', linestyle='solid', label=label2)

    axs[1, 0].set_ylabel(r"$\mathrm{n\ (cm^{-3})}$", fontsize=20)

    axs[1, 0].set_xticklabels([])
    axs[1, 0].tick_params(axis="x", which='minor', length = 4, direction="in")
    axs[1, 0].tick_params(axis="x", which='major', width=2, length=7, direction="in")
    axs[1, 0].tick_params(axis="y", which='major', labelsize=14)


# 4) Radial Velocity
# first we need to find the bulk velocity and subtract it from the vels
for sp in spheres:
    bulk_vel = sp.quantities.bulk_velocity()
    sp.set_field_parameter("bulk_velocity", bulk_vel)

axs[1, 1].semilogx(radial_profiles[0].x.value, radial_profiles[0]["radial_velocity"].in_units("km/s").value, color='b', linestyle='solid', label=label1)
axs[1, 1].semilogx(radial_profiles[1].x.value, radial_profiles[1]["radial_velocity"].in_units("km/s").value, color='r', linestyle='solid', label=label2)

axs[1, 1].set_xlabel(r"$\mathrm{r\ (pc)}$", fontsize=20)
axs[1, 1].set_ylabel(r"$\mathrm{Radial\ velocity\ (km/s)} $", fontsize=20)

axs[1, 1].tick_params(axis="x", which='minor', length=4)
axs[1, 1].tick_params(axis="x", which='major', labelsize=18, width=2, length=7)
axs[1, 1].tick_params(axis="y", which='major', labelsize=18)
axs[1, 1].yaxis.tick_right()
axs[1, 1].yaxis.set_label_position("right")

# save 2 plots in 1 figure
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(14.129921, 6.8661417)
plot_name = 'radial-profile-plot-' + str(x) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)

