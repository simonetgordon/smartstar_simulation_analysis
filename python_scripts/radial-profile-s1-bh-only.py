""""
Radial profiles in one 2x2 figure
Particle position set on position at time of formation and must be found prior
"""
import yt
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

# macros
x = "seed1-bh-only-" # plot number
y = sys.argv[-1] # naming plot
mass_density = True # to create volume-weight radial profiles
star_pos0 = [0.49048772, 0.49467235, 0.50964457] # found with smartstar-find.py

# load data
root_dir = "/work/sc070/sc070/stg/simulations/seed1-bh-only/270msun/replicating-beckmann/1f.RSb01"

# all ds
DD = []
ds1 = yt.load(os.path.join(root_dir, sys.argv[1])) # before particle formation
label1 = "before formation"
DD.append(ds1)
ds2 = yt.load(os.path.join(root_dir, sys.argv[2])) # later time
label2 = "10 kyr"
DD.append(ds2)
if str(sys.argv[3]).startswith('DD'):
    ds3 = yt.load(os.path.join(root_dir, sys.argv[3])) # latest time
    label3 = "30 kyr"
    DD.append(ds3)
if str(sys.argv[4]).startswith('DD'):
    ds4 = yt.load(os.path.join(root_dir, sys.argv[4])) # latest time
    label4 = "100 kyr"
    DD.append(ds4)

# define dark_matter_mass field + add to all ds
def _dm_mass(field, data):
    return data[("gas", "dark_matter_density")].in_units("g/cm**3") * data[("index", "cell_volume")].in_units("cm**3")


for ds in DD:
    ds.add_field(("gas", "dm_mass"), function=_dm_mass, units="g", sampling_type="cell")


""" Make Sphere of Most Massive Halo """

# Make initial sphere centred on the star at the time of formation with radius = 3 * virial radius
spheres = []
r = ds1.quan(1000, "pc") # virial radius
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

fontsize = 8
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }

c = pl.cm.tab20b(np.linspace(0, 1, len(spheres)))

# 1) H2 molecule fraction vs. Radius
axs[0, 0].loglog(radial_profiles_0_rad, radial_profiles_0_h2, color=c[0], linestyle='solid', label=label1)
axs[0, 0].loglog(radial_profiles_1_rad, radial_profiles_1_h2, color=c[1], linestyle='solid', label=label2)
if str(sys.argv[3]).startswith('DD'):
    axs[0, 0].loglog(radial_profiles_2_rad, radial_profiles_2_h2, color=c[2], linestyle='solid', label=label3)
if str(sys.argv[4]).startswith('DD'):
    axs[0, 0].loglog(radial_profiles_3_rad, radial_profiles_3_h2, color=c[3], linestyle='solid', label=label4)

axs[0, 0].set_xlabel(r"r (pc)", fontdict=font)
axs[0, 0].set_ylabel(r"Gas fraction in HII", fontdict=font)

axs[0, 0].tick_params(axis="x", which='minor', length=4, direction="in")
axs[0, 0].tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=7, direction="in")
axs[0, 0].tick_params(axis="y", which='major', labelsize=fontsize)


# 2) Temperature vs. Radius
axs[0, 1].loglog(radial_profiles_0_rad, radial_profiles_0_temp, color=c[0], linestyle='solid', label=label1)
axs[0, 1].loglog(radial_profiles_1_rad, radial_profiles_1_temp, color=c[1], linestyle='solid', label=label2)
if str(sys.argv[3]).startswith('DD'):
    axs[0, 1].loglog(radial_profiles_2_rad, radial_profiles_2_temp, color=c[2], linestyle='solid', label=label3)
if str(sys.argv[4]).startswith('DD'):
    axs[0, 1].loglog(radial_profiles_3_rad, radial_profiles_3_temp, color=c[3], linestyle='solid', label=label4)

axs[0, 1].set_xlabel(r"r (pc)", fontdict=font)
axs[0, 1].set_ylabel(r"T (K)", fontdict=font)

axs[0, 1].tick_params(axis="x", which='minor', length=4, direction="in")
axs[0, 1].tick_params(axis="x", which='major', labelsize=fontsize,  width=2, length=7, direction="in")
axs[0, 1].yaxis.tick_right()
axs[0, 1].yaxis.set_label_position("right")
axs[0, 1].tick_params(axis="y", which='major', labelsize=fontsize)


# 3) Mass Density
if mass_density:
    axs[1, 0].loglog(radial_profiles_vol[0].x.value, radial_profiles_vol[0]["number_density"].value, color=c[0],
                     linestyle='solid', label=label1)
    axs[1, 0].loglog(radial_profiles_vol[1].x.value, radial_profiles_vol[1]["number_density"].value, color=c[1],
                     linestyle='solid', label=label2)
    if str(sys.argv[3]).startswith('DD'):
        axs[1, 0].loglog(radial_profiles_vol[2].x.value, radial_profiles_vol[2]["number_density"].value, color=c[2],
                         linestyle='solid', label=label3)
    if str(sys.argv[4]).startswith('DD'):
        axs[1, 0].loglog(radial_profiles_vol[3].x.value, radial_profiles_vol[3]["number_density"].value, color=c[3],
                         linestyle='solid', label=label4)
    
    axs[1, 0].set_ylabel(r"n (H $\mathrm{cm^{-3}})$", fontdict=font)
    axs[1, 0].set_xlabel("r (pc)", fontdict=font)

    axs[1, 0].legend(loc="upper right", fontsize=fontsize, ncol=1)  # upper/lower
    axs[1, 0].tick_params(axis="x", which='minor', length=4, direction="in")
    axs[1, 0].tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=7, direction="in")
    axs[1, 0].tick_params(axis="y", which='major', labelsize=fontsize)


# 4) Radial Velocity
# first we need to find the bulk velocity and subtract it from the vels
for sp in spheres:
    bulk_vel = sp.quantities.bulk_velocity()
    sp.set_field_parameter("bulk_velocity", bulk_vel)

axs[1, 1].semilogx(radial_profiles[0].x.value, radial_profiles[0]["radial_velocity"].in_units("km/s")[radial_profiles[0].used], color=c[0],
                   linestyle='solid', label=label1)
axs[1, 1].semilogx(radial_profiles[1].x.value, radial_profiles[1]["radial_velocity"].in_units("km/s").value, color=c[1],
                   linestyle='solid', label=label2)
if str(sys.argv[3]).startswith('DD'):
    axs[1, 1].semilogx(radial_profiles[2].x.value, radial_profiles[2]["radial_velocity"].in_units("km/s").value, color=c[2],
                       linestyle='solid', label=label3)
if str(sys.argv[4]).startswith('DD'):
    axs[1, 1].loglog(radial_profiles[3].x.value, radial_profiles[3]["radial_velocity"].in_units("km/s").value, color=c[3],
                     linestyle='solid', label=label4)


axs[1, 1].set_xlabel("r (pc)", fontdict=font)
axs[1, 1].set_ylabel("Radial velocity (km/s)", fontdict=font)

axs[1, 1].tick_params(axis="x", which='minor', length=4, direction="in")
axs[1, 1].tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=7, direction="in")
axs[1, 1].tick_params(axis="y", which='major', labelsize=fontsize)
axs[1, 1].yaxis.tick_right()
axs[1, 1].yaxis.set_label_position("right")

# save 2 plots in 1 figure
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(8.129921, 4.8661417)
plot_name = 'radial-profile-plot-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)

