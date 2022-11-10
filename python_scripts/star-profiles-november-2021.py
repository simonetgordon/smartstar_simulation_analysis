# Bfore supernova plots in one 2x2 figure
# Star position currently set on position at time of formation

import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
from matplotlib.lines import Line2D # for custom legend

# Position of star at formation = [0.49053251 0.4946698  0.50963437]

# Load snapshots and add pop3/metallicity fields+filters

# star forms at t=124.55

# For Rays10
ds_10_1 = yt.load("DD0122/DD0122") # t = 124.55, Time of star formation.
ds_10_2 = yt.load("DD0222/DD0222") # t = 125.55, 1 Myr after star particle inserted
ds_10_3 = yt.load("DD0322/DD0322") # t = 126.55, 2 Myr after star particle inserted

# For Rays7
ds_7_1 = yt.load("DD0422/DD0422") # 127.55. 3 Myr after star particle inserted
ds_7_2 = yt.load("DD0508/DD0508") # 128.41,  3.87 Myr, last dd in stellar lifetime
ds_7_3 = yt.load("DD0520/DD0520") # 128.53 0.1Myr after BH is formed.

# For Rays5.1
ds_5_1 = yt.load("DD0122/DD0122")
ds_5_2 = yt.load("DD0222/DD0222")
ds_5_3 = yt.load("DD0322/DD0322")

# Make array of all data
DD = [ds_10_1, ds_10_2, ds_10_3, ds_7_1, ds_7_2, ds_7_3, ds_5_1, ds_5_2, ds_5_3]

# Define dark_matter_mass field

def _dm_mass(field, data):
    return data[("gas", "dark_matter_density")].in_units("g/cm**3") * data[("index", "cell_volume")].in_units("cm**3")

for ds in DD:
    # Add  dm_mass field to ds
    ds.add_field(("gas",'dm_mass'),function=_dm_mass, units = "g", sampling_type= 'cell') # need to define with units here.

    # Add p2p fields and filters
    #add_p2p_fields(ds)
    #add_p2p_particle_filters(ds)

## Make Sphere of Most Massive Halo ##

# Load merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('/disk01/sgordon/pop3/gas+dm-L3/rockstar_halos/out_0.list')

# Load my_tree and find radius
a1 = ytree.load('/disk01/sgordon/pop3/gas+dm-L3/tree_810/tree_810.h5')

r_halo = a1[0]["virial_radius"].to('pc')
r = ds_10_1.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 5 * virial radius

star_pos0 = [0.49053118, 0.49467636, 0.50964283]
sp10_a = ds_10_1.sphere(star_pos0, 5*r)
sp10_b = ds_10_2.sphere(star_pos0, 5*r)
sp10_c = ds_10_3.sphere(star_pos0, 5*r)

sp7_a = ds_7_1.sphere(star_pos0, 5*r)
sp7_b = ds_7_2.sphere(star_pos0, 5*r)
sp7_c = ds_7_3.sphere(star_pos0, 5*r)

sp5_a = ds_5_1.sphere(star_pos0, 5*r)
sp5_b = ds_5_2.sphere(star_pos0, 5*r)
sp5_c = ds_5_3.sphere(star_pos0, 5*r)
 # need to do individual spheres


## Make Profile Plots ##

    
# Print out Pop III star mass and make sphere centred on current star

#star_pos10_1 = sp10_a['Smartstar', 'particle_position']
#star_pos10_2 = sp10_b['pop3', 'particle_position']
#star_pos10_3 = sp10_c['pop3', 'particle_position']#

#star_pos7_1 = sp7_a['pop3', 'particle_position']
#star_pos7_2 = sp7_b['pop3', 'particle_position']
#star_pos7_3 = sp7_c['pop3', 'particle_position']

#star_pos5_1 = sp5_a['pop3', 'particle_position']
#star_pos5_2 = sp5_b['pop3', 'particle_position']
#star_pos5_3 = sp5_c['pop3', 'particle_position']

# sphere around current star position. Changed to star_pos0 to centre on star at time of formation (for comparison purposes)
sp10_1 = ds_10_1.sphere(star_pos0, 3*r)
sp10_2 = ds_10_2.sphere(star_pos0, 3*r)
sp10_3 = ds_10_3.sphere(star_pos0, 3*r)

sp7_1 = ds_7_1.sphere(star_pos0, 3*r)
sp7_2 = ds_7_2.sphere(star_pos0, 3*r)
sp7_3 = ds_7_3.sphere(star_pos0, 3*r)

sp5_1 = ds_5_1.sphere(star_pos0, 3*r)
sp5_2 = ds_5_2.sphere(star_pos0, 3*r)
sp5_3 = ds_5_3.sphere(star_pos0, 3*r)

all_spheres = [sp10_1 , sp10_2, sp10_2, sp7_1, sp7_2, sp7_3, sp5_1, sp5_2, sp5_3]


# Bin up the data from the sphere into radial profiles.

# Radial profile weighted by default = 'cell_mass'. (WAS RP)
rp10_mass_1 = yt.create_profile(
    sp10_1,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p1_fraction", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "g/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp10_mass_1_temp = rp10_mass_1["temperature"][rp10_mass_1.used]
rp10_mass_1_rad = rp10_mass_1.x[rp10_mass_1.used]

rp10_mass_1_h2 = rp10_mass_1["H_p1_fraction"][rp10_mass_1.used] # changed from H2 TO hii
rp10_mass_1_h1 = rp10_mass_1["H_p0_fraction"][rp10_mass_1.used] # for HI fraction


rp10_mass_2 = yt.create_profile(
    sp10_2,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "g/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp10_mass_2_temp = rp10_mass_2["temperature"][rp10_mass_2.used]
rp10_mass_2_rad = rp10_mass_2.x[rp10_mass_2.used]


# delete zeros for h2 fraction plot
rp10_mass_2_h2 = rp10_mass_2["H_p1_fraction"][rp10_mass_2.used]
rp10_mass_2_h1 = rp10_mass_2["H_p0_fraction"][rp10_mass_2.used] # for HI fraction

rp10_mass_3 = yt.create_profile(
    sp10_3,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "g/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp10_mass_3_temp= rp10_mass_3["temperature"][rp10_mass_3.used]
rp10_mass_3_rad = rp10_mass_3.x[rp10_mass_3.used]

# delete zeros for h2 fraction plot
rp10_mass_3_h2 = rp10_mass_3["H_p1_fraction"][rp10_mass_3.used]
rp10_mass_3_h1 = rp10_mass_3["H_p0_fraction"][rp10_mass_3.used] # for HI fraction


# Radial profile weighted by default = 'cell_mass'.
rp7_mass_1 = yt.create_profile(
    sp7_1,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp7_mass_1_temp = rp7_mass_1["temperature"][rp7_mass_1.used]
rp7_mass_1_rad = rp7_mass_1.x[rp7_mass_1.used]

# delete zeros for h2 fraction plot
rp7_mass_1_h2 = rp7_mass_1["H_p1_fraction"][rp7_mass_1.used]
rp7_mass_1_h1 = rp7_mass_1["H_p0_fraction"][rp7_mass_1.used] # for HI fraction

rp7_mass_2 = yt.create_profile(
    sp7_2,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp7_mass_2_temp = rp7_mass_2["temperature"][rp7_mass_2.used]
rp7_mass_2_rad = rp7_mass_2.x[rp7_mass_2.used]

# delete zeros for h2 fraction plot
rp7_mass_2_h2 = rp7_mass_2["H_p1_fraction"][rp7_mass_2.used]
rp7_mass_2_h1 = rp7_mass_2["H_p0_fraction"][rp7_mass_2.used] # for HI fraction

rp7_mass_3 = yt.create_profile(
    sp7_3,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp7_mass_3_temp = rp7_mass_3["temperature"][rp7_mass_3.used]
rp7_mass_3_rad = rp7_mass_3.x[rp7_mass_3.used]

# delete zeros for h2 fraction plot
rp7_mass_3_h2 = rp7_mass_3["H_p1_fraction"][rp7_mass_3.used]
rp7_mass_3_h1 = rp7_mass_3["H_p0_fraction"][rp7_mass_3.used] # for HI fraction

# Radial profile weighted by default = 'cell_mass'.
rp5_mass_1 = yt.create_profile(
    sp5_1,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)

# delete zeros for temperature plot
rp5_mass_1_temp = rp5_mass_1["temperature"][rp5_mass_1.used]
rp5_mass_1_rad = rp5_mass_1.x[rp5_mass_1.used]

# delete zeros for h2 fraction plot
rp5_mass_1_h2 = rp5_mass_1["H_p1_fraction"][rp5_mass_1.used]
rp5_mass_1_h1 = rp5_mass_1["H_p0_fraction"][rp5_mass_1.used]


rp5_mass_2 = yt.create_profile(
    sp5_2,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)


rp5_mass_2_temp = rp5_mass_2["temperature"][rp5_mass_2.used]
rp5_mass_2_rad = rp5_mass_2.x[rp5_mass_2.used]

# delete zeros for h2 fraction plot
rp5_mass_2_h2 = rp5_mass_2["H_p1_fraction"][rp5_mass_2.used]
rp5_mass_2_h1 = rp5_mass_2["H_p0_fraction"][rp5_mass_2.used]

rp5_mass_3 = yt.create_profile(
    sp5_3,
    "radius",
    ["density", "temperature", "H_p1_fraction", "radial_velocity", "cell_mass", "matter_mass", "cell_volume", "H_p0_fraction"], # 'dm_mass' raises error if put here.
    units={"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "density": "Msun/cm**3"},
    logs={"radius": True} # Set true for more bins
)

rp5_mass_3_temp = rp5_mass_3["temperature"][rp5_mass_3.used]
rp5_mass_3_rad = rp5_mass_3.x[rp5_mass_3.used]

 #delete zeros for h2 fraction plot
rp5_mass_3_h2 = rp5_mass_3["H_p1_fraction"][rp5_mass_3.used]
rp5_mass_3_h1 = rp5_mass_3["H_p0_fraction"][rp5_mass_3.used]



# Radial profile weighted by 'cell_volume'. For densities. (WAS RP2)

# Rays10
rp10_vol_1 = yt.create_profile(
    sp10_1,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density", "number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp10_vol_2 = yt.create_profile(
    sp10_2,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp10_vol_3 = yt.create_profile(
    sp10_3,
    "radius",
    [ "cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)


# Rays7

rp7_vol_1 = yt.create_profile(
    sp7_1,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp7_vol_2 = yt.create_profile(
    sp7_2,
    "radius",
    ["radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp7_vol_3 = yt.create_profile(
    sp7_3,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)


# Ray5

rp5_vol_1 = yt.create_profile(
    sp5_1,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp5_vol_2 = yt.create_profile(
    sp5_2,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)

rp5_vol_3 = yt.create_profile(
    sp5_3,
    "radius",
    ["cell_volume", "radius", "matter_density", "density", "Dark_Matter_Density","number_density"], # 'dm_mass' also raises error here.
    units={"radius": "pc", "Dark_Matter_Density": "g/cm**3", "matter_density": "g/cm**3", "density": "g/cm**3"},
    weight_field='cell_volume',
    accumulation=True,
    logs={"radius": True} # Set true for more bins
)


# Set up fig and 1x2 subplots

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

axs[0].loglog(rp10_mass_1_rad, rp10_mass_1_h2, color = 'b', linestyle = 'solid', label = 'Star forms = 124.55 Myr')
axs[0].loglog(rp10_mass_2_rad, rp10_mass_2_h2, color = 'r', linestyle = 'solid', label = '1Myr post-SF')
axs[0].loglog(rp10_mass_3_rad, rp10_mass_3_h2, color = 'r', linestyle = 'dashed', label = '2Myr post-SF')


axs[0].loglog(rp7_mass_1_rad, rp7_mass_1_h2, color = 'r', linestyle = 'dotted', label = '3Myr post-SF')
axs[0].loglog(rp7_mass_2_rad, rp7_mass_2_h2, color = 'g', linestyle = 'solid', label = '3.8Myr lifetime of 40Msun')
axs[0].loglog(rp7_mass_3_rad, rp7_mass_3_h2, color = 'g', linestyle = 'dashed', label = '0.1Myr post-BH')

#axs[0].legend(loc = "lower right", fontsize=12,  ncol=1)

axs[0].set_xlabel(r"r (pc)",  fontdict=font)
axs[0].set_ylabel(r"Gas fraction in HII", fontdict=font)

axs[0].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[0].tick_params(axis="x", which='major', labelsize = 14, width=2, length=7, direction="in")
axs[0].tick_params(axis="y", which='major', labelsize = 14)


# 2) Temperature vs. Radius


axs[1].loglog(rp10_mass_1_rad, rp10_mass_1_temp, color = 'b', linestyle = 'solid', label = 'Star forms = 124.55 Myr')
axs[1].loglog(rp10_mass_2_rad, rp10_mass_2_temp, color = 'r', linestyle = 'dashed', label = '1Myr post-SF')
axs[1].loglog(rp10_mass_3_rad, rp10_mass_3_temp, color = 'r', linestyle = 'dotted', label = '2Myr post-SF')

axs[1].loglog(rp7_mass_1_rad, rp7_mass_1_temp, color = 'r', linestyle = 'solid', label = '3Myr post-SF')
axs[1].loglog(rp7_mass_2_rad, rp7_mass_2_temp, color = 'g', linestyle = 'dashed', label = '3.8Myr lifetime of 40Msun')
axs[1].loglog(rp7_mass_3_rad, rp7_mass_3_temp, color = 'g', linestyle = 'dotted', label = '0.1Myr post-BH')

axs[1].set_xlabel(r"r (pc)", fontdict=font)
axs[1].set_ylabel(r"T (K)", fontdict=font)

axs[1].legend(loc = "upper left", fontsize=14,  ncol=1)

axs[1].tick_params(axis="x", which='minor', length = 4, direction="in")
axs[1].tick_params(axis="x", which='major', labelsize = 14,  width=2, length=7, direction="in")
axs[1].yaxis.tick_right()
axs[1].yaxis.set_label_position("right")
axs[1].tick_params(axis="y", which='major', labelsize = 14)




# save 2 plots in in 1 figure

fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(14.129921, 6.8661417)
fig.savefig('star-profile-40msun-seed1-bondi1-tonale.pdf', dpi=100)

