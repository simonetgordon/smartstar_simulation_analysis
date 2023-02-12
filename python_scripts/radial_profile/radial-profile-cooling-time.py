"""
Radial profile of cooling time, mass weighted
"""
import yt
import ytree
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

# macros
x = "seed1-bh-only-cooling-time-"
y = sys.argv[-1]  # plot  name
star_pos0 = [0.49030654, 0.50721054, 0.49178297] # found with smartstar-find.py

# load data
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/BF-0.0078125"

# all ds
DD = []
ds1 = yt.load(os.path.join(root_dir, sys.argv[1])) # t = 124.760, before particle formation
label1 = sys.argv[1]
DD.append(ds1)
ds2 = yt.load(os.path.join(root_dir, sys.argv[2])) # later time
label2 = sys.argv[2]
DD.append(ds2)
if str(sys.argv[3]).startswith('DD'):
    ds3 = yt.load(os.path.join(root_dir, sys.argv[3])) # latest time
    label3 = sys.argv[3]
    DD.append(ds3)

# define dark_matter_mass field + add to all ds
def _t_cool_abs(field, data):
    return np.abs(data[("gas", "cooling_time")].in_units("s"))


for ds in DD:
    ds.add_field(("gas", "cooling_time_abs"), function=_t_cool_abs, units="s", sampling_type="cell") # need to define with units here.



""" Make Sphere of Most Massive Halo """

# Load seed1 merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds1.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation with radius = 3 * virial radius
spheres = []
for i, ds in enumerate(DD):
    sphere = ds.sphere(star_pos0, 3*r)
    spheres.append(sphere)


""" Make Profile Plots weighted by cell mass """

all_fields = ["cell_mass", "matter_mass", "cooling_time", "cooling_time_abs"]
all_units = {"radius": "pc", "cell_mass": "Msun", "matter_mass": "Msun", "cooling_time_abs": "s", "cooling_time": "s"}

# Radial profile weighted by = 'cell_mass'
radial_profiles = []
for i, sp in enumerate(spheres):
    radial_profile = yt.create_profile(
        sp,
        "radius",
        all_fields,
        units=all_units,
        weight_field="cell_mass",
        accumulation=False,
        logs={"radius": True} # Set true for more bins
    )
    radial_profiles.append(radial_profile)


""" Plot """
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


plt.loglog(radial_profiles[0].x.value,
           radial_profiles[0][("gas","cooling_time_abs")],
           color='b', linestyle='solid', label=label1)
plt.loglog(radial_profiles[1].x.value,
           radial_profiles[1]["cooling_time_abs"],
           color='r', linestyle='solid', label=label2)
if str(sys.argv[3]).startswith('DD'):
    plt.loglog(radial_profiles[2].x[radial_profiles[2].used],
               radial_profiles[2]["cooling_time"][radial_profiles[2].used],
               color='r', linestyle='solid', label=label3)

plt.xlabel(r"r (pc)", fontdict=font)
plt.ylabel(r"t (s)", fontdict=font)

plt.legend(loc="lower left", fontsize=14,  ncol=1) # upper/lower

plt.tick_params(axis="x", which='minor', length = 4, direction="in")
plt.tick_params(axis="x", which='major', labelsize = 14,  width=2, length=7, direction="in")
# plt.yaxis.tick_right()
# plt.yaxis.set_label_position("right")
# plt.tick_params(axis="y", which='major', labelsize = 14)

plot_name = 'radial-profile-plot-' + str(x) + str(y) + '.pdf'
plt.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)