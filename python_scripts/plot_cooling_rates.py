"""
For producing radial profiles of advective vs radiative cooling in nuclear disc
python plot_cooling_rates.py DD0130/DD0130
"""

import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds
from find_disc_attributes import make_disk

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# add all cooling rate derived fields to ds
add_fields_ds(ds)
ad = ds.all_data()

# define cooling rates
specific_energy = ad['enzo', 'GasEnergy'].to('ergs/g') # g cm^2 s^-2 cm^-3 = g cm^-1 s^-2 = ergs g^-1
print("q_energy: ", specific_energy)

r = ad["index", "radius"].to('pc')
q_sim = ad["enzo", "simulation_cooling_rate"] # [erg g^-1 s^-1] = (g cm^2 s^-2) g^-1 s^-1 = cm^2 s^-3
q_rad = ad["enzo", "radiative_cooling_rate"]
q_adv = ad["enzo", "advective_cooling_rate"]

# radiative cooling rate [erg s-1 cm-3]
print("q_rad: ", q_rad)

# advective cooling rate [erg s-1 cm-3]
print("q_adv: ", q_adv)

# set up figure
fig = plt.figure()
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"
linewidth = 1
plt.rcParams['lines.linewidth'] = linewidth

fontsize = 10
font = {'family': 'serif',
        'weight': 'normal',
        'size': fontsize,
        }
plt.xlabel(r"$Radius \, (pc)$", fontdict=font)
plt.ylabel(r"$Q_{adv} / Q_{rad}$", fontdict=font)

# make profile from disk object
ss_pos, ss_mass, ss_age = ss_properties(ds)
disk = make_disk(ss_pos, 10, 0.2)
profile = yt.create_profile(
    data_source=disk,
    bin_fields=[("index", "radius")],
    fields=[("enzo", "radiative_cooling_rate"), ("enzo", "advective_cooling_rate")],
    n_bins=260,
    units=dict(radius="pc"),
    logs=dict(radius=True),
    weight_field=("gas", "mass"),
)

# plot and save figure
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


plt.loglog(profile.x[profile.used], profile[("enzo", "radiative_cooling_rate")][profile.used] /
           profile[("enzo", "advective_cooling_rate")][profile.used])
plt.axhline(y=1, color='black', linestyle='dashed', lw=linewidth, alpha=1)

# naming plot
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84

plot_name = 'radial-profile-plot-cooling-rate-' + str(root_dir[index:]) + '_' + str(input)[7:] + '.pdf'
plt.savefig('plots/' + plot_name, dpi=100)
print("created ", 'plots/' + plot_name)
