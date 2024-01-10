"""
For producing radial profiles of advective vs radiative cooling in nuclear disc
python -i plot_cooling_rates.py DD0130/DD0130
"""

import yt
import sys
import os
import numpy as np
import pandas as pd
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from derived_fields import add_fields_ds
from helper_functions import _make_disk_L, adaptive_moving_average, extract_colors
from matplotlib import rc
import matplotlib.transforms as mtransforms


# my data 
root_dir = ["/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann",
            "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann", 
            "/Backup00/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann",
            "/Backup00/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2",
            "/Backup00/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2"]
sim = [
    "1B.RSb01-2", 
    "1B.RSb16",
    "1B.m16-4dx",
    "2B.RSb08/2B.RSb08-2",
    "2B.m08-4dx/2B.m16-4dx-2",
    #"2B.m08-4dx"
    ] 
dds = [
    "DD0138/DD0138", # 1 Myr
    "DD0166/DD0166", # 1 Myr
    "DD0202/DD0202", # 1.000 Myr
    "DD0279/DD0279", # 1 Myr
    #"DD0299/DD0299", # 1 Myr in 2B.m08-4dx
    "DD0304/DD0304"] # 1 Myr in 2B.m16-4dx

labels = []
DS = []

# Beckmann data
data_beck = pd.read_csv("cooling_line_beckmann.csv")
radius_beck = data_beck['radius'].tolist()
q_ratio_beck = data_beck['qratio'].tolist()

# set font properties
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18
plt.rcParams["lines.linewidth"] = 2
#font_properties = {'family': 'serif', 'serif': ['Times'], 'weight': 'light'}
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
plt.rcParams["mathtext.default"] = "regular" # for the same font used in regular text.
rc('text', usetex=True)

# colors
c_s1 = extract_colors('viridis', 3, portion="middle")
c_s2 = extract_colors('magma', 2, portion="middle", start=0.4, end=0.6)
c = np.concatenate((c_s1, c_s2))

# set up figure
fig, axs = plt.subplots(1, sharex=True)

# set x and y labels
plt.xlabel(r"Radius (pc)", fontsize=20)
plt.ylabel(r"$\rm Q_{adv} / Q_{rad}$", fontsize=20)
#plt.ylabel(r"$\rm Q_{rad} (\rm erg \, s^{-1} \, cm^{-3})$", fontsize=20)
# set minorticks
min_n_index = -3
max_n_index = 1
a1 = np.logspace(min_n_index, max_n_index, np.abs(min_n_index - max_n_index) + 1)
a2 = np.arange(1, 10, 1)
minorticks = np.outer(a1, a2).flatten()
plt.yticks(minorticks, minor=True)

# plot Beckmann cooling line
plt.loglog(radius_beck, q_ratio_beck, color="grey", label="B18",linewidth=4)

labels = ['1B.b01', '1B.b16', '1B.m16', '2B.b08', '2B.m16']
# Initialize lists to store aggregated data
all_x = []
all_y = []
for j, ds in enumerate(dds):
    # Load the data from the local directory
    ds = yt.load(os.path.join(root_dir[j], sim[j], ds))

    # add all cooling rate derived fields to ds
    add_fields_ds(ds)
    
    # make profile from disk object
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    disk, L = _make_disk_L(ds, ss_pos, 10*yt.units.pc, 0.2*yt.units.pc)
    profile = yt.create_profile(
        data_source=disk,
        bin_fields=[("index", "radius")],
        fields=[("enzo", "radiative_cooling_rate"), ("enzo", "advective_cooling_rate")],
        n_bins=512,
        units=dict(radius="pc"),
        logs=dict(radius=True),
        weight_field=("gas", "mass"),
    )

    # take absolute value of cooling rate
    cooling_ratio = np.abs(profile[("enzo", "radiative_cooling_rate")][profile.used]/profile[("enzo", "advective_cooling_rate")][profile.used])
    n = 10
    x_avg = adaptive_moving_average(profile.x[profile.used], window_size=n)
    y_avg = adaptive_moving_average(cooling_ratio, window_size=n)
    x_avg = x_avg.d
    y_avg = y_avg.d

    # Append the data to the aggregated lists
    all_x.append(x_avg)
    all_y.append(y_avg)

    #plt.loglog(profile.x[profile.used], cooling_ratio)
    plt.loglog(x_avg, y_avg, color=c[j], label=labels[j])
    plt.axhline(y=1, color='grey', linestyle='dashed', linewidth=2, alpha=1)
    plt.xlim(4e-4, 10.2)
    #plt.ylim(2e-4, 1.2e4)
    plt.legend(ncol=1, fontsize=16, frameon=False, handlelength=2, handletextpad=0.8, loc='lower right')

    # annotate with simulation data
    # 1) BH Age in text box
    mass_label = 'BH Mass = {:.0f}'.format(ss_mass)
    print(mass_label)
    #props = dict(boxstyle='square', facecolor='white', ec='black', lw=1, alpha=0.5)
    #plt.text(0.2, 5e-4, label, verticalalignment='top', bbox=props)

    # 2) Denote disc region
    rdiscs = [0.15, 
              #0.04
              ]
    disc_r = rdiscs[0] # pc -from projection plot inspection
    accretion_r = 1.23e-2
    cr = ["lemonchiffon", "lemonchiffon"]
    plt.axvline(accretion_r, color="black", alpha=1)
    trans = mtransforms.blended_transform_factory(axs.transData, axs.transAxes)
    if j == 1:
        plt.fill_between(adaptive_moving_average(profile.x[profile.used], window_size=n), 0, adaptive_moving_average(cooling_ratio, window_size=n).max(), 
                         where=adaptive_moving_average(profile.x[profile.used]) <= disc_r, facecolor=cr[0], transform=trans, alpha=0.5)


# After the loop, compute the trendline for the aggregated data
# Ensure all_x and all_y are flat arrays with numeric types.
# This assumes you want to use floating-point numbers for the fit.
all_x = np.hstack(all_x).astype(float)
all_y = np.hstack(all_y).astype(float)

# Fit a linear trendline to the aggregated data
z = np.polyfit(all_x, all_y, 1)  # Adjust the degree (1 for linear) as needed
p = np.poly1d(z)

# Create a range of x values for plotting the trendline
trend_x = np.linspace(min(all_x), max(all_x), len(all_x))
trend_y = p(trend_x)

# Plot the trendline
#plt.plot(trend_x, trend_y, linestyle="--", color='black', label='Trendline')

# Save the figure
plt.title("Advection/Radiative Cooling Rate Radial Profile", fontsize=20)
plot_name = f'radial-profile-radiative-cooling-ratio-{str(sim[0])}_{str(sim[1])}_window={n}.pdf'
fig.set_size_inches(8, 5.8)
plt.savefig('plots/' + plot_name, dpi=100)
print("created ", 'plots/' + plot_name)
