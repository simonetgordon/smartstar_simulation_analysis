"""
For producing radial profiles of advective vs radiative cooling in nuclear disc
python -i plot_cooling_rates.py DD0130/DD0130
"""

import yt
import shutil
import sys
import os
import numpy as np
import pandas as pd
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from derived_fields import add_fields_ds
from find_disc_attributes import _make_disk_L
from matplotlib import rc
import matplotlib.transforms as mtransforms
from plot_multi_projections import tidy_data_labels
import matplotlib.cm as cm

# plot and save figure
def moving_average(a, n=5) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# function that does the interpolation
def interpolate_data(arr, N=100):
    # interpolate arr over `N` evenly spaced points
    min_val = np.min(arr)
    max_val = np.max(arr)

    t_orig = np.linspace(min_val, max_val, len(arr))
    t_interp = np.linspace(min_val, max_val, N)
    f = interp1d(x=t_orig, y=arr)
    interp_arr = f(t_interp)
    return interp_arr


def extract_colors(cmap_name, n, portion=None, start=None, end=None):
    cmap = cm.get_cmap(cmap_name)

    if start is not None and end is not None:
        values = np.linspace(start, end, n)
    elif portion == "beginning":
        values = np.linspace(0, 0.3, n)
    elif portion == "middle":
        values = np.linspace(0.3, 0.95, n)
    elif portion == "end":
        values = np.linspace(0.7, 1, n)
    elif portion is None:
        values = np.linspace(0, 1, n)
    else:
        raise ValueError("Invalid portion specified. Use 'beginning', 'middle', 'end', or None.")

    colors = cmap(values)
    return colors

# my data 
root_dir = ["/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann",
            "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann", 
            "/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann",
            "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2",
            "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2"]
sim = [
    "1B.RSb01-2", 
    "1B.RSb16",
    "1B.m16-4dx",
    "2B.RSb08/2B.RSb08-2",
    "2B.m08-4dx"]
dds = [
    "DD0138/DD0138", # 1 Myr
    "DD0166/DD0166", # 1 Myr
    "DD0202/DD0202", # 1.000 Myr
    "DD0252/DD0252", # 0.750 Myr
    "DD0291/DD0291"] # 0.934 Myr
labels = []
DS = []

# Beckmann data
data_beck = pd.read_csv("beckmann-data/cooling_line_beckmann.csv")
radius_beck = data_beck['radius'].tolist()
q_ratio_beck = data_beck['qratio'].tolist()

## First check there is a local /data area                                                                                                                                                                 
# if os.path.isdir("/data"):
#     #sys.exit("Error: no /data")

#     ## Second, check we have a directory. If not, then create one.                                                                                                                                             
#     UserName=os.getlogin()
#     LocalDir=os.path.join("/data",UserName)
#     if not os.path.isdir(LocalDir):
#         print("Creating Directory "+LocalDir)
#         os.mkdir(LocalDir)

#     ## Third, check if the data is already there, and if not, copy it over.                                                                                                                                    
#     DataDumpFull = sys.argv[1]
#     DataDump = DataDumpFull.split('/')
#     LocalData = os.path.join(LocalDir,DataDump[0])
#     if not os.path.isdir(LocalData):
#         print("Copying data to "+LocalData)
#         shutil.copytree(os.path.join(root_dir,DataDump[0]),LocalData)
#         print("Done copying data")
#     else:
#         print("Found a local copy in "+LocalData)
    
#     root_dir = LocalDir

# set font properties
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 18
plt.rcParams["lines.linewidth"] = 2
#font_properties = {'family': 'serif', 'serif': ['Times'], 'weight': 'light'}
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
plt.rcParams["mathtext.default"] = "regular" # for the same font used in regular text.
rc('text', usetex=True)

# colors
# c = ['blueviolet','plum', 'green', 'red']
#n = 2
c_s1 = extract_colors('viridis', 3, portion="middle")
c_s2 = extract_colors('magma', 2, portion="middle", start=0.4, end=0.6)
c = np.concatenate((c_s1, c_s2))

# set up figure
fig, axs = plt.subplots(1, sharex=True)

# set x and y labels
plt.xlabel(r"Radius (pc)", fontsize=20)
#plt.ylabel(r"$\rm Q_{adv} / Q_{rad}$", fontsize=20)
plt.ylabel(r"$\rm Q_{rad} (\rm erg \, s^{-1} \, cm^{-3})$", fontsize=20)
# set minorticks
min_n_index = -3
max_n_index = 1
a1 = np.logspace(min_n_index, max_n_index, np.abs(min_n_index - max_n_index) + 1)
a2 = np.arange(1, 10, 1)
minorticks = np.outer(a1, a2).flatten()
plt.yticks(minorticks, minor=True)

# plot Beckmann cooling line
#plt.loglog(radius_beck, q_ratio_beck, color="grey", label="Beckmann_2018",linewidth=4)

labels = ['1B.b01_1Myr', '1B.b16_1Myr', '1B.m16_1Myr', '2B.b08_0.75Myr', '2B.m08_0.93Myr']
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

    # make ds label
    # label = tidy_data_labels(str(sim[j]) + "_" + "{:.2f}".format(ss_age[0]/1e6) + "Myr")
    # DS.append(ds)
    # labels.append(label)

    # take absolute value of cooling rate
    cooling_ratio = np.abs(profile[("enzo", "radiative_cooling_rate")][profile.used]) 
                    # profile[("enzo", "advective_cooling_rate")][profile.used])
    n = 2
    x_avg = interpolate_data(moving_average(profile.x[profile.used], n=n))
    y_avg = interpolate_data(moving_average(cooling_ratio, n=n))

    #plt.loglog(profile.x[profile.used], cooling_ratio)
    plt.loglog(x_avg, y_avg, color=c[j], label=labels[j])
    #plt.axhline(y=1, color='grey', linestyle='dashed', linewidth=2, alpha=1)
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
        plt.fill_between(x_avg, 0, y_avg.max(), where=x_avg <= disc_r, facecolor=cr[0], transform=trans, alpha=0.5)

plt.title("Radiative Cooling Rate Radial Profile", fontsize=20)

plot_name = 'radial-profile-plot-radiative-cooling-' + str(sim[0]) + "_" + str(sim[1]) + '+s2.pdf'
#+ '_' + str(sim[1]) + '_' + str(sim[2]) 

fig.set_size_inches(8, 5.8)
plt.savefig('plots/' + plot_name, dpi=100)
print("created ", 'plots/' + plot_name)
