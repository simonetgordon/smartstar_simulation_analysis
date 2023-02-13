import sys
import matplotlib.pyplot as plt
from smartstar_find import ss_properties
import numpy as np
import seaborn as sns
import yt
import os

y = sys.argv[-1] # naming plot

root_dir = ["~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
            "~/disk14/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann/"]
sim = ["1B.RSb01", "2B.RSb1"]
dds = ["DD0128/DD0128", "DD0198/DD0198"]
labels = []
DS = []
for i, dd in enumerate(dds):
    ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
    j = i + 1
    label = "s" + str(j) + "_" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"
    DS.append(ds)
    labels.append(label)


fig = plt.figure()
fig, axs = plt.subplots(1, 1, sharex=True)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"
linewidth = 2
plt.rcParams['lines.linewidth'] = linewidth

fontsize = 10
font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': fontsize,
        }

c = sns.color_palette("Paired", len(DS*2))
j = 1
alpha = 0.7
for i, ds in enumerate(DS):

    ss_pos = ss_properties(ds, 0.2)[0]

    sp = ds.sphere(ss_pos, (10, "kpc"))

    rp = yt.create_profile(
        sp,
        "radius",
        [("gas", "mass")],
        units={"radius": "pc", ("gas", "mass"): "Msun"},
        logs={"radius": True, ("gas", "mass"): True},
        weight_field=None
    )

    # axs[0].loglog(rp[0].x[rp[0].used], rp[0][("gas", "mass")][rp[0].used],
    #               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    axs.loglog(rp.x.value, rp[("gas", "mass")].value.cumsum(),
               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    j += 2

axs.tick_params(bottom=True, left=True)
axs.tick_params(axis="x", which='minor', length=4, direction="in")
axs.tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=4, direction="in")
axs.tick_params(axis="y", which='major', labelsize=fontsize)
axs.tick_params(axis="y", which='minor', labelsize=fontsize-2)
axs.set_xlabel(r"$r \, (pc)$", fontdict=font)
axs.set_ylabel(r"$M \, (M_{\odot})$", fontdict=font)

axs.legend(loc="lower right", fontsize=fontsize, ncol=2)  # upper/lower
axs.set_title("Mass enclosed at time of BH formation", fontdict=font)

# save plot as pdf
fig = plt.gcf()
# fig.subplots_adjust(wspace=0, hspace=0)
# fig.set_size_inches(4.8, 8)
plot_name = 'radial_profile_mass_enclosed' + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created plots/" + str(plot_name))