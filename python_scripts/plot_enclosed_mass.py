import sys
import matplotlib.pyplot as plt
from smartstar_find import ss_properties
import seaborn as sns
import yt
import os

y = sys.argv[-1] # naming plot

root_dir = ["/work/sc070/sc070/stg/simulations/seed1-bh-only/270msun/replicating-beckmann/",
            "/work/sc070/sc070/stg/simulations/seed2-bh-only/270msun/replicating-beckmann-2/"]
sim = ["1B.RSb01", "2B.RSb01"]
dds = ["DD0128", "DD0197"]
labels = []
DS = []
for i, dd in enumerate(dds):
    ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
    label = "s" + str(i) + "_" + str(float(ds.current_time.to('Myr')))
    DS.append(ds)
    labels.append(label)


fig = plt.figure()
fig, axs = plt.subplots(1, 1, sharex=True)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"
linewidth = 2
plt.rcParams['lines.linewidth'] = linewidth

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }
fontsize = 8

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
        logs={"radius": False},
    )

    mass_plot = axs[0].loglog(rp[0].x[rp[0].used], rp[0][("gas", "mass")][rp[0].used],
                              color=c[j], linestyle='solid', label=labels[i], alpha=alpha)
    j += 2


axs[0].set_xlabel(r"$r \, (pc)$", fontdict=font)
axs[0].set_ylabel(r"$M \, (M_{\odot})$", fontdict=font)

# save plot as pdf
fig = plt.gcf()
# fig.subplots_adjust(wspace=0, hspace=0)
# fig.set_size_inches(4.8, 8)
plot_name = 'radial_profile_mass_enclosed' + str(y) + '.png'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)