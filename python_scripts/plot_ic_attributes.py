import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from smartstar_find import ss_properties
from derived_fields import add_fields_ds
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
    add_fields_ds(ds)
    j = i + 1
    label = "s" + str(j) + "_" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"
    DS.append(ds)
    labels.append(label)


fig = plt.figure()
fig, axs = plt.subplots(5, 1, sharex=True)
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
    ss_pos = ss_properties(ds)[0]
    sp = ds.sphere(ss_pos, (1, "kpc"))

    rp = yt.create_profile(
        sp,
        "radius",
        [("gas", "mass"), ("gas", "temperature"), ("gas", "dynamical_time"), ("gas", "baryon_overdensity_sg")],
        units={"radius": "pc", ("gas", "mass"): "Msun"},
        logs={"radius": True, ("gas", "mass"): True, ("gas", "temperature"): True},
        weight_field=None
    )

    rp2 = yt.create_profile(
        sp,
        "radius",
        [("gas", "temperature"), ("gas", "mass"), ("gas", "dynamical_time"), ("gas", "baryon_overdensity_sg"),
         ("gas", "averaged_density")],
        units={"radius": "pc", ("gas", "mass"): "Msun"},
        logs={"radius": True, ("gas", "temperature"): True}
    )

    axs[0].loglog(rp.x.value, rp[("gas", "mass")].value.cumsum(),
               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                  color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "dynamical_time")][rp2.used].to("yr"),
                  color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    # Omega_b = 0.0449
    axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "baryon_overdensity_sg")][rp2.used],
                  color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    axs[4].loglog(rp2.x[rp2.used], rp2[("gas", "averaged_density")][rp2.used],
                  color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    # t_dyn is larger when using non-mass-weighted densities
    # axs[2].loglog(rp.x[rp.used], rp[("gas", "dynamical_time")][rp.used].to("yr"),
    #               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    j += 2

for i in range(5):
    axs[i].tick_params(bottom=True, left=True)
    axs[i].minorticks_on()
    axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
    axs[i].tick_params(axis="y", which='minor', labelsize=fontsize)
    axs[i].grid(color='r', linestyle='-', linewidth=2)

r_lines = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth),
            Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
#axs[3].legend(r_lines, [r"$r_{200}$"], loc="upper right", fontsize=6, ncol=1)  # upper/lower

axs[4].set_xlabel(r"$r \, (pc)$", fontdict=font)
axs[4].set_ylabel(r"$\rho_{avg} \, (g \, cm^{-3})$", fontdict=font)
axs[3].set_ylabel(r"$\delta_b$", fontdict=font)
axs[3].axhline(y=200, color='grey', linestyle='dashed', lw=linewidth, alpha=1)
#axs[3].axhline(y=500, color='grey', linestyle='dashed', lw=linewidth, alpha=1)
for i in range(3):
    axs[i].axvline(x=70, color='grey', linestyle='dashed', lw=linewidth, alpha=1, label=r"$r_{200}$")
axs[2].set_ylabel(r"$t_{dyn} \, (yr)$", fontdict=font)
axs[0].set_ylabel(r"$M \, (M_{\odot})$", fontdict=font)
axs[1].set_ylabel(r"$T \, (K)$", fontdict=font)
axs[0].legend(loc="lower right", fontsize=fontsize, ncol=2)  # upper/lower
axs[0].set_title("Gas properties at time of BH formation", fontdict=font)

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(6, 8)
plot_name = 'radial_profile_ics_halo_1kpc.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created plots/" + str(plot_name))