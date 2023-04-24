########################################       Plot IC Attributes     ####################################
#
# to run: python plot_ic_attributes.py radial_profile_ics_halo_1kpc.pdf
##########################################################################################################


import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from smartstar_find import ss_properties
from derived_fields import add_fields_ds
import numpy as np
import yt
import os
from matplotlib import rc

y = sys.argv[-1] # naming plot

root_dir = ["~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
            "~/disk14/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann/"]
sim = ["1B.RSb01-2", "2B.RSb1"]
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

# font settings
fontsize = 12 # for projection annotations
linewidth = 2
plt.rcParams['font.size'] = fontsize
plt.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)

plt.rcParams["mathtext.default"] = "regular"
plt.rcParams['lines.linewidth'] = linewidth

# colours for plotting
c = ['firebrick', 'royalblue', 'crimson', 'dodgerblue', 'limegreen', 'salmon', 'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']

fig = plt.figure()
n_subplots = 4
fig, axs = plt.subplots(n_subplots, 1, sharex=True)

# create profiles
j = 1
alpha = 0.7
for i, ds in enumerate(DS):
    ss_pos = ss_properties(ds)[0]
    sp = ds.sphere(ss_pos, (1, "kpc"))

    rp = yt.create_profile(
        sp,
        "radius",
        [("gas", "mass"), ("gas", "temperature"), ("gas", "dynamical_time"), ("gas", "baryon_overdensity_sg"),
         ("gas", "blackhole_freefall_timescale"), ("gas", "theta_vel_dynamical_timescale"), ("deposit", "all_mass") ],
         units={"radius": "pc", ("gas", "mass"): "Msun", ("deposit", "all_mass"): "Msun"},
        logs={"radius": True, ("gas", "mass"): True, ("gas", "temperature"): True, ("deposit", "all_mass"):True},
        weight_field=None
    )

    rp2 = yt.create_profile(
        sp,
        "radius",
        [("gas", "temperature"), ("gas", "mass"), ("gas", "dynamical_time"), ("gas", "baryon_overdensity_sg"),
         ("gas", "number_density"), ("gas", "blackhole_freefall_timescale"), ("gas", "theta_vel_dynamical_timescale")],
        units={"radius": "pc", ("gas", "mass"): "Msun"},
        logs={"radius": True, ("gas", "temperature"): True}
    )

    # plot data
    axs[0].loglog(rp.x.value, rp[("gas", "mass")].value.cumsum(),
               color=c[i], linestyle='solid', label=labels[i] + "_gas", alpha=alpha)
    axs[0].loglog(rp.x.value, rp[("deposit", "all_mass")].value.cumsum(),
            color=c[i+2], linestyle='solid', label= labels[i] + "_DM", alpha=alpha)
    r200 = rp[("deposit", "all_mass")].value.cumsum().max() + rp[("gas", "mass")].value.cumsum().max()
    print("r200 mass: ", r200)

    axs[1].loglog(rp2.x[rp2.used], rp2[("gas", "temperature")][rp2.used],
                  color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

    axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "blackhole_freefall_timescale")][rp2.used].to("Myr"),
                  color=c[i], linestyle='solid', label=r"BH ff time $\propto M$", alpha=alpha)
    # axs[2].loglog(rp2.x[rp2.used], rp2[("gas", "theta_vel_dynamical_timescale")][rp2.used].to("yr"),
    #               color=c[j], linestyle='solid', label=r"$2\pi r / v_\theta $", alpha=alpha)

    # Omega_b = 0.0449
    # axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "baryon_overdensity_sg")][rp2.used],
    #               color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

    axs[3].loglog(rp2.x[rp2.used], rp2[("gas", "number_density")][rp2.used],
                  color=c[i], linestyle='solid', label=labels[i], alpha=alpha)

    # t_dyn is larger when using non-mass-weighted densities
    # axs[2].loglog(rp.x[rp.used], rp[("gas", "dynamical_time")][rp.used].to("yr"),
    #               color=c[j], linestyle='solid', label=labels[i], alpha=alpha)

    j += 2

# set ticks
for i in range(n_subplots):
    axs[i].tick_params(bottom=True, left=True)
    axs[i].minorticks_on()
    axs[i].tick_params(axis="x", which='minor', length=2, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=3, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
    axs[i].tick_params(axis="y", which='minor', length=2)
    #axs[i].grid(color='grey', linestyle='solid', linewidth=0.5, alpha=0.7)

# make lines for legend
r_lines = [Line2D([0], [0], color='grey', linestyle='dashed', lw=linewidth),
            Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]

# set axis labels
axs[n_subplots-1].set_xlabel(r"$\rm Radius \, (pc)$", fontdict=None)
axs[n_subplots-1].set_xlim(8e-3, 1.5e3)
xticks = np.logspace(-2, 3, 6)
axs[3].set_xticks(xticks)
axs[3].set_ylabel(r"$\rm n \, (cm^{-3})$", fontdict=None)
#axs[3].set_ylabel(r"$\rm \delta_b$", fontdict=None)
#axs[3].axhline(y=200, color='grey', linestyle='dashed', lw=linewidth, alpha=1)
for i in range(n_subplots):
    axs[i].axvline(x=70, color='grey', linestyle='dashed', lw=linewidth, alpha=1, label=r"$r_{200}$")
axs[2].set_ylabel(r"$\rm t_{ff} \, (Myr)$", fontdict=None)
#axs[2].set_ylim(200, rp2[("gas", "blackhole_freefall_timescale")][rp2.used].to("Myr").d.max()+100)
axs[1].set_ylabel(r"$\rm T \, (K)$", fontdict=None)
axs[0].set_ylabel(r"$\rm M_{encl} \, (M_{\odot})$", fontdict=None)
axs[0].legend(loc="lower right", fontsize=fontsize-2, ncol=2)  # upper/lower
#axs[0].set_title("Gas properties at time of BH formation", fontdict=None)

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(4.6, 6.2)
plot_name = 'radial_profile_ics_halo_1kpc.pdf'
fig.savefig('plots/' + plot_name, bbox_inches='tight')
print("created plots/" + str(plot_name))