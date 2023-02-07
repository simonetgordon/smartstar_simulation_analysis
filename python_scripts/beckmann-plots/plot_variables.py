import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np
import seaborn as sns

##########################################################################################################
#                                           Plot BHL Variables
#
# to run: python plot_variables.py [csv1] [csv2] [csv3] [output_plotname]
##########################################################################################################

x = "s1-" # plot number
y = sys.argv[-1] # naming plot

fig = plt.figure()
fig, axs = plt.subplots(5, 1, sharex=True)
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

c1 = ['#d46a7e', 'lightblue', 'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']
# purple/lavender, greeny-blue, cameo green, dusky-pink/darker pink, mustard yellow, classic blues, fuschia/burgundy
c2 = ['#856798', '#0a888a', '#56ae57', '#ba6873', '#ffc512', '#436bad', '#9d0759', '#9d0216', '#7bc8f6', '#d0c101',
      '#c4387f','#7bb274', '#06b48b', '#6c3461']
c = sns.color_palette("Paired", len(bhl_object_list*2))

l = bhl_object_labels
j = 1
alpha = 0.7
for i, BHL in enumerate(bhl_object_list):
    BHL.average_times = BHL.average_times/1e6
    BHL.accrate_times = np.array(BHL.accrate_times)/1e6
    BHL.mass_times = np.array(BHL.mass_times) / 1e6
    BHL.hl_times = np.array(BHL.hl_times) / 1e6

    # 1) BH Mass
    axs[0].plot(BHL.mass_times, BHL.mass, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 2) Accretion Rates
    axs[1].plot(BHL.accrate_times, BHL.accrates, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 3) Densities
    axs[2].plot(BHL.average_times, BHL.average_density, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 4) Velocities
    axs[3].plot(BHL.average_times, BHL.average_vinfinity, color=c[j], linestyle='solid', label=l[i]+'-vinf', alpha=alpha)
    axs[3].plot(BHL.average_times, BHL.average_cinfinity, color=c[j-1], linestyle='dotted', label=l[i]+'-cinf', alpha=alpha)

    # 5) HL radius
    axs[4].plot(BHL.hl_times, BHL.hl_radius, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    j += 2

xlim = 1

for i in range(5):
    axs[i].tick_params(axis="x", which='minor', length=4, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=4, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
    axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)
    axs[i].set_yscale('log')
    axs[i].set_xlim([0, xlim]) # for truncated view

# format plots
axs[0].set_ylabel(r"$M_{BH} \, (M_{\odot})$", fontdict=font)
axs[1].set_ylabel(r"$\dot{M} \, (M_{\odot}/yr)$", fontdict=font)
axs[2].set_ylabel(r"$n \, (H cm^{-3})$", fontdict=font)
axs[3].set_ylabel(r"$\nu \, (km/s)$", fontdict=font)
axs[4].set_ylabel(r"$r_{HL} \, (pc)$", fontdict=font)
axs[3].set_yscale('linear')
axs[4].axhline(y=0.015, color=c[0], linestyle='solid', label="dx = 0.015 pc")
axs[4].axhline(y=0.0031, color=c[2], linestyle='solid', label="dx = 3.1e-03 pc")
axs[4].axhline(y=0.0015, color=c[4], linestyle='solid', label="dx = 1.5e-03 pc")
axs[4].axhline(y=0.00077, color=c[6], linestyle='solid', label="dx = 7.7e-04 pc")
axs[4].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)

if xlim == 0.2:
    axs[0].set_ylim([250, 2500])
    axs[1].set_ylim([1e-4, 4e-2])
if xlim == 1:
    axs[0].set_ylim([250, 2500])
    axs[1].set_ylim([7e-6, 2e-2])
    axs[2].set_ylim([8e4, 2e8])
    axs[3].set_ylim([0, 25])


# legends
dx_lines = [Line2D([0], [0], color=c[0], lw=linewidth),
            Line2D([0], [0], color=c[2], lw=linewidth),
            Line2D([0], [0], color=c[4], lw=linewidth),
            Line2D([0], [0], color=c[6], lw=linewidth),]
dx = ["dx = 0.015 pc", "dx = 3.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
axs[4].legend(dx_lines, dx, loc="upper right", fontsize=6, ncol=2)
vel_lines = [Line2D([0], [0], color='grey', lw=linewidth),
            Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
axs[3].legend(vel_lines, ["vinf", "cinf"], loc="upper left", fontsize=8, ncol=2)  # upper/lower
#axs[2].legend(loc="upper right", fontsize=6, ncol=2)  # upper/lower
#axs[1].legend(loc="upper right", fontsize=6, ncol=2)  # upper/lower
axs[0].legend(loc="lower right", fontsize=8, ncol=2)  # upper/lower

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(4.8, 8)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)
