import sys
import matplotlib.pyplot as plt
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np

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
plt.rcParams['lines.linewidth'] = 2

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 6,
        }

c1 = ['#d46a7e', 'lightblue', 'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']
# purple/lavender, greeny-blue, cameo green, dusky-pink/darker pink, mustard yellow, classic blues, fuschia/burgundy
c = ['#856798', '#0a888a', '#56ae57', '#ba6873', '#ffc512', '#436bad', '#9d0759', '#9d0216', '#7bc8f6', '#d0c101',
      '#c4387f','#7bb274', '#06b48b', '#6c3461']
l = bhl_object_labels

for i, BHL in enumerate(bhl_object_list):
    BHL.average_times = BHL.average_times/1e6
    BHL.accrate_times = np.array(BHL.accrate_times)/1e6
    BHL.mass_times = np.array(BHL.mass_times) / 1e6
    BHL.hl_times = np.array(BHL.hl_times) / 1e6

    # 1) BH Mass
    axs[0].plot(BHL.mass_times, BHL.mass, color=c[i], linestyle='solid', label=l[i], alpha=0.7)

    # 2) Accretion Rates
    axs[1].plot(BHL.accrate_times, BHL.accrates, color=c[i], linestyle='solid', label=l[i], alpha=0.7)

    # 3) Densities
    axs[2].plot(BHL.average_times, BHL.average_density, color=c[i], linestyle='solid', label=l[i], alpha=0.7)

    # 4) Velocities
    axs[3].plot(BHL.average_times, BHL.average_vinfinity, color=c[i], linestyle='solid', label=l[i]+'-vinf', alpha=0.7)
    axs[3].plot(BHL.average_times, BHL.average_cinfinity, color=c[::-1][i], linestyle='dotted', label=l[i]+'-cinf', alpha=0.7)

    # 5) HL radius
    axs[4].plot(BHL.hl_times, BHL.hl_radius, color=c[i], linestyle='solid', label=l[i], alpha=0.7)


# format plots
axs[0].set_ylabel(r"Black Hole Mass ($\mathrm{M_{\odot}}$)", fontdict=font)
axs[0].legend(loc="lower right", fontsize=6, ncol=1)  # upper/lower
axs[1].set_ylabel(r"Accretion Rate ($\mathrm{M_{\odot}/yr}$)", fontdict=font)
axs[2].set_ylabel(r"Number Density ($\mathrm{cm^{-3}}$)", fontdict=font)
axs[3].set_ylabel(r"Velocities (km/s)", fontdict=font)
axs[4].set_ylabel(r"HL Radius (pc)", fontdict=font)
axs[4].axhline(y=0.015, color=c[::-1][0], linestyle='solid', label="dx = 0.015 pc")
axs[4].axhline(y=0.0031, color=c[1], linestyle='solid', label="dx = 3.1e-03 pc")
axs[4].axhline(y=0.0015, color=c[::-1][2], linestyle='solid', label="dx = 1.5e-03 pc")
axs[4].axhline(y=0.00077, color=c[::-1][3], linestyle='solid', label="dx = 7.7e-04 pc")
axs[4].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)

for i in range(5):
    axs[i].tick_params(axis="x", which='minor', length=4, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=6, width=2, length=4, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=6)
    axs[i].tick_params(axis="y", which='minor', labelsize=6)
    axs[i].set_yscale('log')
    axs[i].set_xlim([0, 0.200]) # for truncated view

axs[4].legend(loc="upper right", fontsize=4, ncol=2)
axs[3].set_yscale('linear')
axs[3].legend(loc="upper right", fontsize=4, ncol=2)  # upper/lower
#axs[2].legend(loc="upper right", fontsize=6, ncol=2)  # upper/lower
#axs[1].legend(loc="upper right", fontsize=6, ncol=2)  # upper/lower
axs[0].legend(loc="upper left", fontsize=6, ncol=2)  # upper/lower

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(4.8, 8)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)
