import sys
import matplotlib.pyplot as plt
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
import numpy as np

##########################################################################################################
#                                           Plot BHL Variables
#
# to run: python plot_variables.py [output_plotname]
##########################################################################################################

x = "s1-" # plot number
y = sys.argv[-1] # naming plot

fig = plt.figure()
fig, axs = plt.subplots(4, 1, sharex=True)
plt.rcParams["font.family"] = "serif"
plt.rcParams['lines.linewidth'] = 2

font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 8,
        }

c = ['#d46a7e', 'lightblue', 'lightgreen', 'salmon', 'khaki', 'plum', 'seagreen', 'steelblue']
l = bhl_object_labels

for i, BHL in enumerate(bhl_object_list):
    BHL.average_times = BHL.average_times/1e6
    BHL.accrate_times = np.array(BHL.accrate_times)/1e6

    # 1) Accretion Rates
    axs[0].plot(BHL.accrate_times, BHL.accrates, color=c[i], linestyle='solid', label=l[i])

    # 2) Velocities
    axs[1].plot(BHL.average_times, BHL.average_vinfinity, color=c[i], linestyle='solid', label=l[i])
    axs[1].plot(BHL.average_times, BHL.average_cinfinity, color=c[::-1][i], linestyle='dotted', label=l[i])

    # 3) HL radius
    axs[2].plot(BHL.average_times, BHL.hl_radius, color=c[i], linestyle='solid', label=l[i])
    a = np.empty(len(BHL.average_times))
    a.fill(0.015)
    axs[2].plot(BHL.average_times, a, color='#9a0200', linestyle='solid', label="dx = 0.015 pc")

    # 4) BH Mass
    axs[3].plot(BHL.average_times, BHL.mass, color=c[i], linestyle='solid', label=l[i])

# format plots
axs[0].set_ylabel(r"Accretion Rate (Msun/yr)", fontdict=font)
axs[1].set_ylabel(r"Velocities (km/s)", fontdict=font)
axs[2].set_ylabel(r"Radii (pc)", fontdict=font)
axs[3].set_ylabel(r"Black Hole Mass (Msun)", fontdict=font)
axs[3].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)

for i in range(4):
    axs[i].tick_params(axis="x", which='minor', length=4, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize=8, width=2, length=7, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize=8)
    axs[i].tick_params(axis="y", which='minor', labelsize=6)
    axs[i].set_yscale('log')
    axs[i].legend(loc="upper left", fontsize=8, ncol=2)  # upper/lower

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(5, 10)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)
