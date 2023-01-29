import sys
import matplotlib.pyplot as plt
import BHL_variables as BHL
import numpy as np

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

font2 = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 14,
        }
c = ['#d46a7e', 'lightblue', 'lightgreen', 'salmon', 'khaki', 'plum', 'seagreen', 'steelblue']
l = ['Accretion Rate', 'vInfinity', 'cInfinity', 'Hoyle-Lyttleton Radius', 'Cell Width', 'BH Mass']

BHL.avg_times = BHL.avg_times/1e6
BHL.accrate_times = np.array(BHL.accrate_times)/1e6
# 1) Accretion Rates
axs[0].plot(BHL.accrate_times, BHL.accrates, color=c[0], linestyle='solid', label=l[0])
#axs[0, 0].set_xlabel(r"r (pc)", fontdict=font)
axs[0].set_ylabel(r"Accretion Rate (Msun/yr)", fontdict=font)

for i in range(4):
    axs[i].tick_params(axis="x", which='minor', length = 4, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize = 8, width=2, length=7, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize = 8)
    axs[i].tick_params(axis="y", which='minor', labelsize=6)

# 2) Velocities
axs[1].plot(BHL.avg_times, BHL.avg_vinfinities, color=c[1], linestyle='solid', label=l[1])
axs[1].plot(BHL.avg_times, BHL.avg_cinfinities, color=c[2], linestyle='solid', label=l[2])
axs[1].set_ylabel(r"Velocities (km/s)", fontdict=font)

# 3) HL radius
a = np.empty(len(BHL.avg_times))
a.fill(0.015)
axs[2].plot(BHL.avg_times, BHL.hl_radii, color=c[3], linestyle='solid', label=l[3])
axs[2].plot(BHL.avg_times, a, color=c[4], linestyle='solid', label=l[4])
axs[2].set_ylabel(r"Radii (pc)", fontdict=font)

# 4) BH Mass
axs[3].plot(BHL.avg_times, BHL.bh_masses, color=c[5], linestyle='solid', label=l[5])
#axs[3].set_xscale('log')
axs[3].set_ylabel(r"Black Hole Mass (Msun)", fontdict=font)
axs[3].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)
for i in range(4):
    #axs[i].set_xscale('log')
    axs[i].set_yscale('log')
    axs[i].legend(loc="upper left", fontsize=8,  ncol=2) # upper/lower

# axs[3].tick_params(axis="x", which='minor', length=4, direction="in")
# axs[3].tick_params(axis="x", which='major', labelsize=8, width=2, length=7, direction="in")
# axs[3].tick_params(axis="y", which='major', labelsize=8)
#axs[3].yaxis.tick_right()
#axs[3].yaxis.set_label_position("left")

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(5, 10)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)