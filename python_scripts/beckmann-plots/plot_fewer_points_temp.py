import sys
import matplotlib.pyplot as plt
import BHL_variables as BHL
import numpy as np
from rdp import rdp

x = "s1-" # plot number
y = sys.argv[-1] # naming plot

fig = plt.figure()
fig, axs = plt.subplots(2, 1, sharex=True)
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

acc = np.column_stack((BHL.accrate_times.reshape(-1,1), BHL.accrates.reshape(-1,1)))
print("acc len  = ", len(acc))

mask = rdp(acc, epsilon=0.004, algo="rec")
acc_reduced = mask
print("len acc[mask] is ", len(acc_reduced))

# 1) Accretion Rates
axs[0].plot(acc_reduced[:, 0], acc_reduced[:, 1], color=c[0], linestyle='solid', label=l[0])
#axs[0, 0].set_xlabel(r"r (pc)", fontdict=font)
axs[0].set_ylabel(r"Accretion Rate (Msun/yr)", fontdict=font)
axs[0].set_yscale('log')

for i in range(1):
    axs[i].tick_params(axis="x", which='minor', length = 4, direction="in")
    axs[i].tick_params(axis="x", which='major', labelsize = 8, width=2, length=7, direction="in")
    axs[i].tick_params(axis="y", which='major', labelsize = 8)
    axs[i].tick_params(axis="y", which='minor', labelsize=6)


# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(5, 5)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)