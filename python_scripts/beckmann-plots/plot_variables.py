import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from read_arrays_from_csv_2 import bhl_object_list, bhl_object_labels
import numpy as np
import seaborn as sns

##########################################################################################################
#                                           Plot BHL Variables
#
# to run: python plot_variables_2.py [csv1] [csv2] [csv3] [output_plotname]
##########################################################################################################

x = "s1-270msun-" # plot number
y = sys.argv[-1] # naming plot
xlim = 3

fig = plt.figure()
fig, axs = plt.subplots(6, 1, sharex=True)
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.default"] = "regular"
linewidth = 1.5
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
#c = sns.color_palette("Paired", len(bhl_object_list*8))
c0 = ['#9d0216', '#E66101', '#436bad', '#A6611A']
#c = sns.diverging_palette(140, 300, s=100, l=50, sep=5, center="light", n=len(bhl_object_list*2))
#c = sns.color_palette("spectral", len(bhl_object_list*2), desat=1)
c = sns.blend_palette(["#2CA02C", '#17BECF', "#e377c2", "#76549A"], len(bhl_object_list*2))

l = bhl_object_labels
j = 1
alpha = 0.9
for i, BHL in enumerate(bhl_object_list):
    # convert ages from yrs to Myrs
    BHL.ages = np.array(BHL.ages) / 1e6

    # 1) BH Mass
    axs[0].plot(BHL.ages, BHL.mass, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 2) Accretion Rates
    axs[1].plot(BHL.ages, BHL.accrates, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 3) Densities
    axs[2].plot(BHL.ages, BHL.average_density, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 4) Velocities
    axs[3].plot(BHL.ages, BHL.average_vinfinity, color=c[j], linestyle='solid', label=l[i]+'-vinf', alpha=alpha)
    axs[3].plot(BHL.ages, BHL.average_cinfinity, color=c[j], linestyle='dotted', label=l[i]+'-cinf', alpha=alpha)

    # 5) HL radius
    axs[4].plot(BHL.ages, BHL.hl_radius, color=c[j], linestyle='solid', label=l[i], alpha=alpha)

    # 6) Jeans length
    axs[5].plot(BHL.ages, BHL.jeans_length, color=c[j], linestyle='solid', label=l[i]+'-jeans-length', alpha=alpha)

    j += 2

for i in range(6):
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
axs[5].set_ylabel(r"$r_{jeans} \, (pc)$", fontdict=font)
axs[3].set_yscale('linear')
axs[0].set_yscale('linear')
axs[0].set_title(str(x) + str(y), fontdict=font)
for i in [4, 5]:
    if x == "s1-270msun-":
        dx = [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04]
    elif x == "s1-40msun-":
        dx = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
    elif x == "s2-270msun-":
        dx = [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
    elif x == "s2-40msun-":
        dx = [8.4e-03, 4.2e-03, 2.1e-03]
    elif x == "s2-40msun-2-":
        dx = [1.3e-03, 5.2e-04, 1.3e-04]
    c1 = c2 = c3 = 'black'
    l1 = 'dashed'
    l2 = 'dashdot'
    l3 = 'dotted'
    axs[i].axhline(y=dx[0], color=c1, linestyle=l1, lw=linewidth*0.5,  label="dx = " + str(dx[0]) + "pc", alpha=1)
    axs[i].axhline(y=dx[1], color=c2, linestyle=l2, lw=linewidth*0.5, label="dx = " + str(dx[1]) + "pc", alpha=1)
    axs[i].axhline(y=dx[2], color=c3, linestyle=l3, lw=linewidth*0.5, label="dx = " + str(dx[2]) + "pc", alpha=1)
    #axs[i].axhline(y=0.00077, color=c[6], linestyle='solid', label="dx = 7.7e-04 pc")
axs[5].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)

if xlim == 0.2:
    axs[0].set_ylim([250, 2500])
    axs[1].set_ylim([1e-4, 2e-2])
if xlim == 3:
    axs[0].set_ylim([0, 8000])
    #axs[1].set_ylim([5e-5, 2e-2])
    #axs[2].set_ylim([8e3, 3e8])
    #axs[3].set_ylim([0, 25])


# legends
dx_lines = [Line2D([0], [0], color=c1, lw=linewidth, linestyle=l1),
            Line2D([0], [0], color=c2, lw=linewidth, linestyle=l2),
            Line2D([0], [0], color=c3, lw=linewidth, linestyle=l3),
            #Line2D([0], [0], color=c0[3], lw=linewidth, linestyle='dashed'),
            ]
if x == "s1-270msun-":
    dx = ["dx = 0.012 pc", "dx = 3.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
elif x == "s1-40msun-":
    dx = ["dx = 2.5e-02 pc", "dx = 1.2e-02 pc", "dx = 3.1-03 pc", "dx = 7.7e-04 pc"]
elif x == "s2-270msun-":
    dx = ["dx = 8.3e-03 pc", "dx = 2.1e-03 pc", "dx = 1.5e-03 pc", "dx = 7.7e-04 pc"]
elif x == "s2-40msun-":
    dx = ["dx = 8.4e-03 pc", "dx = 4.2e-03 pc", "dx = 2.1e-03 pc"]
elif x == "s2-40msun-2":
    dx = ["dx = 1.3e-03 pc", "dx = 5.2e-04 pc", "dx = 1.3e-04  pc"]
axs[4].legend(dx_lines, dx, loc="upper right", fontsize=6, ncol=2)
vel_lines = [Line2D([0], [0], color='grey', lw=linewidth),
            Line2D([0], [0], color='grey', linestyle='dotted', lw=linewidth)]
axs[3].legend(vel_lines, ["vinf", "cinf"], loc="upper left", fontsize=6, ncol=2)  # upper/lower
axs[0].legend(loc="lower right", fontsize=6, ncol=2)  # upper/lower

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0, hspace=0)
fig.set_size_inches(4.8, 8)
plot_name = 'time-' + str(x) + str(y) + '.pdf'
fig.savefig('plots/' + plot_name, dpi=100)
print("created ", plot_name)
