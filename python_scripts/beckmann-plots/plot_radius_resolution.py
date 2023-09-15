import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import sys
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import tidy_data_labels, first_index, interpolate_data, movingaverage

##########################################################################################################
#                                     Plot radius resolution vs time
#
# to run: python plot_radius_resolution.py [csv1] [csv2] [csv3] [output_plotname e.g MF-BHL]
# for 2x2 update: list MF runs first, then BHL runs. Name: MF+BHL
##########################################################################################################

# set x-axis extent in Myr and simulation set
xlim = 0.4 #0.225 - 0.4 for 1S, 1 for 1B
sim = "s1-10.8msun-"
atol = 1e-4 # 7e-2 for 1B.m, 1e-4 otherwise
accretion = sys.argv[-1] # naming plot with accretion scheme

# set y axis limits (optional)
min_hl = max_hl = min_bondi = max_bondi = 1
# min_hl = 0.08
max_hl = 50
# min_bondi = 1
# max_bondi = 6000

# font settings
linewidth = 1.5
fontsize = 12
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams['font.size'] = fontsize
plt.rcParams['lines.linewidth'] = linewidth

# initialise figure
fig = plt.figure()
num_subplots = 2
fig, axs = plt.subplots(num_subplots, 2, sharex=True)
c = ['indigo', 'blueviolet', 'turquoise', 'limegreen', 'darkgreen'] # line colours

# tidy data labels
l = tidy_data_labels(bhl_object_labels)

# allocate cell widths to be shown based on input
if sim == "s1-270msun-":
    dx = [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04]
elif sim == "s1-10.8msun-":
    dx = [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04]
elif sim == "s2-270msun-":
    dx = [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04]
elif sim == "s2-10.8msun-":
    dx = [0.00832, 0.00416, 0.002079909, 0.001297054, 0.0005188221, 0.0001297054]
elif sim == "s2-10.8msun-2-":
    dx = [1.3e-03, 5.2e-04, 1.3e-04]

# parameters
alpha = 0.9
time_cutoff = xlim  # Myrs
i_start = 1
window_size = 10
i_mf = np.arange(int(len(bhl_object_list)/2))
i_bhl = np.arange(int(len(bhl_object_list)/2), len(bhl_object_list))

for i, BHL in enumerate(bhl_object_list):
    # convert ages from yrs to Myrs
    BHL.ages = np.array(BHL.ages) / 1e6

    # find index of age that matches end age of time limit
    i_age = first_index(BHL.ages[i_start:], time_cutoff, rtol=1e-8, atol=atol)

    # calculate age and radius moving averages
    age = movingaverage(BHL.ages[i_start:i_age], window_size)
    hl_radius = movingaverage(BHL.hl_radius[i_start:i_age], window_size)
    bondi_radius = movingaverage(BHL.bondi_radius[i_start:i_age], window_size)

    if i in i_mf:

        # calculate how many cells it's resolving the hl radius by
        dx_res_hl = hl_radius / dx[i]
        dx_res_bondi = bondi_radius / dx[i]

        # 1) HL radius resolution in cell widths
        axs[0, 0].plot(age, dx_res_hl, color=c[i], linestyle='solid', label=l[i], alpha=alpha)

        # 2) Bondi radius resolution in cell widths
        axs[1, 0].plot(age, dx_res_bondi, color=c[i], linestyle='solid', label=l[i], alpha=alpha)
    else:

        # calculate how many cells it's resolving the hl radius by
        dx_res_hl = hl_radius / dx[i-int(len(bhl_object_list)/2)]
        dx_res_bondi = bondi_radius / dx[i-int(len(bhl_object_list)/2)]

        # 1) HL radius resolution in cell widths
        axs[0, 1].plot(age, dx_res_hl, color=c[i-int(len(bhl_object_list)/2)], linestyle='solid', label=l[i], alpha=alpha)

        # 2) Bondi radius resolution in cell widths
        axs[1, 1].plot(age, dx_res_bondi, color=c[i-int(len(bhl_object_list)/2)], linestyle='solid', label=l[i], alpha=alpha)

    # update y-axis limits
    min_hl = min(dx_res_hl.min(), min_hl)
    max_hl = max(dx_res_hl.max(), max_hl)
    min_bondi = min(dx_res_bondi.min(), min_bondi)
    max_bondi = max(dx_res_bondi.max(), max_bondi)


################################ Format Plot ####################################

# label axes
for i in range(num_subplots):
    for j in range(num_subplots):
        axs[i, j].set_xlabel(r"BH Age (Myr)", fontdict=None)

        # x-axis ticks
        axs[i, j].set_xticks(np.arange(0, time_cutoff-0.02, 0.1))
        axs[i, j].minorticks_on()
        axs[i, j].xaxis.set_minor_locator(plt.MultipleLocator(0.02))
        axs[i, j].tick_params(axis="x", which='minor', length=2, direction="in")
        axs[i, j].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=4, direction="in")
        axs[i, j].set_yscale('log')
        axs[i, j].legend()

        # set no cells = 1 line
        axs[i, j].axhline(y=1, color="grey", linestyle='dashdot', lw=linewidth, alpha=alpha)

# format individual subplots
axs[0, 0].set_ylabel(r"$R_{\rm HL}$ Resolution (cell widths)", fontdict=None)
axs[0, 0].set_title(str(sim) + str(accretion)[:2], fontdict=None)
axs[0, 0].set_ylim(min_hl - min_hl*0.1, max_hl + max_hl*0.2)
axs[1, 0].set_ylabel(r"$R_{\rm Bondi}$ Resolution (cell widths)", fontdict=None)
axs[1, 0].set_ylim(min_bondi - min_bondi*0.1, max_bondi + max_bondi*0.25)
axs[0, 1].set_title(str(sim) + str(accretion)[3:6], fontdict=None)
axs[0, 1].tick_params(axis="y", which='minor', length=2, direction="in")
axs[0, 1].tick_params(axis="y", which='major', width=1, length=4, direction="in")
axs[0, 1].set_yticklabels([])
axs[1, 1].tick_params(axis="y", which='minor', length=2, direction="in")
axs[1, 1].tick_params(axis="y", which='major', width=1, length=4, direction="in")
axs[1, 1].set_yticklabels([])
axs[0, 0].set_ylim(min_hl - min_hl*0.1, max_hl + max_hl*0.2)
axs[0, 1].set_ylim(min_hl - min_hl*0.1, max_hl + max_hl*0.2)
axs[1, 0].set_ylim(min_bondi - min_bondi*0.1, max_bondi + max_bondi*0.25)
axs[1, 1].set_ylim(min_bondi - min_bondi*0.1, max_bondi + max_bondi*0.25)

# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0.02, hspace=0.02)
fig.set_size_inches(10, 7)
plot_name = 'time-dx_res-' + str(sim) + str(accretion) + '.pdf'
fig.savefig('plots/' + plot_name, bbox_inches='tight')
print("created plots/", plot_name)

