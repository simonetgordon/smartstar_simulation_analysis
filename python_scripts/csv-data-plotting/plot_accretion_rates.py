import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
import sys
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import tidy_data_labels, first_index, movingaverage

"""
Produce all accretion rates for 1S group
python -i plot_accretion_rates.py [csv1] [csv2] [csv3] [output_plotname e.g mass-flux-x4]
List MF first, then BHL
BEWARE: added 1st row of data-1S.RSm04.csv to data-1S.RSmf8.csv and data-1S.RSmf4.csv. -> need to all start from 10.8 msun.
"""


# text format
linewidth = 2
fontsize = 14
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)
plt.rcParams["mathtext.default"] = "regular"
plt.rcParams['font.size'] = fontsize
plt.rcParams['lines.linewidth'] = linewidth

# line colours
c = ['blueviolet', 'turquoise', 'limegreen', 'darkgreen']

# set up figure
fig = plt.figure()
num_subplots = 2
fig, axs = plt.subplots(1, num_subplots, sharey=True)

# misc parameters
l = tidy_data_labels(bhl_object_labels)
l = ['10.8msun.UR8', '10.8msun.UR4', '10.8msun.R01', '10.8msun.R04', '10.8msun.UR8', '10.8msun.UR4', '10.8msun.R01', '10.8msun.R04']
j = 0
alpha = 0.9
time_cutoff = 0.42 # Myrs
i_start = 0
window_size = 1

# for MF and BHL separation into subplots
i_mf = np.arange(int(len(bhl_object_list)/2))
i_bhl = np.arange(int(len(bhl_object_list)/2), len(bhl_object_list))

for i, BHL in enumerate(bhl_object_list):

    # convert ages from yrs to Myrs
    BHL.ages = np.array(BHL.ages) / 1e6

    # find index of age that matches end age of time limit
    i_age = first_index(BHL.ages[i_start:], time_cutoff, rtol=1e-5, atol=4e-3) 

    # set age and mass 
    age = movingaverage(BHL.ages[i_start:i_age], window_size)
    mass =  movingaverage(BHL.mass[i_start:i_age], window_size)

    # 1) plot BH mass
    if i in i_mf:
        axs[0].plot(age, mass, color=c[i-int(len(bhl_object_list)/2)], linestyle='solid', label=l[i], alpha=alpha)

    else:
        axs[1].plot(age, mass, color=c[i-int(len(bhl_object_list)/2)], linestyle='solid', label=l[i], alpha=alpha)


# label axes
for j in range(num_subplots):
    # set same ylim for both
    axs[j].set_ylim(9.5, 2000)

    # x-axis ticks
    axs[j].set_xticks(np.arange(0, time_cutoff, 0.1))
    axs[j].axhline(20, linestyle='dashed', color='grey', alpha=0.5)
    axs[j].minorticks_on()
    axs[j].xaxis.set_minor_locator(plt.MultipleLocator(0.02))
    axs[j].tick_params(axis="x", which='minor', length=2, direction="in")
    axs[1].tick_params(axis="y", which='minor', length=2, direction="in")
    axs[1].tick_params(axis="y", which='major', width=1, length=4, direction="in")
    axs[j].tick_params(axis="x", which='major', labelsize=fontsize, width=1, length=4, direction="in")
    axs[j].set_yscale('log')
    

# format individual subplots
axs[0].set_title("Mass-Flux Accretion", fontdict=None)
axs[0].set_xlabel("BH Age (Myr)", fontdict=None)
axs[1].set_title("BHL Accretion", fontdict=None)
axs[1].set_xlabel("BH Age (Myr)", fontdict=None)
axs[1].text(0.26, 270, 'High Resolution', fontsize=12, color=c[-1])
axs[1].text(0.26, 1200, 'Low Resolution', fontsize=12, color=c[0])
# Add an arrow
arrow_props = dict(facecolor='red', arrowstyle='->')
axs[1].annotate('', xy=(0.34, 380), xytext=(0.34, 1100),
            arrowprops=arrow_props)
# set  y-label for both
axs[0].set_ylabel(r"BH Mass ($\rm M_\odot$)", fontdict=None)
axs[1].legend(loc='upper left')


# save plot as pdf
fig = plt.gcf()
fig.subplots_adjust(wspace=0.02, hspace=0.15)
fig.set_size_inches(9, 4)
plot_name = 'time-bhmass-1S.png'
fig.savefig('plots/' + plot_name, bbox_inches='tight')
plot_name = 'plots/time-bhmass-1S.pdf'
fig.savefig(plot_name, bbox_inches='tight')
print("created ", plot_name)
