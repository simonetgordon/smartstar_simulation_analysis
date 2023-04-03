import string

import yt
import sys
import os
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
#from smartstar_find import ss_properties
from read_arrays_from_csv import bhl_object_list, bhl_object_labels

########################################       Plot Residuals     ########################################
#
# to run: python plot_residuals.py [csv1 - baseline] [csv2] [csv3] [output_plotname]
##########################################################################################################

# function that does the interpolation
def interpolate_data(arr, N=1000):
    # interpolate arr over `N` evenly spaced points
    min_val = np.min(arr)
    max_val = np.max(arr)

    t_orig = np.linspace(min_val, max_val, len(arr))
    t_interp = np.linspace(min_val, max_val, N)
    f = interp1d(x=t_orig, y=arr)
    interp_arr = f(t_interp)
    return interp_arr


def first_index(a, val, rtol=0.001, atol=10):
    return next(i for i, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


# def movingaverage(arr, window_size):
#     window = np.ones(int(window_size))/float(window_size)
#     return np.convolve(arr, window, 'valid')


def movingaverage(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

def running_mean(l, N):
    sum = 0
    result = list( 0 for x in l)

    for i in range(0, N):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range(N, len(l)):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result



if __name__ == "__main__":

    # tidy up data labels
    data_labels = [i.replace("-2", "") for i in bhl_object_labels]
    data_labels = [i.replace("RS", "") for i in data_labels]

    # grac accretion scheme from data labels
    if "m" in data_labels[1]:
        acc_scheme = "Mass-Flux"
    else:
        acc_scheme = "BHL"

    # initialise figure
    baseline_sim = sys.argv[1][16:24]
    plot_name = "residuals-" + acc_scheme + "-plot-baseline-" + str(baseline_sim) + ".pdf"
    n_subplots = 3
    fig = plt.figure()
    fig, axs = plt.subplots(n_subplots, 1, sharex=True)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["mathtext.default"] = "regular"
    linewidth = 1.5
    plt.rcParams['lines.linewidth'] = linewidth

    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 12
            }
    fontsize = 12

    c = ['dodgerblue', 'limegreen', 'crimson', 'salmon', 'lightgreen', 'khaki', 'plum', 'seagreen', 'steelblue', 'salmon']

    # define baseline age and accrate lines + number of data points
    i_start = 100 # starting index
    baseline_age_raw = bhl_object_list[0].ages[i_start:]
    n_data_max = baseline_age_raw.shape[0]
    N = int(n_data_max/10)
    baseline_age = interpolate_data(bhl_object_list[0].ages[i_start:], N=N)
    baseline_accrate = interpolate_data(bhl_object_list[0].accrates[i_start:], N=N)

    # iterate over non-baseline data and generate arrays of same size as the baseline
    # using interpolation
    for i, BHL in enumerate(bhl_object_list[1:]):

        # print("original final age: ", BHL.ages[-1])
        # print("target final age: ", baseline_age[-1])

        # True
        len(BHL.ages) == len(BHL.accrates)

        # find index of age that matches end age of baseline
        i_age = first_index(BHL.ages[i_start:], baseline_age[-1], atol=10)

        print("new final age accuracy: ", 1-np.abs(BHL.ages[i] - baseline_age[-1])/baseline_age[-1])

        # truncate the array at this point
        trunc_age = BHL.ages[i_start:i_age]
        trunc_accrate = BHL.accrates[i_start:i_age]

        # try moving average
        window_size = 100
        move_avg_accrate = movingaverage(trunc_accrate, window_size)
        move_avg_age = movingaverage(trunc_age, window_size)

        # interpolate this array over N evenly spaced points
        interp_age = interpolate_data(move_avg_age, N=N)
        interp_age_raw = interpolate_data(trunc_age, N=n_data_max)
        interp_accrate = interpolate_data(move_avg_accrate, N=N)
        interp_accrate_raw = interpolate_data(trunc_accrate, N=n_data_max)

        print("----------------------------------------")

        # define residual
        residual = (interp_accrate - baseline_accrate)*yt.units.msun/yt.units.yr
        residual_normed = np.nan_to_num((interp_accrate - baseline_accrate)/baseline_accrate)

        #axs[0].plot(interp_age/1e6, interp_accrate, color=c[i], label=data_labels[i+1])
        axs[0].plot(interp_age_raw / 1e6, interp_accrate_raw, color=c[i], label=data_labels[i + 1], alpha=0.2)
        axs[0].plot(move_avg_age / 1e6, move_avg_accrate, color=c[i], label=data_labels[i + 1])
        axs[1].plot(interp_age/1e6, residual, color=c[i], label=data_labels[i+1])
        axs[2].plot(interp_age/1e6, residual_normed, color=c[i], label=data_labels[i+1])

    # plot baselines in grey
    axs[0].plot(baseline_age / 1e6, baseline_accrate, color='grey', label=data_labels[0])
    axs[1].axhline(y=0, color='grey', linestyle='--', lw=linewidth, alpha=1)

    # format plots
    axs[0].set_title(acc_scheme + " with " + str(baseline_sim) + " baseline", fontdict=font)
    axs[0].set_ylabel(r"$\dot{M} \, (M_{\odot}/yr)$", fontdict=font)
    axs[0].set_yscale('log')
    axs[1].set_ylabel(r"$\Delta \dot{M} \, (M_{\odot}/yr)$", fontdict=font)
    #axs[0].set_xlabel(r"Time Since Formation (yr)", fontdict=font)
    axs[2].legend()
    axs[2].set_ylabel(r"$\Delta \dot{M} \, / \, \dot{M}_{m16}$", fontdict=font)
    axs[2].set_xlabel(r"Time Since Formation (Myr)", fontdict=font)
    for i in range(n_subplots):
        axs[i].tick_params(axis="x", which='minor', length=4, direction="in")
        axs[i].tick_params(axis="x", which='major', labelsize=fontsize, width=2, length=4, direction="in")
        axs[i].tick_params(axis="y", which='major', labelsize=fontsize)
        axs[i].tick_params(axis="y", which='minor', labelsize=fontsize-2)
        #axs[i].set_yscale('log')
        #axs[i].set_xlim([0, xlim]) # for truncated view

    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(8, 6)
    fig.savefig('plots/' + plot_name, dpi=100)
    print("created ", plot_name)