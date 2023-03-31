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


def first_index(a, val, rtol=0.0001, atol=1):
    return next(i for i, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


if __name__ == "__main__":

    # initialise figure
    plot_name = "residuals-plot-baseline-" + str(sys.argv[1][27:34])
    fig = plt.figure()
    fig, axs = plt.subplots(2, 1, sharex=True)
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

    # define baseline age and accrate lines + number of data points
    baseline_age = bhl_object_list[0].ages
    baseline_accrate = bhl_object_list[0].accrates
    n_data = baseline_age.shape[0]

    # iterate over non-baseline data and generate arrays of same size as the baseline
    # using interpolation
    for i, BHL in enumerate(bhl_object_list[1:]):

        # print("original final age: ", BHL.ages[-1])
        # print("target final age: ", baseline_age[-1])

        # True
        len(BHL.ages) == len(BHL.accrates)

        # find index of age that matches end age of baseline
        i = first_index(BHL.ages, baseline_age[-1], atol=10)

        print("new final age accuracy: ", 1-np.abs(BHL.ages[i] - baseline_age[-1])/baseline_age[-1])

        # truncate the array at this point
        trunc_age = BHL.ages[:i]
        trunc_accrate = BHL.accrates[:i]

        # interpolate this array over N evenly spaced points
        N = n_data
        interp_age = interpolate_data(trunc_age, N=N)
        interp_accrate = interpolate_data(trunc_accrate, N=N)

        print("----------------------------------------")






