import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import stats
import numpy as np
import sys
from read_arrays_from_csv import bhl_object_list, bhl_object_labels
from plot_variables import tidy_data_labels, first_index, interpolate_data, moving_average, eddington_rate
from matplotlib.ticker import FixedLocator

# Define the line fit function
def log_line_fit(x, slope, intercept):
    return 10 ** (np.log10(x) * slope + intercept)

def determine_simulation_type():
    datafile = sys.argv[2]
    for sub in ['2S', '1S', '2B', '1B']:
        if sub in datafile:
            print("The datafile contains one of the specified substrings.")
            if sub == '1S':
                halo = 1
                bh_mass = 10.8
                sim = 's1-10.8msun-'
            elif sub == '2S':
                halo = 2
                bh_mass = 10.8
                sim = 's2-10.8msun-'
            elif sub == '1B':
                halo = 1
                bh_mass = 270
                sim = 's1-270msun-'
            elif sub == '2B':
                halo = 2
                bh_mass = 270
                sim = 's2-270msun-'
            return bh_mass, halo, sim

def setup_fonts(fontsize):
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = 1.5

def setup_plot(fontsize, linewidth):
    """ Set up the plot with font and line settings """
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth

def get_dx_values(sim):
    dx_mapping = {
        "s1-270msun-": [1.229791e-02, 3.074475e-03, 1.537645e-03, 7.692833e-04],
        "s1-10.8msun-": [2.459867e-02, 1.229940e-02, 3.074829e-03, 7.692833e-04],
        "s2-270msun-": [8.298311e-03, 2.074568e-03, 1.537645e-03, 7.687095e-04, 3.8435475e-04, 1.296596e-04],
        "s2-10.8msun-": [0.00832, 0.00416, 0.002079909, 0.001297054, 0.0005188221, 0.0001297054],
        "s2-10.8msun-2-": [1.3e-03, 5.2e-04, 1.3e-04]
    }
    return dx_mapping.get(sim, [])

def process_bhl_data(bhl_object_list, dx, time_cutoff, atol, eddington, window_size, c, l, axs, i_start=5, N=1000):
    lst_accrate_bhl, lst_res_hl_bhl, lst_res_bondi_bhl = [], [], []
    lst_accrate_mf, lst_res_hl_mf, lst_res_bondi_mf = [], [], []

    # Determine N evenly spaced points for the interpolated data
    baseline_age_raw = bhl_object_list[0].ages[i_start:]
    n_data_max = baseline_age_raw.shape[0]
    N = int(n_data_max/100)

    # Determine the indices of the BHL and MF objects
    i_mf = np.arange(int(len(bhl_object_list) / 2))
    i_bhl = np.arange(int(len(bhl_object_list) / 2), len(bhl_object_list))
    alpha = 0.5
    min_hl = max_hl = min_bondi = max_bondi = 1

    for i, BHL in enumerate(bhl_object_list):
        # [Add your existing loop logic here, including plotting]
        # Append the processed data to the respective lists

        # convert ages from yrs to Myrs
        BHL.ages = np.array(BHL.ages) / 1e6

        # find index of age that matches end age of time limit
        i_age = first_index(BHL.ages[i_start:], time_cutoff, rtol=1e-8, atol=atol) # 7e-2 for 1B.m, 1e-4 otherwise

        # calculate mass and hl_radius moving averages
        accrate = moving_average(BHL.accrates[i_start:i_age], window_size)
        eddrate = eddington_rate(moving_average(BHL.mass[i_start:i_age], window_size))
        if eddington:
            accrate = accrate/eddrate
        hl_radius = moving_average(BHL.hl_radius[i_start:i_age], window_size)
        bondi_radius = moving_average(BHL.bondi_radius[i_start:i_age], window_size)

        # calculate how many cells it's resolving the hl radius by
        if i in i_mf:
            dx_res_hl = hl_radius/dx[i]
            dx_res_bondi = bondi_radius / dx[i]
        else: 
            dx_res_hl = hl_radius / dx[i - int(len(bhl_object_list) / 2)]
            dx_res_bondi = bondi_radius / dx[i - int(len(bhl_object_list) / 2)]

        # interpolate this array over N evenly spaced points
        interp_accrate = interpolate_data(accrate, N=N)
        interp_eddrate = interpolate_data(eddrate, N=N)
        interp_hl = interpolate_data(dx_res_hl, N=N)
        interp_bondi = interpolate_data(dx_res_bondi, N=N)
        accrate = interp_accrate
        eddrate = interp_eddrate
        dx_res_hl = interp_hl
        dx_res_bondi = interp_bondi

        if i in i_mf:

            # add the ages and radii to lst (for Pearson coeff calculation)
            lst_accrate_mf.append(accrate)
            lst_res_hl_mf.append(dx_res_hl)
            lst_res_bondi_mf.append(dx_res_bondi)

            # 1) HL radius resolution in cell widths
            axs[0, 0].scatter(dx_res_hl, accrate, color=c[i], linestyle='solid', label=l[i], alpha=alpha)

            # 2) Bondi radius resolution in cell widths
            axs[0, 1].scatter(dx_res_bondi, accrate, color=c[i], linestyle='solid', label=l[i], alpha=alpha)

        else:

            # add the ages and radii to lst (for Pearson coeff calculation)
            lst_accrate_bhl.append(accrate)
            lst_res_hl_bhl.append(dx_res_hl)
            lst_res_bondi_bhl.append(dx_res_bondi)

            # 1) HL radius resolution in cell widths
            axs[1, 0].scatter(dx_res_hl, accrate, color=c[i - int(len(bhl_object_list) / 2)],
                              linestyle='solid', label=l[i], alpha=alpha)

            # 2) Bondi radius resolution in cell widths
            axs[1, 1].scatter(dx_res_bondi, accrate, color=c[i - int(len(bhl_object_list) / 2)],
                              linestyle='solid', label=l[i], alpha=alpha)

        # update y-axis limits
        min_hl = min(dx_res_hl.min(), min_hl)
        max_hl = max(dx_res_hl.max(), max_hl)
        min_bondi = min(dx_res_bondi.min(), min_bondi)
        max_bondi = max(dx_res_bondi.max(), max_bondi)

    return lst_accrate_bhl, lst_res_hl_bhl, lst_res_bondi_bhl, lst_accrate_mf, lst_res_hl_mf, lst_res_bondi_mf

def perform_analysis_and_plot_trendlines(lst_accrate_bhl, lst_res_hl_bhl, lst_res_bondi_bhl, lst_accrate_mf, lst_res_hl_mf, lst_res_bondi_mf, axs, fontsize):
    lst_res_hl_bhl = np.concatenate(lst_res_hl_bhl)
    lst_res_hl_mf = np.concatenate(lst_res_hl_mf)
    lst_accrate_bhl = np.concatenate(lst_accrate_bhl)
    lst_accrate_mf = np.concatenate(lst_accrate_mf)
    lst_res_bondi_mf = np.concatenate(lst_res_bondi_mf)
    lst_res_bondi_bhl = np.concatenate(lst_res_bondi_bhl)
    
    # perform Pearson correlation coefficient analysis
    # 1) HL radius resolution in cell widths
    corr_hl_bhl, pv_hl_bhl = stats.pearsonr(np.log10(lst_res_hl_bhl), np.log10(lst_accrate_bhl))
    corr_hl_mf, pv_hl_mf = stats.pearsonr(np.log(lst_res_hl_mf), np.log10(lst_accrate_mf))
    print("HL radius resolution in cell widths:")
    print("BHL: corr = {:.2f}, p = {:.3e}".format(corr_hl_bhl, pv_hl_bhl))
    print("MF: corr = {:.2f}, p = {:.3e}".format(corr_hl_mf, pv_hl_mf))

    # 2) Bondi radius resolution in cell widths
    corr_bondi_bhl, pv_bondi_bhl = stats.pearsonr(np.log10(lst_res_bondi_bhl), np.log10(lst_accrate_bhl))
    corr_bondi_mf, pv_bondi_mf = stats.pearsonr(np.log10(lst_res_bondi_mf),np.log10(lst_accrate_mf))
    print("Bondi radius resolution in cell widths:")
    print("BHL: corr = {:.2f}, p = {:.3e}".format(corr_bondi_bhl, pv_bondi_bhl))
    print("MF: corr = {:.2f}, p = {:.3e}".format(corr_bondi_mf, pv_bondi_mf))

    # perform linear regression on aggregate data per accretion scheme and radius resolution
    # 1) BHL
    slope_hl, intercept_hl, r_value_hl, p_value_hl, std_err_hl = stats.linregress(np.log10(lst_res_hl_bhl), np.log10(lst_accrate_bhl))
    slope_bondi, intercept_bondi, r_value_bondi, p_value_bondi, std_err_bondi = stats.linregress(np.log10(lst_res_bondi_bhl), np.log10(lst_accrate_bhl))
    line_str_pearson_hl_bhl = "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}$".format(r_value_hl**2, corr_hl_bhl)
    line_str_pearson_bondi_bhl = "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}*$".format(r_value_bondi**2, corr_bondi_bhl) if sim == 's2-270msun-' and xlim > 0.1 else "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}$".format(r_value_bondi**2, corr_bondi_bhl)
    p_hl_bhl = [slope_hl, intercept_hl]
    p_bondi_bhl = [slope_bondi, intercept_bondi]

    # 2) MF
    slope_hl_2, intercept_hl_2, r_value_hl_2, p_value_hl_2, std_err_hl_2 = stats.linregress(np.log10(lst_res_hl_mf), np.log10(lst_accrate_mf))
    slope_bondi_2, intercept_bondi_2, r_value_bondi_2, p_value_bondi_2, std_err_bondi_2 = stats.linregress(np.log10(lst_res_bondi_mf), np.log10(lst_accrate_mf))
    line_str_pearson_hl_mf = "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}$*".format(r_value_hl_2 ** 2, corr_hl_mf) if sim == 's1-10msun-' else "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}$".format(r_value_hl_2 ** 2, corr_hl_mf)
    line_str_pearson_bondi_mf = "$R^2$ = ${:.2f}$\n $P$ = ${:.2f}$".format(r_value_bondi_2 ** 2, corr_bondi_mf)
    p_hl_mf = [slope_hl_2, intercept_hl_2]
    p_bondi_mf = [slope_bondi_2, intercept_bondi_2]

    # Plot trendlines
    linear_fit_color = "yellow"
    for i in range(num_subplots):
        for j in range(num_subplots):
            axs[i, j].set_yscale('log')
            axs[i, j].set_xscale('log')
    axs[0, 0].plot(lst_res_hl_mf, log_line_fit(lst_res_hl_mf, *p_hl_mf), linear_fit_color, linestyle='solid')
    axs[0, 1].plot(lst_res_bondi_mf, log_line_fit(lst_res_bondi_mf, *p_bondi_mf), linear_fit_color, linestyle='solid')
    axs[1, 0].plot(lst_res_hl_bhl, log_line_fit(lst_res_hl_bhl, *p_hl_bhl), linear_fit_color, linestyle='solid')
    axs[1, 1].plot(lst_res_bondi_bhl, log_line_fit(lst_res_bondi_bhl, *p_bondi_bhl), linear_fit_color, linestyle='solid')

    return p_hl_bhl, p_bondi_bhl, p_hl_mf, p_bondi_mf, line_str_pearson_hl_bhl, line_str_pearson_bondi_bhl, line_str_pearson_hl_mf, line_str_pearson_bondi_mf


def format_plots(axs, num_subplots, fontsize, linewidth, alpha, eddington, line_str_pearson_bondi_mf, line_str_pearson_hl_bhl, line_str_pearson_bondi_bhl):
    # Set up the plot with font and line settings
    for i in range(num_subplots):
        for j in range(num_subplots):
            # axs[i, j].set_yscale('log')
            # axs[i, j].set_xscale('log')
            #axs[i, j].set_yticks(yticks)
            #axs[i, j].yaxis.set_minor_locator(FixedLocator(yticks_minor))
            axs[i, j].tick_params(axis="both", which='both', length=2, direction="in", labelsize=fontsize)
            axs[i, j].axvline(x=1, color="grey", linestyle='--', lw=linewidth, alpha=alpha)
            axs[i, j].set_xlim(8e-4, 1.5e2)
            axs[i, j].set_ylim(3e-5,3e3)
            

    # Specific axis labels and settings
    axs[1, 0].set_xlabel(r"$R_{\rm HL}$ Resolution (cell widths)")
    axs[1, 1].set_xlabel(r"$R_{\rm Bondi}$ Resolution (cell widths)")
    ylabel = r"$\rm \dot{M} / \dot{M}_{Edd}$" if eddington else r"$\rm \dot{M} (M_\odot/yr) $"
    axs[1, 0].set_ylabel(ylabel)
    axs[0, 0].set_ylabel(ylabel)
    
    # Additional text and legend
    xtext = 0.705 if sim == 's2-270msun-' else 0.69
    ytext = 0.085 if sim == 's2-270msun-' else 0.085
    xtext2 = ytext if sim == 's2-270msun-' else 0.05
    #text1 = f"{line_str_pearson_hl_mf}\n{line2}"
    axs[0, 1].text(1.02, 0.5, 'Mass-Flux', transform=axs[0, 1].transAxes, ha='left', va='center')
    axs[1, 1].text(1.02, 0.5, 'BHL', transform=axs[1, 1].transAxes, ha='left', va='center')
    axs[0, 0].text(xtext, ytext, line_str_pearson_hl_mf, fontsize=fontsize, bbox=dict(facecolor='white', edgecolor='lightgrey'), transform=axs[0, 0].transAxes)
    axs[0, 1].text(xtext2, ytext, line_str_pearson_bondi_mf, fontsize=fontsize, bbox=dict(facecolor='white', edgecolor='lightgrey'), transform=axs[0, 1].transAxes)
    axs[1, 0].text(xtext, ytext, line_str_pearson_hl_bhl, fontsize=fontsize, bbox=dict(facecolor='white', edgecolor='lightgrey'), transform=axs[1, 0].transAxes)
    axs[1, 1].text(xtext2, ytext, line_str_pearson_bondi_bhl, fontsize=fontsize, bbox=dict(facecolor='white', edgecolor='lightgrey'), transform=axs[1, 1].transAxes)
    axs[0, 1].legend(fontsize=fontsize - 2, ncol=2, loc="upper left")
    axs[1, 1].legend(fontsize=fontsize - 2, ncol=2, loc="upper left")

    # Remove y-tick labels from top row
    axs[0, 1].set_yticklabels([])
    axs[1, 1].set_yticklabels([])

    if sim == 's1-10.8msun-':
        # limit y-axis range in top row
        axs[0, 0].set_ylim(8e-2, 3e3)
        axs[0, 1].set_ylim(8e-2, 3e3)

    if sim == 's1-270msun-':
        # limit x-axis range in all plots
        axs[0, 0].set_xlim(6e-2, 2e3)

        # limit y-axis range in top row
        axs[0, 0].set_ylim(8e-1, 5e3)
        axs[0, 1].set_ylim(8e-1, 5e3)

        # limit y-axis range in bottom row
        axs[1, 0].set_ylim(8e-1, 5e3)
        axs[1, 1].set_ylim(8e-1, 5e3)

    if sim == 's2-10.8msun-':
        # limit y-axis range in top row
        axs[0, 0].set_ylim(8e-2, 3e3)
        axs[0, 1].set_ylim(8e-2, 3e3)

    if sim == 's2-270msun-':
        # limit x-axis range in top row
        axs[0, 0].set_xlim(7e-3, 8e2)

        # limit y-axis range in top row
        axs[0, 0].set_ylim(8e-2, 3e3)
        axs[0, 1].set_ylim(8e-2, 3e3)

        # limit y-axis range in bottom row
        axs[1, 0].set_ylim(2, 5e3)
        axs[1, 1].set_ylim(2, 5e3)



if __name__ == "__main__":

    ##########################################################################################################
    #                               Plot Accretion Rate vs Radius Resolution vs 
    #
    # to run: python -i plot_radius_relationship.py [csv1] [csv2] [csv3] [output_plotname e.g MF-BHL]
    # for 2x2 update: list MF runs first, then BHL runs. Name: MF+BHL (in order of low res to high res)
    ##########################################################################################################

    # Define the plotting parameters
    xlim = 0.4
    atol = 5e-4
    eddington = True
    accretion = sys.argv[-1]
    fontsize = 12
    window_size = 1
    linewidth = 1.5
    alpha = 0.5

    bh_mass, halo, sim = determine_simulation_type()
    print("Simulation: ", sim)
    print("Black hole mass: ", bh_mass)
    print("Halo: ", halo)

    # Define the colors and labels
    colors = ['blueviolet', 'turquoise', 'limegreen', 'darkgreen'] if bh_mass == 270 else ['goldenrod', 'chocolate', 'orangered', 'violet', 'deeppink', 'blueviolet',] # line colours s2S
    colors = [ 'violet', 'deeppink', 'blueviolet','turquoise'] if bh_mass == 10.8 and halo == 1 else colors

    # Set up the plot
    setup_fonts(fontsize)
    num_subplots = 2
    fig, axs = plt.subplots(num_subplots, num_subplots, sharex=True, sharey=False)
    labels = tidy_data_labels(bhl_object_labels)
    dx = get_dx_values(sim)

    # Process BHL data
    lst_accrate_bhl, lst_res_hl_bhl, lst_res_bondi_bhl, lst_accrate_mf, lst_res_hl_mf, lst_res_bondi_mf = process_bhl_data(
        bhl_object_list, dx, xlim, atol, eddington, window_size, colors, labels, axs
    )

    # Perform analysis and plot trendlines
    p_hl_bhl, p_bondi_bhl, p_hl_mf, p_bondi_mf, line_str_pearson_hl_bhl, line_str_pearson_bondi_bhl, line_str_pearson_hl_mf, line_str_pearson_bondi_mf = perform_analysis_and_plot_trendlines(
        lst_accrate_bhl, lst_res_hl_bhl, lst_res_bondi_bhl, lst_accrate_mf, lst_res_hl_mf, lst_res_bondi_mf, axs, fontsize
    )

    # Format and plot trendlines
    min_hl = max_hl = min_bondi = max_bondi = 1
    accrate_total_mf, accrate_total_bhl = np.concatenate(lst_accrate_mf), np.concatenate(lst_accrate_bhl)
    yticklabels = [r"$\rm 10^{-3}$", r"$\rm 10^{-2}$", r"$\rm 10^{-1}$", r"$\rm 10^{0}$", r"$\rm 10^{1}$", r"$\rm 10^{2}$", r"$\rm 10^{3}$"]
    format_plots(axs, num_subplots, fontsize, linewidth, alpha, eddington, line_str_pearson_bondi_mf, line_str_pearson_hl_bhl, line_str_pearson_bondi_bhl)

    # save plot as pdf
    #fig = plt.gcf()
    fig.subplots_adjust(wspace=0.005, hspace=0.005)
    fig.set_size_inches(7, 5)
    fig.suptitle(f"Halo {halo} {bh_mass} $M_\odot$ Black Hole, t = {xlim} Myr", fontsize=fontsize + 5, y=0.95)
    fig.tight_layout()
    plot_name = 'plots/accrate-dx_res-' + str(sim) + str(accretion) + "-" + str(xlim) + 'Myr_2.pdf'
    fig.savefig(plot_name)
    print("created ", plot_name)