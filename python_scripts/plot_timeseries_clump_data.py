"""
Plots a 2 panels, one with a scatter of clump distance vs. age with clump mass bins, and the other with total clump mass vs. age.
There are extra yaxes for BH mass on the first panel and number of clumps on the second panel.
"""

import matplotlib.pyplot as plt
import numpy as np
import yt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator, LogFormatterExponent
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.transforms import IdentityTransform
import matplotlib.font_manager as fm
from helper_functions import extract_simulation_name
from csv_data_plotting.clump_analysis import setup_plot_env, bin_clump_mass, determine_num_bins, extract_simulation_name

def filter_clump_data(df_clump_distances, simulation_name, min_density_value, max_clump_distance, max_mass, max_dist):
    """
    Filter clump data based on user-specified parameters.
    """
    filtered_clump_masses = df_clump_distances[
        (df_clump_distances['simulation'] == simulation_name) &
        (df_clump_distances['clump_distance'] <= max_clump_distance) &
        (df_clump_distances['min_density'] == float(min_density_value))
    ]

    to_remove = filtered_clump_masses[
        (filtered_clump_masses['clump_mass'] > max_mass) & 
        (filtered_clump_masses['clump_distance'] < max_dist)
    ]
    print("Clumps to be removed by the filtering condition:")
    print(to_remove)

    filtered_clump_masses = filtered_clump_masses[
        ~((filtered_clump_masses['clump_mass'] > max_mass) & 
          (filtered_clump_masses['clump_distance'] < max_dist))
    ]

    return filtered_clump_masses

def calculate_total_clump_mass(filtered_clump_masses):
    """
    Calculate the total clump mass for each age.
    """
    total_clump_mass = filtered_clump_masses.groupby('age_myr')['clump_mass'].sum().reset_index()
    total_clump_mass.columns = ['age_myr', 'total_clump_mass']
    return total_clump_mass

def calculate_max_clump_mass(filtered_clump_masses):
    """
    Calculate the maximum clump mass for each age.
    """
    max_clump_mass = filtered_clump_masses.groupby('age_myr')['clump_mass'].max().reset_index()
    max_clump_mass.columns = ['age_myr', 'max_clump_mass']
    return max_clump_mass

def calculate_number_of_clumps(filtered_clump_masses):
    """
    Calculate the number of clumps for each age.
    """
    number_of_clumps = filtered_clump_masses.groupby('age_myr').size().reset_index(name='no_clumps')
    return number_of_clumps

def prepare_filt_data(df_clump_data, simulation_name, total_clump_mass):
    """
    Prepare the filtered data by merging total clump mass and sorting by age.
    """
    filt_data = pd.merge(df_clump_data[df_clump_data['simulation'] == simulation_name], total_clump_mass, on='age_myr', how='left')
    filt_data = filt_data.sort_values(by='age_myr')
    return filt_data

def plot_clump_analysis(filtered_clump_masses, total_clump_mass, max_clump_mass, filt_data, simulation_name, max_clump_distance, alpha,
                        ylim1=None, ylim2=None, ylim3=None, ylim4=None, ylim5=None, xlim=None):
    """
    Plot the clump analysis including Clump Distance vs. Age and Clump Mass.
    """
    fig, (ax1, ax4) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Determine number of bins for clump mass
    number_of_clumps = calculate_number_of_clumps(filtered_clump_masses)
    overall_max_clump_mass = max_clump_mass.max()['max_clump_mass']
    print("Max clump mass:", overall_max_clump_mass)
    num_bins = determine_num_bins(overall_max_clump_mass)
    print("Number of bins:", num_bins)
    palette = ["#FAC64C", "#8CC65E",  "#3fcab5", "#3A96D0","#ad6cd0", "#c34177", "#FF5A1F", "#BB2020", 
               "#5B203E","#3D0000","#C70303", "#3D0000"][:num_bins]
    bin_order = ['0-10', '11-100', '101-200', '201-300', '301-400', '401-500', '501-1000', '1001-2000', '2000+'][:num_bins]

    # Plot Clump Distance vs Age with Clump Mass bins
    sns.scatterplot(data=filtered_clump_masses, x='age_myr', y='clump_distance', hue='clump_mass_bin', hue_order=bin_order, palette=palette, alpha=alpha, 
                ax=ax1, s=100, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Age (Myr)')
    ax1.set_ylabel('Clump Distance (pc)')
    ax1.set_title(f'{simulation_name} (Clump Distance $\leq$ {max_clump_distance} pc)')
    ax1.grid(True)
    ax1.legend(title=r'Clump Mass ($\rm M_\odot$)', bbox_to_anchor=(1.15, 1), loc='upper left')

    # Secondary y-axis for BH Mass
    ax2 = ax1.twinx()
    ax2.set_yscale('log')
    sns.lineplot(data=filt_data, x='age_myr', y='bh_mass', color='grey', ax=ax2)
    ax2.set_ylabel('BH Mass ($M_\\odot$)')

    # Set ylim for ax1 and ax2 if provided
    if ylim1:
        ax1.set_ylim(*ylim1)
    if ylim2:
        ax2.set_ylim(*ylim2)

    # # Top x-axis for total number of clumps
    # ax3 = ax1.twiny()
    # ax3.set_xlim(ax1.get_xlim())
    # xticks = number_of_clumps['age_myr'].values
    # ax3.set_xticks(xticks)
    # ax3.set_xticklabels(number_of_clumps['no_clumps'])
    # ax3.set_xlabel('Number of Clumps')

    # Plot Total Clump Mass vs Age
    sns.lineplot(data=total_clump_mass, x='age_myr', y='total_clump_mass', ax=ax4, marker='o', label='Total Clump Mass')
    sns.lineplot(data=max_clump_mass, x='age_myr', y='max_clump_mass', ax=ax4, marker='x', color='red', label='Max Clump Mass')
    ax4.set_xlabel('Age (Myr)')
    ax4.set_ylabel(r'Clump Mass ($\rm M_\odot$)')
    ax4.legend(loc='upper left')
    ax4.grid(True)

    # Set ylim for ax4 if provided
    if ylim4:
        ax4.set_ylim(*ylim4)

    # Plot Number of Clumps
    ax5 = ax4.twinx()
    sns.lineplot(data=number_of_clumps, x='age_myr', y='no_clumps', ax=ax5, marker='D', color='seagreen', label='Number of Clumps')
    ax5.set_ylabel('Number of Clumps')
    ax5.legend(loc='upper right')

    # Set ylim for ax5 if provided
    if ylim5:
        ax5.set_ylim(*ylim5)

    # Set xlim if provided
    if xlim:
        ax1.set_xlim(*xlim)
        ax4.set_xlim(*xlim)

    fig.tight_layout()

    # Save the figure
    figname = f'plots/clump_scatter_2panels_{simulation_name}_max_dist_{max_clump_distance}pc.png'
    plt.savefig(figname, dpi=300)
    print(f"Saved plot to {figname}")
    plt.show()

if __name__ == "__main__":
    # User-set parameters
    fontsize = 12
    linewidth = 2
    window_size = 100
    setup_plot_env(fontsize, linewidth)

    filter_method = "complete_name"  # Choose "prefix" or "complete_name"
    simulation_name = "1B.resim.th.b04-eta-0.1"
    min_density_value = 597.5500448162534
    max_clump_distance = 2
    max_mass = 500
    max_dist = 2
    alpha = 0.9
    title = "Accretion Event 2: 2.50 Myr - 2.80 Myr"

    # Optional limits for axes
    ylim1 = None        # Clump distance
    # ylim2 = (3e3, 6e3)  # BH mass
    # ylim3 = (0, 5)      # Total number of clumps
    # ylim4 = (0, 100)    # Total Clump Mass
    # ylim5 = (0, 5)      # Total number of clumps
    # xlim = (2.2, 3.2)   # x-axis limits for ax1 and ax4
    ylim2 = None    # BH mass
    ylim3 = None    # Total number of clumps
    ylim4 = None    # Total Clump Mass
    ylim5 = None    # Total number of clumps
    xlim = None     # x-axis limits for ax1 and ax4

    # Load the data
    df_clump_distances = pd.read_csv('csv_data_plotting/data_misc/clump_distances.csv')
    df_clump_data = pd.read_csv('csv_data_plotting/data_misc/clump_data.csv')

    # Apply custom binning
    df_clump_distances['clump_mass_bin'] = df_clump_distances['clump_mass'].apply(bin_clump_mass)
    df_clump_distances = df_clump_distances.round({'age_myr': 2, 'clump_distance': 2})

    # Perform filtering and data preparation
    filtered_clump_masses = filter_clump_data(df_clump_distances, simulation_name, min_density_value, max_clump_distance, max_mass, max_dist)

    # Further filter the data for the bottom plots
    filtered_clump_masses_bottom = filtered_clump_masses[
        (filtered_clump_masses['clump_mass'] <= max_mass) &
        (filtered_clump_masses['clump_distance'] <= max_dist)
    ]

    # Apply xlim filtering if xlim is defined
    if xlim is not None:
        filtered_clump_masses_bottom = filtered_clump_masses_bottom[
            (filtered_clump_masses_bottom['age_myr'] >= xlim[0]) &
            (filtered_clump_masses_bottom['age_myr'] <= xlim[1])
        ]

    # Remove duplicate rows 
    filtered_clump_masses_bottom = filtered_clump_masses_bottom.drop_duplicates()

    # Recalculate total clump mass and max clump mass for the filtered data
    total_clump_mass = calculate_total_clump_mass(filtered_clump_masses_bottom)
    max_clump_mass = calculate_max_clump_mass(filtered_clump_masses_bottom)

    # Prepare data for plotting BH mass and number of clumps
    filt_data = prepare_filt_data(df_clump_data, simulation_name, total_clump_mass)

    # Plot the data
    plot_clump_analysis(filtered_clump_masses_bottom, total_clump_mass, max_clump_mass, filt_data, simulation_name, max_clump_distance, alpha,
                        ylim1=ylim1, ylim2=ylim2, ylim3=ylim3, ylim4=ylim4, ylim5=ylim5, xlim=xlim)



