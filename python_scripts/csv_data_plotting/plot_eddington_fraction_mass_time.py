import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from plot_variables import *
from helper_functions import setup_plot_env
from matplotlib.ticker import MaxNLocator
from colors import colors

def extract_simulation_name_from_csv(fp):
    """
    Extract the simulation name from a .csv data file path.
    e.g. 'data_files/data-1B.resim.th.b01.csv' -> '1B.resim.th.b01'
    e.g. 'data_files/data-1B.resim.th.b04-eta-0.1-fb-r-7dx.csv' -> '1B.resim.th.b04-eta-0.1-fb-r-7dx'

    Parameters:
    fp (str): The file path.

    Returns:
    str: The extracted simulation name.
    """
    # Extract the base filename without directory and extension
    base_name = os.path.basename(fp)  
    name_without_extension = os.path.splitext(base_name)[0] 
    
    # Extract everything after the first dash
    match = re.search(r'data-(.+)', name_without_extension)
    
    if match:
        return match.group(1)  # Return everything after the dash

    print("No match found")
    return None

def compute_fractions(f_edd, mass_gained, time_intervals):
    # Define the Eddington regimes
    regimes = {
        r"$f_{\rm Edd} < 0.01$": (0, 0.01),
        r"$0.01 \leq f_{\rm Edd} \leq 1$": (0.01, 1),
        r"$f_{\rm Edd} > 1$": (1, np.inf)
    }
    
    mass_fraction = {}
    time_fraction = {}
    
    total_mass = np.sum(mass_gained)
    total_time = np.sum(time_intervals)
    
    for regime_name, (low, high) in regimes.items():
        # Mask the data based on the regime
        mask = (f_edd >= low) & (f_edd < high)
        
        # Calculate the fraction of mass and time in this regime
        mass_fraction[regime_name] = np.sum(mass_gained[mask]) / total_mass
        time_fraction[regime_name] = np.sum(time_intervals[mask]) / total_time
    
    return mass_fraction, time_fraction

def load_data(data_file, n=10, time_cutoff=None):
    # Load the CSV file into a DataFrame
    df = pd.read_csv(data_file)

        # Apply time cutoff if specified
    if time_cutoff is not None:
        df = df[df['age'] <= time_cutoff*1e6]

    # Extract the columns you're interested in
    age = df['age'].values / 1e6  # Convert age to Myr
    bh_mass = df['BH mass'].values
    accrate = df['accrate'].values

    # Apply adaptive smoothing
    age_av     = adaptive_moving_average(age, window_size=n)
    bh_mass_av = adaptive_moving_average(bh_mass, window_size=n)
    accrate_av = adaptive_moving_average(accrate, window_size=n)
    eddrate_av = eddington_rate(bh_mass_av)
    f_edd = accrate_av / eddrate_av

    # Calculate the mass gained between consecutive time intervals
    mass_gained = np.diff(bh_mass_av)
    time_intervals = np.diff(age_av)

    print(f"Starting mass: {bh_mass_av[0]:.2e} Msun")
    print(f"Ending mass: {bh_mass_av[-1]:.2e} Msun")
    print(f"Total mass gained: {np.sum(mass_gained):.2e} Msun")
    print(f"Starting time: {age_av[1]:.2f} Myr")
    print(f"Ending time: {age_av[-1]:.2f} Myr")
    print(f"Total time: {age_av[-1] - age_av[1]:.2f} Myr")

    return f_edd[:-1], mass_gained, time_intervals

def plot_for_data_files(prefix, color='coral', alpha=0.8, single_simulation=False,time_cutoff=None):
    # Set up the plot environment
    setup_plot_env(fontsize=14, linewidth=1.5)

    # Ensure output directory exists
    output_dir = "plots/barplots_edd_fraction"
    os.makedirs(output_dir, exist_ok=True)

    if single_simulation:
        # Use the exact prefix for a single file
        file_list = [f"data_files/data-{prefix}.csv"]
    else:
        # Use wildcard to match multiple files with the given prefix
        file_pattern = f"data_files/data-{prefix}*.csv"
        file_list = glob.glob(file_pattern)

    for data_file in file_list:
        # Extract simulation name from file
        sim_name = extract_simulation_name_from_csv(data_file)
        
        # Load data
        f_edd, mass_gained, time_intervals = load_data(data_file, time_cutoff=time_cutoff)

        # Compute fractions
        mass_fraction, time_fraction = compute_fractions(f_edd, mass_gained, time_intervals)
        
        # Plotting
        labels = [r"$f_{\rm Edd} < 0.01$", r"$0.01 \leq f_{\rm Edd} \leq 1$", r"$f_{\rm Edd} > 1$"]
        mass_values = [mass_fraction[label] for label in labels]
        time_values = [time_fraction[label] for label in labels]
        
        x = np.arange(len(labels))  # the label locations
        width = 0.4  # the width of the bars
        
        _, ax = plt.subplots(figsize=(5, 5))
        
        # Overlap the bars by setting them to the same positions
        ax.bar(x, mass_values, width, label='Mass Gained', color=colors[color], alpha=alpha)
        ax.bar(x, time_values, width, label='Time', color=colors['lilac'], hatch='//', alpha=0.3)
        
        # Add title and labels
        ax.set_title(f'{sim_name} ({np.sum(time_intervals):.2f} Myr)')
        ax.set_ylabel('Fraction')
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylim(0, 1)
        ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax.legend()
        
        # Save the plot
        filename = f"{output_dir}/barplot_edd_fraction_{sim_name}_{time_cutoff}Myr.pdf" if time_cutoff is not None else f"{output_dir}/barplot_edd_fraction_{sim_name}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        plt.close()  # Close the figure to save memory
        print(f"Saved plot to {filename}")

def plot_combined_for_data_files(simulation_prefixes, color_lst, alpha=0.8, time_cutoff=None):
    """
    Generate combined plots for multiple simulations in a single row of panels.

    Parameters:
    simulation_prefixes (list of str): List of simulation prefixes to match the data files.
    color_lst (list of str): List of colors to use for each simulation.
    alpha (float): Transparency level for the bars.
    time_cutoff (float, optional): Maximum time cutoff for the simulations (in Myr).
    """
    # Set up the plot environment
    setup_plot_env(fontsize=14, linewidth=1.5)

    # Ensure output directory exists
    output_dir = "plots/barplots_edd_fraction_combined"
    os.makedirs(output_dir, exist_ok=True)

    # Number of subplots (one for each simulation)
    n_plots = len(simulation_prefixes)

    # Create figure and subplots
    fig, axes = plt.subplots(1, n_plots, figsize=(4 * n_plots, 4), sharey=True)
    # Reduce horizontal space between subplots
    plt.subplots_adjust(wspace=0.05) 


    # Iterate over the simulations
    for i, prefix in enumerate(simulation_prefixes):
        # Find the data file
        file_pattern = f"data_files/data-{prefix}.csv"
        file_list = glob.glob(file_pattern)
        print(f"Found {len(file_list)} files for prefix: {prefix}")

        if not file_list:
            print(f"No files found for prefix: {prefix}")
            continue

        # Use the first matching file
        data_file = file_list[0]
        sim_name = extract_simulation_name_from_csv(data_file)

        # Load data
        f_edd, mass_gained, time_intervals = load_data(data_file, time_cutoff=time_cutoff)

        # Compute fractions
        mass_fraction, time_fraction = compute_fractions(f_edd, mass_gained, time_intervals)

        # Calculate the final BH mass and format it
        bh_mass_final = np.sum(mass_gained)
        bh_mass_text = r"Mass gained: {bh_mass_final:.2e} $M_\odot$"

        # Plotting on the corresponding axis
        ax = axes[i] if n_plots > 1 else axes
        labels = [r"$f_{\rm Edd} < 0.01$", r"$0.01 \leq f_{\rm Edd} \leq 1$", r"$f_{\rm Edd} > 1$"]
        mass_values = [mass_fraction[label] for label in labels]
        time_values = [time_fraction[label] for label in labels]

        x = np.arange(len(labels))  # the label locations
        width = 0.5  # the width of the bars

        # Overlap the bars by setting them to the same positions
        ax.bar(x, mass_values, width, label='Mass Gained', color=colors[color_lst[i]], alpha=alpha)
        ax.bar(x, time_values, width, label='Time', color=colors['lilac'], hatch='//', alpha=0.3)

        # Add BH mass text to the plot
        ax.text(0.5, 0.95, bh_mass_text, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', horizontalalignment='center')
        
        # Add legend on first and last subplot
        if i == 0:
            ax.legend()
        if i == n_plots - 1:
            ax.legend()

        # Add title and labels
        ax.set_title(f'{sim_name}') # ({np.sum(time_intervals):.2f} Myr)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=12)
        ax.set_ylim(0, 1)

        # Only show the y-axis label on the first subplot
        if i == 0:
            ax.set_ylabel('Fraction')

    # Save the combined plot
    filename = f"{output_dir}/combined_barplot_edd_fraction_{'_'.join(simulation_prefixes)}_{time_cutoff}Myr.pdf" if time_cutoff is not None else f"{output_dir}/combined_barplot_edd_fraction_{'_'.join(simulation_prefixes)}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()  # Close the figure to save memory
    print(f"Saved combined plot to {filename}")

# Example usage: Combine multiple simulations into one row of panels
simulation_prefixes = ["1B.resim.th.b01-eps-0.001-2", 
                       "1B.resim.th.b04-eps-0.001-2", 
                       "1B.resim.th.b04-eta-0.1-fb-r-7dx", 
                       "1B.resim.th.b04-eta-0.1"]
color_lst = ['pale_green', 'coral', 'coral', 'coral']

plot_combined_for_data_files(simulation_prefixes, color_lst, alpha=0.7, time_cutoff=2.9)


# Example usage: Generate plots for all files with the specified prefix
# Event 1 1B.resim.th.b04*: coral 2.8 Myr cutoff
# Event 2 1B.resim.th.b01-3*: cornflower_blue 31.9 Myr cutoff
# Event 3 1B.resim.th.b01-e*: pale_green
# No-feedback: light_pink 
# 1B.resim.th.b01-eps-0.001-2, 1B.resim.th.b04-eps-0.001-2 , 1B.resim.th.b04-eta-0.1-fb-r-7dx, 1B.resim.th.b04-eta-0.1
#plot_for_data_files("1B.resim.th.b04-eta-0.1", color='coral', alpha=0.7, single_simulation=True, time_cutoff=None)
