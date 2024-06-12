import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Load the data from the CSV file
df = pd.read_csv('clump_data.csv', comment='#')

# Apply filters (optional)
df['future_bound'] = df['future_bound'].map({' True': True, ' False': False})
df = df.round({'age_myr': 2, 'max_clump_mass': 1, 'total_clump_mass': 1, 'bh_mass':1})
filt = df[df['future_bound'] == True]  # Keep only future bound clumps
filt = filt[filt['min_density'] == '597550.0448162535']

# Sort the DataFrame by 'age_myr' to ensure the plots are in chronological order
filt = filt.sort_values('age_myr')

# Unique simulations and BH Age values for grouping
simulations = filt['simulation'].unique()
ages = filt['age_myr'].unique()
properties = ["max_clump_mass", "total_clump_mass", "bh_mass"]

# Generate a color for each simulation using the colormap
color_map = plt.cm.viridis
colors = color_map(np.linspace(0, 1, len(simulations)))
simulation_colors = {simulation: color for simulation, color in zip(simulations, colors)}

# Create subplots - one for each property
fig, axs = plt.subplots(len(properties), 1, figsize=(7, 4 * len(properties)))

# Define the bar width and offset increment
bar_width = 0.2  # Adjust as needed to fit your data
offset_increment = bar_width * len(simulations)

# Create a dictionary to store whether the legend label has been added for each simulation
added_legend = {simulation: False for simulation in simulations}

# Iterate over properties to create each subplot
for i, property in enumerate(properties):
    # Iterate over the BH Age points for consistent x-coordinates
    for age_idx, age in enumerate(ages):
        # For each age, iterate over simulations to get the respective y-coordinate
        for j, simulation in enumerate(simulations):
            # Filter the data for each simulation at the current age
            data = filt[(filt['simulation'] == simulation) & (filt['age_myr'] == age)]
            # Get the color for the current simulation
            color = simulation_colors[simulation]
            # Check if there is data to plot
            if not data.empty:
                # Calculate the bar position for the current simulation
                bar_position = age_idx + j * bar_width
                # Plot the bar if data is available
                bar = axs[i].bar(bar_position, data[property].values[0], width=bar_width, color=color, label=simulation if not added_legend[simulation] else "_nolegend_")
                # Add the simulation to the legend only if it hasn't been added before
                if not added_legend[simulation]:
                    added_legend[simulation] = True

    # Set subplot titles and labels
    axs[i].set_title(property)
    if property == 'bh_mass':
        axs[i].set_ylabel(r'$M_{BH} (M_{\odot})$')
        axs[i].set_ylim([60100, 64000])
    elif property == 'max_clump_mass':
        axs[i].set_ylabel('Max Clump Mass ($M_{\odot}$)')
       #axs[i].set_ylim([100, 3000])
    elif property == 'total_clump_mass':
        axs[i].set_ylabel('Total Clump Mass ($M_{\odot}$)')
        #axs[i].set_ylim([500, 6000])

# Configure x-axis for all subplots
labels = list(added_legend.keys())  # This will ensure the order of the keys is preserved
for ax in axs:
    ax.set_xticks(np.arange(len(ages)) + offset_increment / 2 - bar_width / 2)
    ax.set_xticklabels(ages, rotation=45)  # Rotate labels if they overlap

# Create the legend manually
handles = []
labels = []
for simulation in simulations:
    handles.append(plt.Rectangle((0, 0), 1, 1, color=simulation_colors[simulation]))
    labels.append(simulation)

# Add the legend to the last subplot
axs[-1].legend(handles, labels, loc='upper left', title='Simulation')

# Set the x-axis label for the last subplot
axs[-1].set_xlabel('BH Age [Myr]')

# Adjust the layout
plt.tight_layout()

# Save the plot
timestamp = datetime.now().strftime("%y-%m-%d-%H-%M")
plt.savefig(f'plots/clump_properties_time_series_barchart_{timestamp}.png')
print(f'Plot saved as plots/clump_properties_time_series_barchart_{timestamp}.png')
