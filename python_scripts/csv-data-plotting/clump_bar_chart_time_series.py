import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
This script makes a time series bar chart for clump properties using data from a CSV file.
The data is grouped by simulation and BH age, and each property is plotted as a bar chart over time.
The color of the bars is determined by the simulation, and a legend of simulations is added to the plot.
"""


# Load the data from the CSV file
df = pd.read_csv('clump_data.csv', comment='#')

# Sort the DataFrame by 'age_myr' to ensure the plots are in chronological order
df = df.sort_values('age_myr')

# Unique simulations and BH Age values for grouping
simulations = df['simulation'].unique()
ages = df['age_myr'].unique()
properties = ["max_clump_mass", "total_clump_mass", "bh_mass"]  # Replace with your actual properties

# Create a colormap object for 'viridis'
color_map = plt.cm.viridis

# Generate a color for each simulation using the colormap
colors = color_map(np.linspace(0, 1, len(simulations)))
simulation_colors = {simulation: color for simulation, color in zip(simulations, colors)}

# Create subplots - one for each property
fig, axs = plt.subplots(len(properties), 1, figsize=(7, 4 * len(properties)))

# Define the bar width and offset increment
bar_width = 0.2  # Adjust as needed to fit your data
offset_increment = bar_width * len(simulations)

# Create a dictionary to store whether the legend label has been added for each simulation
# Add legend to the top subplot using the handles we created
# Create a list of handles and labels for the legend, ordered by added_legend keys
handles = []
added_legend = {simulation: False for simulation in simulations}

# Iterate over properties to create each subplot
for i, property in enumerate(properties):
    # Set log scale for 'bh_mass' property
    if property == 'bh_mass':
        axs[i].set_yscale('log')

    # Iterate over the BH Age points for consistent x-coordinates
    for age_idx, age in enumerate(ages):
        # For each age, iterate over simulations to get the respective y-coordinate
        for j, simulation in enumerate(simulations):
            # Filter the data for each simulation at the current age
            data = df[(df['simulation'] == simulation) & (df['age_myr'] == age)]
            # Get the color for the current simulation
            color = simulation_colors[simulation]
            # Check if there is data to plot
            if not data.empty:
                # Calculate the bar position for the current simulation
                bar_position = age_idx + j * bar_width
                # Plot the bar if data is available
                bar = axs[i].bar(bar_position, data[property].values[0], width=bar_width, color=color)
                # Add the simulation to the legend only if it hasn't been added before
                if not added_legend[simulation]:
                    bar.set_label(simulation)
                    added_legend[simulation] = True

    # Set subplot titles and labels
    axs[i].set_title(property)
    axs[i].set_ylabel(r'$M_{BH} (M_{\odot})$') if property == 'bh_mass' else None
    axs[i].set_ylabel(r'Max Clump Mass $(M_{\odot})$') if property == 'max_clump_mass' else None
    axs[i].set_ylabel(r'Total Clump Mass $(M_{\odot})$') if property == 'total_clump_mass' else None
    axs[i].set_ylim([300, 1150]) if property == 'max_clump_mass' else None
    axs[i].set_ylim([1000, 3250]) if property == 'total_clump_mass' else None

# Configure x-axis for all subplots
labels = list(added_legend.keys())  # This will ensure the order of the keys is preserved
for ax in axs:
    ax.set_xticks(np.arange(len(ages)) + offset_increment / 2 - bar_width / 2)
    ax.set_xticklabels(ages, rotation=45)  # Rotate labels if they overlap
    # Force x-axis ticks to appear on all subplots
    for label in labels:
        for container in ax.containers:
            if container.get_label() == label:
                handles.append(container)


# Create the legend at the top subplot
axs[-1].legend(handles, labels, loc='upper left', title='Simulation')

# Set the x-axis label for the last subplot
axs[-1].set_xlabel('BH Age [Myr]')

# Adjust the layout
plt.tight_layout()

# Save the plot
plt.savefig('plots/clump_properties_time_series_barchart.png')
