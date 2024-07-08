import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from helper_functions import remove_values, grow_black_hole, redshift_to_time


def load_mass_vs_redshift_data(file_path):
    """
    Load mass vs redshift data from a file and organize it into separate DataFrames.

    The function reads a file containing mass vs redshift data, organized into sections
    identified by comments starting with '#'. Each section is stored in a separate DataFrame.

    Parameters:
    file_path (str): The path to the input file containing mass vs redshift data.
                     The file should be formatted as follows:
                     - Lines starting with '#' denote the start of a new section.
                     - Data lines contain comma-separated values of redshift and mass.

    Returns:
    dict: A dictionary where the keys are section names (str) and the values are DataFrames (pd.DataFrame)
          containing the mass vs redshift data for each section. Each DataFrame has two columns:
          'z' (redshift) and 'mass'.
    """
    # Dictionary to hold dataframes
    dataframes = {}
    
    # Read the file content to identify the sections
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Temporary storage for data
    temp_data = []
    current_section = None
    
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            if current_section is not None and temp_data:
                print(f"Creating DataFrame for section: {current_section}")
                dataframes[current_section] = pd.DataFrame(temp_data, columns=['z', 'mass'])
                temp_data = []
            current_section = line[2:].strip()
            print(f"New section: {current_section}")
        elif line and not line.startswith('#'):
            if current_section is not None:
                z, mass = map(float, line.split(','))
                temp_data.append((z, mass))
                print(f"Adding data to section {current_section}: z={z}, mass={mass}")
    
    # Add the last section
    if current_section is not None and temp_data:
        print(f"Creating DataFrame for section: {current_section}")
        dataframes[current_section] = pd.DataFrame(temp_data, columns=['z', 'mass'])
    
    return dataframes


# Customize the plot
plt.xlabel('redshift')
plt.ylabel(r'$\log(M_{BH}/M_{\odot})$')
ax1_lim = 5.6, 25
plt.xlim(ax1_lim)
plt.ylim(1, 10)
plt.gca().invert_xaxis()  # Flip the x-axis

# Secondary x-axis for time
ax2 = plt.twiny()
ax2.set_xlim(ax1_lim)
ax2.set_xticks([18.509, 13.880, 11.285, 9.589, 8.376,  6.738, 5.666]) # 0.2, 0.3, 0.4,.. 7.46=0.7, 6.153=0.9
ax2_ticks = np.round(np.arange(0.2, 1.0001, 0.1, dtype='float16'), decimals=1)
to_remove = [0.7, 0.9]  
ax2_ticks = remove_values(ax2_ticks, to_remove)
ax2.set_xticklabels(ax2_ticks)
ax2.set_xlabel('t (Gyr)')
ax2.invert_xaxis()  # Flip the secondary x-axis

## Plot Edd rates ##
end_time_gyr = 1.01

# Define scenarios for light and heavy black holes
scenarios = [
    {'label': 'Light ($f_{Edd} = 1$)', 'color': 'moccasin', 'masses': [1e1, 1e2], 'fedd': 1, 'start_time': redshift_to_time(25)},
    {'label': 'Light ($f_{Edd} = 2$)', 'color': 'lightsalmon', 'masses': [1e1, 1e2], 'fedd': 2, 'start_time': redshift_to_time(25)},
    {'label': 'Heavy ($f_{Edd} = 1$)', 'color': 'powderblue', 'masses': [1e4, 1e5], 'fedd': 1, 'start_time': redshift_to_time(17.5)},
    {'label': 'Heavy ($f_{Edd} = 2$)', 'color': 'aquamarine', 'masses': [1e4, 1e5], 'fedd': 2, 'start_time': redshift_to_time(17.5)}
]

for scenario in scenarios:
    start_mass1, start_mass2 = scenario['masses']
    redshifts, masses1 = grow_black_hole(start_mass1, scenario['start_time'], end_time_gyr, time_steps=1000, fedd=scenario['fedd'])
    redshifts, masses2 = grow_black_hole(start_mass2, scenario['start_time'], end_time_gyr, time_steps=1000, fedd=scenario['fedd'])
    #plt.plot(redshifts, masses2, color=scenario['color'])

    plt.fill_between(redshifts, masses1, masses2, color=scenario['color'], alpha=0.5, label=scenario['label'])


## Observations ##

# Colors and markers for each dataset
styles = {
    "JWST": {"color": "crimson", "marker": "D"},
    "ALMA": {"color": "seagreen", "marker": "s"},
    "GN-z11": {"color": "darkviolet", "marker": "*"},
    "UHZ1": {"color": "magenta", "marker": "o"},
    "GHZ9": {"color": "black", "marker": "+"}
}

file_path = 'mass_vs_redshift_data.csv'  # JWST, ALMA, GN-z11, UHZ1, GHZ9 datapoints
dataframes = load_mass_vs_redshift_data(file_path)
for key, df in dataframes.items():
    print(f"DataFrame: {key}")
    plt.scatter(df['z'], df['mass'], label=key, color=styles[key]["color"], marker=styles[key]["marker"])


# Save
plt.legend()
figname = "plots/mass_vs_redshift_observations.png"
plt.savefig(figname)
print(f"saved to {figname}")
