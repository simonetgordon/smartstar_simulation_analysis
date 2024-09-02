import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.units import pc, Msun, Myr
from yt.visualization.color_maps import register_yt_colormaps_from_cmyt
from helper_functions import ss_properties, setup_plot_env
from csv_data_plotting.helper_functions import extract_simulation_name
import logging
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from hist_scatter import scatter_hist2d

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Function to load datasets
def load_datasets(original_time_path, no_feedback_path, feedback_path):
    ds_og = yt.load(original_time_path)
    ds_no_feedback = yt.load(no_feedback_path)
    ds_with_feedback = yt.load(feedback_path)
    return ds_og, ds_no_feedback, ds_with_feedback

# Function to create region around the black hole
def define_region(ds, ss_pos, region_radius):
    return ds.sphere(ss_pos, region_radius)

# Function to create multi-panel plot
def create_multi_panel_plot(region, sim_name, age, region_radius, title_prefix):
    # Setup the plot environment
    setup_plot_env()
    fig = plt.figure(figsize=(10, 10))
    gs = GridSpec(3, 3, width_ratios=[2, 8, 0.4], height_ratios=[1, 5, 1])
    
    # # Main phase plot (density vs temperature vs mass) as matplotlib plot
    ax_main = fig.add_subplot(gs[1, 1])

    # Assuming region is a YT data object with the fields of interest
    temperature = region["gas", "temperature"].to("K").v
    density = region["gas", "H_nuclei_density"].to("cm**-3").v
    mass = region["gas", "cell_mass"].to("Msun").v

    # Define the number of bins and the limits for the colorbar (mass)
    num_bins = 200
    vmin = 1e-5
    vmax = 1e1

    # Logarithmic bin edges
    x_bin_edges = np.logspace(np.log10(density.min()), np.log10(density.max()), num_bins)
    y_bin_edges = np.logspace(np.log10(temperature.min()), np.log10(temperature.max()), num_bins)
    # Define LogNorm with the desired limits
    norm = LogNorm(vmin=vmin, vmax=vmax)

    # Plot the phase plot
    register_yt_colormaps_from_cmyt()
    plot = scatter_hist2d(density, temperature, bins=[x_bin_edges, y_bin_edges], weights=mass.T, s=5, 
                          cmap='arbre', ax=ax_main, norm=norm)

    # Add a logarithmically scaled colorbar
    cbar = fig.colorbar(plot, ax=ax_main, norm=norm, pad=0, aspect=40)  # Adjust the colorbar
    cbar.set_label(r"Cell Mass [$M_{\odot}$]")

    # Set ax_main's labels
    # ax_main.set_xlabel(r"H Nuclei Density [$\rm cm^{-3}$]")
    # ax_main.set_ylabel(r"Temperature [K]")

    # Set logarithmic scales
    ax_main.set_xscale('log')
    ax_main.set_yscale('log')

    # Set limits
    ax_main.set_xlim(2, 2e8)
    ax_main.set_ylim(1e1, 2e8)

    # Add title and annotations
    text = f"Radius = {region_radius[0]:.2f} pc"
    ax_main.text(0.95, 0.95, text, transform=ax_main.transAxes, fontsize=12, verticalalignment='top', horizontalalignment='right', 
                 bbox=dict(facecolor='grey', alpha=0.3))
    ax_main.set_title(f"{title_prefix}: {sim_name} at t = {age/1e6:.2f} Myr", loc='center')

    # Additional plot: mass vs temperature (left of y-axis)
    ax_left = fig.add_subplot(gs[1, 0], sharey=ax_main)
    ax_left.scatter(mass, temperature, marker='o', s=2, c='grey')
    ax_left.set_xlabel(r"Mass [$M_{\odot}$]")
    ax_left.set_ylabel(r"Temperature [K]")
    ax_left.set_xscale('log')
    ax_left.set_yscale('log')
    #ax_left.set_xlim(1e-6, 1e1)
    ax_left.grid(True)

    # Additional plot: mass vs H nuclei density (below x-axis)
    ax_bottom = fig.add_subplot(gs[2, 1], sharex=ax_main)
    ax_bottom.scatter(density, mass, marker='o', s=2, c='grey')
    ax_bottom.set_xlabel(r"H Nuclei Density [$\rm cm^{-3}$]")
    ax_bottom.set_ylabel(r"Mass [$M_{\odot}$]")
    ax_bottom.set_xscale('log')
    ax_bottom.set_yscale('log')
    ax_bottom.grid(True)
    aspect = 0.223 if title_prefix == "Feedback" else 'auto'
    ax_bottom.set_aspect(aspect, anchor= 'SW')  # Set the aspect of the axes and anchor the LHS to the bottom

    # Adjust the layout
    pos_main = ax_main.get_position()
    pos_bottom = ax_bottom.get_position()
    pos_left = ax_left.get_position()

    # Align the bottom plot's top edge with the main plot's bottom edge
    ax_bottom.set_position([pos_bottom.x0, pos_main.y0 - 2*pos_bottom.height, pos_bottom.width, pos_bottom.height])

    # Align the left plot's right edge with the main plot's left edge
    ax_left.set_position([pos_main.x0 - 2*pos_left.width, pos_left.y0, pos_left.width, pos_left.height])
    ax_bottom.set_aspect(aspect, anchor= 'SW')
    
    fig.tight_layout(w_pad=0, h_pad=0)
    return fig


# Function to save the figure
def save_plot(fig, filename):
    fig.savefig(filename)
    logger.info(f"Saved plot to {filename}")

# Filepaths for the datasets
original_time_path = "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0445/DD0445/"
no_feedback_path = "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0447/DD0447/"
feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0695/DD0695/"
#feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/fb-radius=7dx/DD0695/DD0695/"
#feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0695/DD0695/"
#feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695/"
feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2783/DD2783/"
feedback_path = "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0652/DD0652/"

# Simulation names
sim_no_fb = extract_simulation_name(no_feedback_path, custom_name=None)
sim_with_fb = extract_simulation_name(feedback_path, custom_name="1B.resim.th.b01-3-eps-0.0001")
sim_og = extract_simulation_name(original_time_path, custom_name=None)

# Print simulation names
logger.info("Simulation Names:")
logger.info(f"Original: {sim_og}")
logger.info(f"No Feedback: {sim_no_fb}")
logger.info(f"With Feedback: {sim_with_fb}")

# Load the datasets
ds_og, ds_no_feedback, ds_with_feedback = load_datasets(original_time_path, no_feedback_path, feedback_path)

# Extract bh properties
ss_pos_no_fb, ss_mass, ss_age = ss_properties(ds_no_feedback)
ss_pos_og, _, ss_age_og = ss_properties(ds_og)
ss_pos_fb, _, ss_age_fb = ss_properties(ds_with_feedback)

# Define the region radius
region_radius = (0.5, 'pc')

# Define regions around the black hole
region_no_feedback = define_region(ds_no_feedback, ss_pos_no_fb, region_radius)
region_with_feedback = define_region(ds_with_feedback, ss_pos_fb, region_radius)
region_og = define_region(ds_og, ss_pos_og, region_radius)

# Create multi-panel plots
fig_og = create_multi_panel_plot(region_og, sim_og, ss_age_og[0], region_radius, "Original")
fig_no_feedback = create_multi_panel_plot(region_no_feedback, sim_no_fb, ss_age[0], region_radius, "No Feedback")
fig_with_feedback = create_multi_panel_plot(region_with_feedback, sim_with_fb, ss_age_fb[0], region_radius, "Feedback")

# Save the plots
save_plot(fig_og, f'paper_2_other/plot_gas_phases_og_{region_radius[0]:.1f}pc_{ss_age_og[0]/1e6:.2f}Myr_{sim_og}.png')
save_plot(fig_no_feedback, f'paper_2_other/plot_gas_phases_no_feedback_{region_radius[0]:.1f}pc_{ss_age[0]/1e6:.2f}Myr_{sim_no_fb}.png')
save_plot(fig_with_feedback, f'paper_2_other/plot_gas_phases_with_feedback_{region_radius[0]:.1f}pc_{ss_age_fb[0]/1e6:.2f}Myr_{sim_with_fb}.png')
