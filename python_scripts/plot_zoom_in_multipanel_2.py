import os
import sys
import yt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from yt.utilities.math_utils import ortho_find
import re
from smartstar_find import ss_properties
from plot_multipanel_time_2 import create_axes_grid, get_min_max_values, configure_font, set_ticks_and_labels, configure_font
from plot_multi_projections import tidy_data_labels
from plot_disc_projections import _make_disk_L
from plot_zoom_in_multipanel import format_sci_notation


def field_from_projection_plot(field, ds, center, width_pc, north, dir, npixels=2048):
    """
    Compute field from a projection plot of a dataset.
    """
    dx = ds.index.get_smallest_dx().in_units('cm')
    p = yt.ProjectionPlot(ds, dir, ("gas", field), center=center, width=width_pc, north_vector=north)
    frb = p.frb[field]/dx
    return frb

# Make data array
root_dir = "/Backup00/sgordon/disk14/cirrus-runs-rsync/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim_names = ["1B.RSb01-2", "1B.RSb04", "1B.RSb16"]  # Simulation names
dds_paths = ["DD0138/DD0138", "DD0138/DD0138", "DD0167/DD0167"]  # Data dump paths for each simulation

# Set the field cmap limits
min_temp = 20
max_temp = 2000
min_density = 2e4
max_density = 1e10
npixels = 800

# Prepare the list for full simulation paths
DS = []
labels = []
print(f"Loading datasets for simulations: {', '.join(sim_names)}")
for sim_name, dds_path in zip(sim_names, dds_paths):
    ds_full_path = os.path.join(root_dir, sim_name, dds_path)
    print(f"Loading dataset {ds_full_path}")
    ds = yt.load(ds_full_path)
    label = tidy_data_labels(sim_name)
    DS.append(ds)
    labels.append(label)

# Define the figure size
configure_font()
fig = plt.figure(figsize=(11, 8.7))
nrows = 3
ncols = 4
gap = 0.038
plt.rcParams['text.usetex'] = True # use LaTeX for all text
fontsize = 12
 
# Define the size of each subplot
total_gap_width = (ncols - 1) * gap
size = 1.0 / nrows  # Since we have 3 rows, the height of each subplot is 1/3 of the figure height
im_width = (1 - total_gap_width) / ncols  # Since we have 4 columns, the width of each subplot is 1/4 of the figure width

# Loop over the datasets (rows) and zoom levels (columns)
for row in range(nrows):

    # Load the dataset and label
    ds = DS[row]
    label = labels[row]

    # find smallest cell width
    dx = ds.index.get_smallest_dx().in_units('pc')
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos

    # find L 
    disc_r_pc = disc_h_pc = 0.5
    _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

    # set north vector and v
    orient = "face-on"
    vecs = ortho_find(L)
    north = vecs[1] if orient=="face-on" else vecs[2]
    dir = vecs[0] if orient=="face-on" else vecs[1]
    
    for col in range(ncols): 
        # Calculate the position of the bottom left corner of each subplot
        left = col * (im_width + gap)
        bottom = 1 - (row + 1) * size  # Subtract from 1 because the y-axis starts from the top

        # Create an axis for each subplot
        ax = fig.add_axes([left+0.05, bottom, im_width, size-0.01], frame_on=True)
        ax = fig.add_axes([left+0.1, bottom, im_width, size-0.01], frame_on=True) if col == 2 or col == 3 else print("temperature cols")

        # col 1
        if col == 0:
            field = "temperature" 
            cmap = "RED TEMPERATURE"
            width = 3*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            # line integral of the temperature field along the line of sight, weighted by density
            temp= yt.off_axis_projection(disk, center, dir, width, npixels, ("gas", "temperature"), weight=("gas", "density")) 
            im1 = ax.imshow(temp, cmap=cmap, origin="lower", norm=LogNorm())
            im1.set_clim(min_temp, max_temp)

            # Cell width and simulation label
            im1.axes.set_title("dx = {} pc".format(format_sci_notation(float(dx))))
            im1.axes.text(0.07, 0.08, str(label), color='black', size=fontsize+4, transform=im1.axes.transAxes, 
                          bbox=dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="white", linewidth=3, alpha=0.5))

        # col 2    
        elif col == 1:
            field = "temperature" 
            cmap = "RED TEMPERATURE"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            # line integral of the temperature field along the line of sight, weighted by density
            temp= yt.off_axis_projection(disk, center, dir, width, npixels, ("gas", field), weight=("gas", "density")) 
            im2 = ax.imshow(temp, cmap=cmap, origin="lower", norm=LogNorm())
            im2.set_clim(min_temp, max_temp)

            # BH age label
            im2.axes.set_title("BH Age = {:.2f} Myr".format(ss_age[0] / 1e6))

        # col 3
        elif col == 2:
            field = "number_density" 
            cmap = "viridis"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            # line integral of the temperature field along the line of sight, weighted by density
            dens = yt.off_axis_projection(disk, center, dir, width, npixels, ("gas", field), weight=("gas", "density")) 
            im3 = ax.imshow(dens, cmap=cmap, origin="lower", norm=LogNorm())
            im3.set_clim(min_density, max_density)

            # BH mass and BH position label
            im3.axes.set_title(r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)))
            ax.scatter(center[0], center[1], color='black', marker="x", s=15, alpha=0.8, linewidths=0.2)

        # col 4
        elif col == 3:
            field = "cooling_time" # change to cooling length resolution
            cmap = "magma"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            cool = yt.off_axis_projection(disk, center, dir, width, npixels, ("gas", field), weight=("gas", "density")) 
            im4 = ax.imshow(cool, cmap=cmap, origin="lower", norm=LogNorm())
            #im4.set_clim(min_cooling, max_cooling)

        # Set ticks and ticklabels
        ax.set_xticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks
        ax.set_yticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks
        ticks = np.linspace(-float(width.d)/2, float(width.d)/2, num=5)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if row == 2:
            tick_labels = ['-1.5', '-0.75', '0.0', '0.75', '1.5']
            tick_labels2 = ['-0.15', '-0.075', '0.0', '0.075', '0.15']
            ax.set_xticklabels(tick_labels) if col == 0 else ax.set_xticklabels(tick_labels2)
            ax.set_xlabel("(pc)")

        # Set tick parameters to create 'inner' ticks, except on the bottom axis of last row
        ax.tick_params(which='minor', length=1)
        ax.tick_params(axis='x', which='both', direction='inout', bottom=True, top=True, length=2, width=1, colors='black')


# Add the colorbar to the left side of the figure
cax_width = 0.02
cax_length = 0.90
cax_left = fig.add_axes([0.01, 0.04, cax_width, cax_length])  # [left, bottom, width, height]
cbar_left = fig.colorbar(im1, cax=cax_left)  # Assuming im1 is one of your imshow panels with a colormap
cbar_left.ax.set_ylabel('Temperature (K)', fontsize=14, rotation=270, labelpad=15)
cbar_left.ax.yaxis.set_ticks_position('left') # Move the ticks to the left
cbar_left.ax.yaxis.set_label_position('left') # Move the label to the left

# Add the colorbar to the right side of the figure
cax_right = fig.add_axes([1.07, 0.04, cax_width, cax_length])
cbar_right = fig.colorbar(im3, cax=cax_right)  
cbar_right.ax.set_ylabel(r'Number Density \big($\rm \frac{1}{cm^{3}}$\big)', fontsize=14)

# Save the figure
plot_name_prefix = f"multiplot_axesgrid_zoom-in_3_fields"
plt.savefig('plots/' + plot_name_prefix + '.pdf', bbox_inches='tight')
print("created plots/" + str(plot_name_prefix)+ '.pdf')
plt.close()
