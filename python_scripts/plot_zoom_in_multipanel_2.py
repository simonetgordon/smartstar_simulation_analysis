import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from yt.utilities.math_utils import ortho_find
import re
from smartstar_find import ss_properties
from plot_multipanel_time_2 import create_axes_grid, get_min_max_values, configure_font, set_ticks_and_labels, make_projection_plot
from plot_multi_projections import tidy_data_labels
from plot_disc_projections import _make_disk_L


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
min_density = 1
max_density = 2e3
min_temp = 9e3
max_temp = 1e10
npixels = 1024

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
fig = plt.figure(figsize=(12, 12))
nrows = 3
ncols = 4
plt.rcParams['text.usetex'] = True # use LaTeX for all text
fontsize = 12
 
# Define the size of each subplot
size = 1.0 / 3  # Since we have 3 rows, the height of each subplot is 1/3 of the figure height
im_width = 1.0 / 4  # Since we have 4 columns, the width of each subplot is 1/4 of the figure width

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
        left = col * im_width
        bottom = 1 - (row + 1) * size  # Subtract from 1 because the y-axis starts from the top

        # Create an axis for each subplot
        ax = fig.add_axes([left, bottom, im_width, size-0.01], frame_on=True)

        # col 1
        if col == 0:
            field = "temperature" 
            cmap = "RED TEMPERATURE"
            width = 3*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            #temp = field_from_projection_plot(field, ds, center, width, north, L, npixels=npixels)
            temp_image = yt.off_axis_projection(disk, center, L, width, npixels, ("gas", "temperature"))
            slc =yt.SlicePlot(ds, dir, ("index", "dx"), center=center, width=(float(width.d), 'pc'), north_vector=north, data_source=disk)
            slc_frb = slc.data_source.to_frb((float(width.d), "pc"), npixels)
            dx_image = slc_frb[("index", "dx")].to("cm")
            temp = temp_image/dx_image
            im1 = ax.imshow(temp, cmap=cmap, origin="lower", norm=LogNorm())
            im1.set_clim(min_temp, max_temp)

            # Cell width and simulation label
            im1.axes.set_title("dx = {} pc".format(format_sci_notation(float(dx))))
            im1.axes.text(0.07, 0.08, str(label), color='black', size=fontsize+2, 
                        transform=im1.axes.transAxes, bbox=dict(boxstyle="square,pad=0.3", 
                        facecolor="white", edgecolor="white", linewidth=3, alpha=0.5))
        elif col == 1:
            field = "temperature" 
            cmap = "RED TEMPERATURE"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            temp = field_from_projection_plot(field, ds, center, width, north, L, npixels=npixels)
            im2 = ax.imshow(temp, cmap=cmap, origin="lower", norm=LogNorm())
            im2.set_clim(min_temp, max_temp)
            df2 = disk.to_dataframe([("gas", "number_density"), ("gas", "temperature")])

            # BH age label
            im2.axes.set_title("BH Age = {:.2f} Myr".format(ss_age[0] / 1e6))

        elif col == 2:
            field = "number_density" # change to cooling length resolution
            cmap = "viridis"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            dens = field_from_projection_plot(field, ds, center, width, north, L, npixels=npixels)
            im3 = ax.imshow(dens, cmap=cmap, origin="lower", norm=LogNorm())
            im3.set_clim(min_density, max_density)

            # BH mass label
            im3.axes.set_title(r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)))

        elif col == 3:
            field = "cooling_time" # change to cooling length resolution
            cmap = "magma"
            width = 0.35*yt.units.pc
            height = width/2
            disk = ds.disk(center, L, width, height)
            cool = field_from_projection_plot(field, ds, center, width, north, L, npixels=npixels)
            im4 = ax.imshow(cool, cmap=cmap, origin="lower", norm=LogNorm())
            #im4.set_clim(min_cooling, max_cooling)


# Add the colorbar to the left side of the figure
cax_left = fig.add_axes([0.05, 0.1, 0.02, 0.8])  # [left, bottom, width, height]
cbar_left = fig.colorbar(im1, cax=cax_left)  # Assuming im1 is one of your imshow panels with a colormap
cbar_left.ax.set_ylabel('Temperature (K)', rotation=270, labelpad=15)

# Add the colorbar to the right side of the figure
cax_right = fig.add_axes([0.93, 0.1, 0.02, 0.8])  # Adjust these values as needed
cbar_right = fig.colorbar(im3, cax=cax_right)  # Assuming im1 is one of your imshow panels with a colormap
cbar_right.ax.set_ylabel('Number Density (1/cm^3)')

# Save the figure
plt.tight_layout()
plot_name_prefix = f"multiplot_axesgrid_zoom-in_3_fields"
plt.savefig('plots/' + plot_name_prefix + '.pdf', bbox_inches='tight')
print("created plots/" + str(plot_name))
plt.close()
