"""
Plot projection of simulation in density, temperature or dark matter. Call like:
python -i plot_multiproj_2x3.py
"""

import yt
import shutil
import sys
import os
import re # complex str searches
from matplotlib import ticker
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from plot_multipanel_time_2 import configure_font, make_projection_plot, create_axes_grid, set_axis_labels_and_colorbar, \
    configure_projection_plot, get_min_max_values, set_ticks_and_labels, make_slice_plot
from plot_multi_projections import tidy_data_labels, find_north_vector
from yt.utilities.math_utils import ortho_find
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time
import numpy as np
import matplotlib as mpl


def plot_cbar_on_right(fig, grid, k, p, cmap):
    # Manually adding another colorbar on the right
    # Get position of the existing colorbar
    cbar_position = grid.cbar_axes[k].get_position()

    # Define new position for the second colorbar
    new_cbar_position = [cbar_position.x1 + 0.05, cbar_position.y0, 0.05, cbar_position.height]

    # Add new axes for the second colorbar
    cax = fig.add_axes(new_cbar_position)

    # Create the second colorbar
    cbar = plt.colorbar(p, cax=cax, orientation='vertical', cmap=cmap)


def main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n_factor=2e5, max_n_factor=0.30, cmap="viridis", slice=False):

    configure_font()
    fig = plt.figure()
    grid = create_axes_grid(fig, nrows=3, ncols=2, dim=(0.1, 0.12, 0.435, 0.77), cbar_size="2%")

    # get min and max values for colorbar
    min_n, max_n, min_n_index, max_n_index = get_min_max_values(root_dir, sim, dds2, min_n_factor=min_n_factor, max_n_factor=max_n_factor)


    # load datasets into list DS and labels into list LABEL
    DS, LABEL = [], []
    for l in range(len(sim)):
        for s, _ in enumerate(dds2):
            print("s = ", s)
            if "m04" in sim[l]:
                dd = dds2[s]
                print("dd = ", dd)
            elif "m08" in sim[l]:
                dd = dds3[s]
            else:
                dd = dds2[s]
            ds = yt.load(os.path.join(root_dir[l], sim[l], dd)) # ds = yt.load(os.path.join(root_dir, sim, dd))
            label = tidy_data_labels(sim[l])
            DS.append(ds)
            LABEL.append(label)
    print("DS = ", DS)

    # loop over datasets 4 and make 2 projection plots per dataset
    for i, ds in enumerate(DS):
        print("ds = ", ds)
        print("i = ", i)

        # grab bh properties and make sphere centered on BH
        ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
        center = ss_pos

        # make small disk data container and define angular momentum vector L and north vector
        disc_r_pc = disc_h_pc = 0.05
        _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
        vecs = ortho_find(L)
        dir = vecs[0]
        north = vecs[1]

        for j in [0,1]:

            # make projection plot and configure basic parameters. remake bigger disk for projection
            disc_r_pc = disc_h_pc = 2.1
            disk = ds.disk(center, L, disc_r_pc, disc_h_pc)
            disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
            if (slice == True):
                if (j == 0):
                    cmap='viridis'
                    field = 'number_density'
                    min_n=9e3
                    max_n=4e9
                    p = yt.SlicePlot(ds, dir, ("gas", field), center=center, width=(width_pc, "pc"), north_vector=north)
                    p.set_cmap(field, cmap)
                    p.set_zlim(field, min_n, max_n)
                    #p.hide_colorbar()
                    if (i == 0):
                        p.annotate_title("Number Density (cm$^{-3}$)")
                else:
                    cmap = 'magma' # 'octarine' 
                    field = 'velocity_cylindrical_theta' # 'velocity_cylindrical_radius'
                    field = 'velocity_cylindrical_z'
                    min_n = -20 # -5 for theta
                    max_n = 20
                    p = yt.SlicePlot(ds, dir, ("gas", field), center=center, width=(width_pc, "pc"), north_vector=north)
                    p.set_log(("gas", "velocity_cylindrical_z"), False)
                    p.set_unit("velocity_cylindrical_z", "km/s")
                    p.set_cmap(field, cmap)
                    p.set_zlim(field, min_n, max_n)
                    if (i == 0):
                        p.annotate_title("Disk z-Velocity (km/s)")
                    #p.hide_colorbar()
            else:
                p = make_projection_plot(ds, width_pc, disk, L, field, vecs, v, north, min_n, max_n, fontsize, cmap=cmap, center=center)
                p = configure_projection_plot(p, field, cmap, min_n, max_n, fontsize=14)
            
            # define axis grid index
            k = i * 2 + j
            print("k = ", i * 2 + j)

            # ... pre-render plot customization ...

            # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            # plot = p.plots[("gas", field)]
            # plot.figure = fig
            # plot.axes = grid[k].axes
            # if k == 1:
            #     print("plotting colorbar on right")
            #     mappable = plot.image
            #     cax = grid.cbar_axes[k]
            #     cax.colorbar(mappable)
            #     print("cax = {}, mappable = {}".format(cax, mappable))

            # age, mass and simulation label in first and third columns
            if (j == 0):
                p.annotate_text((0.07, 0.88), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), 
                    coord_system="axis", text_args={"color": "white", "fontsize": 10})

            # simulation label in first row and every second column only
            if ((i == 0) and (j == 0)):
                p.annotate_text((0.07, 0.1), "{}".format(LABEL[i]), 
                    coord_system="axis", text_args={"color": "yellow", "fontsize": 10}, 
                    inset_box_args={
                                    "boxstyle": "square,pad=0.3",
                                    "facecolor": "black",
                                    "linewidth": 3,
                                    "edgecolor": "white",
                                    "alpha": 0.5,
                                })
            # mark BH position
            p.annotate_marker(center, coord_system="data", color="black", marker="x", plot_args={"s": 10, "alpha": 0.8,'linewidths': 0.1, 'edgecolors': 'black'})

            # keep colorbar activation here, right before rendering - otherwise it doesn't show!
            plot = p.plots[("gas", field)]
            plot.figure = fig
            plot.axes = grid[k].axes
            print("plotting colorbar on right")
            plot.cax = grid.cbar_axes[k]

            # Assuming you're using two columns
            if j == 0:
                cmap_column0 = mpl.cm.viridis  # replace with the actual colormap of column 0
                norm_column0 = mpl.colors.LogNorm(vmin=min_n, vmax=max_n)  # replace with actual normalization
            else:
                cmap_column1 = mpl.cm.magma  # replace with the actual colormap of column 1
                norm_column1 = mpl.colors.Normalize(vmin=min_n, vmax=max_n)  # replace with actual normalization

            # render the plot
            p.render()

            # ... post render customization ...

            # Modify colorbar and axes properties **after** p.render() so that they are not overwritten.
            # if k == 5:
            #     # make colorbar
            #     print("making colorbar")
            #     #plt.colorbar(im, cax=cax)  # im = p.plots[("gas", field)]._plot_obj
            #     plot.cax = grid.cbar_axes[k]
            #     grid.cbar_axes[k].minorticks_on()
            #     # move colorbar and ticks to left hand side
            #     #grid.cbar_axes[k].yaxis.set_label_position('left')
            #     #grid.cbar_axes[k].yaxis.tick_left()
            #     grid.cbar_axes[k].tick_params(labelsize=14)
            #     grid.cbar_axes[k].set_ylabel(r'Disk Radial Velocity (km/s)')

            #     #min_n_index, max_n_index = np.log10(min_n), np.log10(max_n)
            #     # set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, 
            #     #                              max_n, fontsize, ticker)

            # ticks + ticklabels
            set_ticks_and_labels(grid, k, xticks, plot, no_xticklabels=False)

            # x and y extremal axis labels and colorbar settngs
            if j == 0:
                grid[k].axes.set_ylabel("{:.2f} Myr".format((ss_age/1e6)[0]))
            
        # After the loop, create and position the colorbars
        cbar_ax1 = fig.add_axes([0.1, 0.92, 0.21, 0.02])  # Adjust these values
        cbar_ax2 = fig.add_axes([0.32, 0.92, 0.21, 0.02])  # Adjust these values

        # Create colorbars using the stored colormaps and normalizations
        cb1 = mpl.colorbar.ColorbarBase(cbar_ax1, cmap=cmap_column0, norm=norm_column0, orientation='horizontal', ticklocation='top')
        cb2 = mpl.colorbar.ColorbarBase(cbar_ax2, cmap=cmap_column1, norm=norm_column1, orientation='horizontal', ticklocation='top')

        font_properties = {'size': '11', 'weight': 'normal', 'family': 'sans-serif'}  # Adjust as necessary

        # Set labels for the colorbars
        # cb1.set_label(r'Number Density $cm^{-3}$', **font_properties, labelpad=-40)  # Replace 'Density' with the label appropriate for the first column
        # cb2.set_label(r'Disk z-Velocity $km/s$', **font_properties, labelpad=-40)  # Replace 'Velocity' with the label appropriate for the second column



    # save
    plot_name = 'projection-res-' + str(LABEL[i]) + '-' + str(field) + '-s1' + '-zvel2.pdf'
    p.save('plots/' + plot_name,)
    print("created plots/" + str(plot_name))
    plt.show()


if __name__ == "__main__":

    root_dir = [#"/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
                #"/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/"
                "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/"

                ]
    sim = [#"2B.RSm04", 
           #"2B.m08-4dx"
           "1B.m16-4dx"
           ] 
    dds1=None
    dds2 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.29, 079, 1.01 Myr for m01, - always need to define dds2
    #xdds3 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.201, 0.801, 1.002 Myr for b01
    #dds2 = ["DD0229/DD0229", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.8 Myr for 2B.m04, 
    dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
    dds2 = ["DD0167/DD0167", "DD0178/DD0178", "DD0188/DD0188"]  # 0.3, 0.69, 0.8 Myr for 1B.m16,
    fontsize = 14
    width_pc = 0.18 # must change xticks if this changes
    xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    xticks = [-0.1, -0.05, 0.0, 0.05, 0.1]
    xticks = [-0.05, 0.0, 0.05]
    min_n=2e5
    max_n=0.1
    field = "number_density"
    cmap = 'viridis'
    slice = True

    main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n, max_n, cmap, slice)