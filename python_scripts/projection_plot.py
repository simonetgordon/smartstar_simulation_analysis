"""
Plot projection of simulation in density, temperature or dark matter. Call like:
python -i projection_plot.py DD0133/DD0133
"""

import yt
import shutil
import sys
import os
import re # complex str searches
from matplotlib import ticker
from smartstar_find import ss_properties
from find_disc_attributes import _make_disk_L
from plot_multipanel_time_2 import configure_font, make_projection_plot, create_axes_grid, set_axis_labels_and_colorbar, \
    configure_projection_plot, get_min_max_values, set_ticks_and_labels
from plot_multi_projections import tidy_data_labels, find_north_vector
from yt.utilities.math_utils import ortho_find
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


def main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n_factor=2e5, max_n_factor=0.30, cmap="viridis"):

    configure_font()
    fig = plt.figure()
    grid = create_axes_grid(fig, nrows=3, ncols=2, dim=(0.1, 0.1, 0.53, 0.79), cbar_size="2%")

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
                dd = dds1[s]
            ds = yt.load(os.path.join(root_dir[l], sim[l], dd))
            label = tidy_data_labels(sim[l])
            DS.append(ds)
            LABEL.append(label)
    print("DS = ", DS)

    # loop over datasets 4 and make 2 projection plots per dataset
    for i, ds in enumerate(DS):
        print("ds = ", ds)
        print("i = ", i)

        # grab bh properties and make sphere centered on BH
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = ss_pos
        print("center = ", center)

        # make small disk data container and define angular momentum vector L and north vector
        disc_r_pc = disc_h_pc = 0.15
        _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
        vecs = ortho_find(L)
        
        for j in [0,1]:

            # set north vector and v
            north = vecs[0] if (j == 0) else vecs[2]
            v = 2 if (j == 0) else 0

            # make projection plot and configure basic parameters. remake bigger disk for projection
            width_pc = 30.24 * ds.units.pccm
            disc_r_pc = disc_h_pc = 2.1
            disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
            p = make_projection_plot(ds, width_pc, disk, L, field, vecs, v, north, min_n, max_n, fontsize, cmap=cmap, center=center)
            p = configure_projection_plot(p, field, min_n, max_n, fontsize)

            # ... pre-render plot customization ...

            # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[("gas", field)]
            plot.figure = fig
            print("k = ", i * 2 + j)
            k = i * 2 + j
            plot.axes = grid[k].axes

            # age, mass and simulation label in first and third columns
            if ((j == 0)):
                p.annotate_text((0.05, 0.9), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), 
                    coord_system="axis", text_args={"color": "white", "fontsize": 8})
                
                # make colorbar
                plot.cax = grid.cbar_axes[k]

            # simulation label in first row and every second column only
            if ((i == 0) and (j == 0)):
                p.annotate_text((0.07, 0.1), "{}".format(LABEL[i]), 
                    coord_system="axis", text_args={"color": "yellow", "fontsize": 8}, 
                    inset_box_args={
                                    "boxstyle": "square,pad=0.3",
                                    "facecolor": "black",
                                    "linewidth": 3,
                                    "edgecolor": "white",
                                    "alpha": 0.5,
                                })
            # mark BH position
            p.annotate_marker(center, coord_system="data", color="black", marker="x", plot_args={"s": 10, "alpha": 0.8,'linewidths': 0.1, 'edgecolors': 'black'})

            # render the plot
            p.render()

            # ... post render customization ...

            # Modify colorbar and axes properties **after** p.render() so that they are not overwritten.

            # ticks + ticklabels
            set_ticks_and_labels(grid, k, xticks, plot, no_xticklabels=False)

            # x and y extrenal axis labels and colorbar settngs
            set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, 
                                         max_n, fontsize, ticker)

    # save
    plot_name = 'projection-res-' + str(LABEL[i]) + '-' + str(field) + '-' + 's2' + '.pdf'
    p.save('plots/' + plot_name,)
    print("created plots/" + str(plot_name))
    plt.show()


if __name__ == "__main__":

    root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
                #"/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/"
                ]
    sim = ["2B.RSm04", 
           #"2B.m08-4dx"
           ]
    dds1=None
    #dds2 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.29, 079, 1.01 Myr for m01, 
    #xdds3 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.201, 0.801, 1.002 Myr for b01
    dds2 = ["DD0229/DD0229", "DD0268/DD0268", "DD0278/DD0278"]  # 0.3, 0.69, 0.79 Myr for 2B.m04, 
    dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0278/DD0278"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
    fontsize = 14
    width_pc = 1.5 * yt.units.pc # must change xticks if this changes
    xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    min_n_factor=2e5
    max_n_factor=0.1
    field = "number_density"
    cmap = 'viridis'

    main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n_factor, max_n_factor)


        # # Gas density
        # if field == "density":
        #     field = "number_density"
        #     p = yt.ProjectionPlot(ds, 'x', ("gas", field), width=width, center=center, data_source=sp,
        #                             weight_field='density')
        #     p.set_cmap(field, 'viridis')
        #     p.set_font_size(fontsize)
        #     p.set_background_color(("gas", field))
        #     p.set_axes_unit('pc')

        #     # annotate
        #     #p.annotate_scale(corner='lower_right')
        #     #p.annotate_streamlines(("gas", "velocity_y"), ("gas", "velocity_z"), density = 0.7, linewidth=0.6, color='yellow')
        #                             #field_color=("gas", "velocity_x")
        #                             # )
        #     p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=False)
        #     p.annotate_text((0.55, 0.94), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
        #                     text_args={"color": "white"})
        #     p.annotate_grids(min_level=18, cmap='Spectral') # not supported in OffAxisProjection
        #     p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
        #     p.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
        #     #p.annotate_cell_edges(line_width=0.00002, alpha=0.7, color='white')
        #     #p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
        #     p.annotate_title("High Resolution Region")

        #     # save
        #     plot_name = 'density-' + str(sim) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
        #     p.save('density_plots/' + plot_name)
        #     print("created density_plots/" + str(plot_name))
        #     field = "dm"
        #     plt.show()

        # # Temperature
        # elif field == "temperature":
        #     p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width=width, center=center, data_source=sp,
        #                             weight_field='density')
        #     p1.set_cmap('temperature', 'RED TEMPERATURE')
        #     plot = p1.plots[list(p1.plots)[0]]
        #     ax = plot.axes
        #     # nicen up the plot by setting the background color to the minimum of the colorbar
        #     p1.set_background_color(("gas", "temperature"))
        #     # hide the axes, while still keeping the background color correct:
        #     p1.hide_axes(draw_frame=True)
        #     p1.set_font_size(28)
        #     p1.set_axes_unit('pc')
        #     p1.annotate_scale(corner='lower_left')
        #     plot_name = 'temperature-' + str(root_dir[70:]) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
        #     p1.save('temperature_plots/' + plot_name)

        # elif field == "dm":
        #     plt.show()
        #     plt.figure()

        #     #p1 = yt.ProjectionPlot(ds, 'x', "all_cic", width=width, center=center, data_source=sp)
        #     p2 = yt.ParticlePlot(ds, ("all", "particle_position_y"), ("all", "particle_position_z"), ("all", "particle_mass"), 
        #                         width=width, center=center, data_source=sp)
        #     #p1.set_cmap("all_cic", 'viridis')
        #     p2.set_font_size(fontsize)
        #     p2.set_axes_unit('pc')
        #     p2.set_unit(("all", "particle_mass"), "Msun")

        #     p2.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
        #     p2.annotate_text((0.55, 0.94), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
        #                     text_args={"color": "black"})

        #     p2.annotate_marker(center, coord_system="data", color="black")  # mark bh position
        #     p1.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
        #     p2.annotate_title("DM in High Resolution Region")

        #     # save
        #     plot_name = 'dm-particles' + str(sim) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
        #     p2.save('density_plots/' + plot_name)
        #     print("created density_plots/" + str(plot_name))