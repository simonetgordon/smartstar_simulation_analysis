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
    grid = create_axes_grid(fig, nrows=3, ncols=2, dim=(0.1, 0.1, 0.52, 0.79), cbar_size="2%")

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

        # make small disk data container and define angular momentum vector L and north vector
        disc_r_pc = disc_h_pc = 0.1
        _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
        vecs = ortho_find(L)
        
        for j in [0,1]:

            # set north vector and v
            north = vecs[1] if (j == 0) else vecs[0]
            v = 0 if (j == 0) else 2

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
                    coord_system="axis", text_args={"color": "white", "fontsize": 10})
                
                # make colorbar
                plot.cax = grid.cbar_axes[k]

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
    dds2 = ["DD0229/DD0229", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.8 Myr for 2B.m04, 
    dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
    fontsize = 14
    width_pc = 1.5 * yt.units.pc # must change xticks if this changes
    xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    min_n_factor=2e5
    max_n_factor=0.1
    field = "number_density"
    cmap = 'viridis'
    slice

    main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n_factor, max_n_factor)