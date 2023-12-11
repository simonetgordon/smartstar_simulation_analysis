import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from yt.utilities.math_utils import ortho_find
import re
from smartstar_find import ss_properties
from plot_multipanel_time_2 import create_axes_grid, get_min_max_values, configure_font, make_projection_plot, set_ticks_and_labels
from plot_disc_projections import _make_disk_L
from contextlib import redirect_stdout

def tidy_data_labels(labels, custom_name=None):
    if custom_name:
        return custom_name
    
    # Function to apply replacements to a single label
    def apply_replacements(label):
        label = label.replace("-2", "")
        label = label.replace("RS", "")
        label = label.replace("-4dx", "")
        return label

    # Check if labels is a list or a single label
    if isinstance(labels, list):
        # Process each label in the list
        data_labels = [apply_replacements(label) for label in labels]
    elif isinstance(labels, str):
        # Process a single label
        data_labels = apply_replacements(labels)
    else:
        raise TypeError("labels should be a string or a list of strings")

    return data_labels


def first_index(a, val, rtol=0.1, atol=10):
    return next(j for j, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


def format_sci_notation(x):
    a, b = '{:.2e}'.format(x).split('e')
    return r'$\rm {} \times 10^{{{}}}$'.format(a, int(b))


def main(root_dir, sim, dds, field, k, widths_pccm, fontsize, min_n_factor, max_n_factor, orient="face-on", cmap="viridis", seed="s1"):
    configure_font()
    fig = plt.figure()
    grid = create_axes_grid(fig, nrows=3, ncols=1, dim=(0, 0, 0.3, 0.85), 
                            axes_pad=0.275, share_all=False, cbar_location="left", cbar_size="3%", cbar_pad="0%")

    min_n, max_n, min_n_index, max_n_index = get_min_max_values(root_dir, sim, dds, field, min_n_factor=min_n_factor, max_n_factor=max_n_factor)

    for i, dd in enumerate(dds):

        # load and label ds
        ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
        label = tidy_data_labels(sim[i])

        # find bh properties
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = ss_pos

        # define width
        widths = widths_pccm * ds.units.pccm
        k = int(sys.argv[-1])
        width = widths[k]

        # make disk data (with width) container and define angular momentum vector L
        disc_r_pc = disc_h_pc = 0.15
        disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

        # set north vector and v
        vecs = ortho_find(L)
        north = vecs[1] if orient=="face-on" else vecs[2]
        v = 0 if orient=="face-on" else 0
        # if north[-1] < 0:
        #     north *= [1,1,-1]

        print("min_n: {}, max_n: {}".format(min_n, max_n))
        # make projection plot
        disc_r_pc = disc_h_pc = 250
        disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
        p = make_projection_plot(ds, width, disk, L, field, vecs, v, north, min_n, max_n, fontsize, cmap)

        # find smallest cell width
        dx = ds.index.get_smallest_dx().in_units('pc')

        # This section takes care of annotating the projections based on its column
        if k == 0:
            p.annotate_title("dx = {} pc".format(format_sci_notation(float(dx))))
            p.annotate_text([0.07, 0.08], str(label), coord_system="axis", text_args={"color": "black", "size": fontsize+2},
                            inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                            "edgecolor": "white", "alpha": 0.5})
        elif k == 1:
            p.annotate_title("BH Age = {:.2f} Myr".format(ss_age[0] / 1e6))
        elif k == 2:
            #p.annotate_marker(center, coord_system="data", color="black")
            p.annotate_title(r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)))

        # ... Code related to scalebars, plotting and formatting ...

        # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[("gas", field)]
        plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

        # redraw the plot
        p.render()

        # axes ticks
        xticks = grid[0].axes.get_xticks()[1:-1]
        set_ticks_and_labels(grid, i, xticks, plot, no_xticklabels=True)

        # ticklabels
        if k == 0:
            grid[i].axes.set_xticklabels([str(int(x)) for x in xticks])
        elif k == 2:
            grid[i].axes.set_xticklabels(["-0.15", "-0.1", "-0.05", "0", "0.05", "0.1", "0.15"])
        elif k == 1:
            grid[i].axes.set_xticklabels(["-1", "-0.5", "0", "0.5", "1"])

        grid[i].axes.set_ylabel("")
        grid[i].axes.set_xlabel("(pc)")

        # Modify the scalebar, colorbar and axes properties **after** p.render() so that they are not overwritten.
        if (i == 2) and (k != 0): # i=1 -> i=2 temp
            plot.axes.get_yaxis().set_units('pc')

        if k == 0: 
            if field == "temperature":
                grid.cbar_axes[i].set_ylabel(r'Temperature \big($\rm K$\big)', fontsize=14)
                grid.cbar_axes[i].minorticks_on()
                grid.cbar_axes[i].tick_params(labelsize=14)

                # move colorbar and ticks to left hand side
                grid.cbar_axes[i].yaxis.set_label_position('left')
                grid.cbar_axes[i].yaxis.tick_left()
            else:
                grid.cbar_axes[i].set_ylabel(r'Number Density \big($\rm \frac{1}{cm^{3}}$\big)', fontsize=14)
                grid.cbar_axes[i].minorticks_on()

                # make minorticks
                start_value = np.log10(min_n)
                a1 = np.logspace(min_n_index, max_n_index, np.abs(min_n_index - max_n_index)+1)
                a2 = np.arange(1, 10, 1)
                minorticks = np.outer(a1, a2).flatten()
                atol = 10**(max_n_index-1)
                end_i = first_index(minorticks, float(max_n.d), rtol=0.1, atol=atol)
                minorticks = minorticks[:end_i]
                grid.cbar_axes[i].yaxis.set_minor_locator(ticker.AutoMinorLocator())
                grid.cbar_axes[i].set_yticks(minorticks, minor=True)
                grid.cbar_axes[i].tick_params(labelsize=14)

    plot_name_prefix = f"multiplot_axesgrid_{widths.d[k]}pccm_{field}_{orient}_{seed}"
    plt.savefig('plots/' + plot_name_prefix + '.pdf', bbox_inches='tight')
    print("created plots/{}.pdf".format(plot_name_prefix))

if __name__ == "__main__":

    # call like: python plot_zoom_in_multipanel.py 0/1/2 [column number]

    # Initialisation
    #root_dir = ["/cephfs/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"] * 3
    root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
                "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
                "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/",]
    sim = ["1B.RSb01-2", "1B.RSb04", "1B.RSb16"]# s1
    #sim = ["2B.RSb01", "2B.RSb04", "2B.RSb08"]# s2
    #sim = ["2B.RSm01", "2B.RSm04", "2B.m08-4dx"]# s2
    dds = ["DD0138/DD0138", "DD0138/DD0138", "DD0167/DD0167"] # s1
    #dds = ["DD0138/DD0138", "DD0138/DD0138", "DD0534/DD0534"] # s2 BHL at t = 1 Myr
    #dds = ["DD0206/DD0206", "DD0278/DD0278", "DD0277/DD0277"] # s2 MF at t = 0.79 Myr

    # column width in pc
    k = int(sys.argv[-1])
    widths_pccm = [6500, 60, 7] # 8500, 80, 10 for s1, 7 for s2
    fontsize = 12

    # set field, cmap
    field = "temperature" # "temperature" or "number_density
    cmap = "RED TEMPERATURE" # "RED TEMPERATURE" or "hot" or "inferno" or "magma" (temp) or "viridis" (density)

    # set to reduce colorbar limits
    min_n_factor = 50 # 50 for temp, 1e3 for density in s1, 4e4 for density,  in s2
    max_n_factor = 0.03 # 0.1 for temp+density s1, 0.00009 for dens s2

    main(root_dir, sim, dds, field, k, widths_pccm, fontsize, min_n_factor=min_n_factor, max_n_factor=max_n_factor, orient="face-on", cmap=cmap, seed="s2")
