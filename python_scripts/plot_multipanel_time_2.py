import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import pyplot, colors, ticker
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.utilities.math_utils import ortho_find
from smartstar_find import ss_properties
from plot_multi_projections import tidy_data_labels, first_index
from plot_disc_projections import _make_disk_L


def load_datasets(root_dir, sim, dds1, dds2, dds3):
    DS, LABEL = [], []
    for l in range(len(sim)):
        for s, _ in enumerate(dds1):
            if "2B.RSb" in sim[l]:
                dd = dds3[s]
            elif "2B.RSm" in sim[l]:
                dd = dds2[s]
            else:
                dd = dds1[s]
            ds = yt.load(os.path.join(root_dir[l], sim[l], dd))
            label = tidy_data_labels(sim[l])
            DS.append(ds)
            LABEL.append(label)
    return DS, LABEL


def configure_plots():
    pyplot.rcParams['font.size'] = 14
    pyplot.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"


def make_projection_plot(ds, width_pc, disk, L, field, vecs, v, north, min_n, max_n, fontsize):
    p = yt.ProjectionPlot(ds, vecs[v], ("gas", field), weight_field=("gas", "density"), north_vector=north, 
                          center=disk.center, width=width_pc, data_source=disk)
    p.set_axes_unit('pc')
    p.set_font_size(fontsize)
    p.set_cmap(field, 'viridis')
    p.set_zlim(("gas", field), min_n, max_n)
    p.hide_colorbar()
    return p


def create_axes_grid(fig, dds1, sim):
    return AxesGrid(
        fig,
        (0.01, 0.01, 0.76, 1.152),
        nrows_ncols=(len(dds1 * 2), len(sim)),
        axes_pad=0,
        label_mode="L",
        aspect=False,
        share_all=True,
        cbar_location="right",
        cbar_mode="single",
        cbar_size="2%",
        cbar_pad="0%",
    )


def get_min_max_values(root_dir, sim, dds2, min_n_factor=2e5, max_n_factor=0.30):
    field = "number_density"
    ds_hr = yt.load(os.path.join(root_dir[-1], sim[-1], dds2[-1]))
    min_n = ds_hr.r[("gas", field)].min() * 2e5
    max_n = ds_hr.r[("gas", field)].max() * 0.30
    return min_n, max_n, int(np.log10(min_n)), int(np.log10(max_n))


def map_to_single_value(i, j):
    if i < 6 and (j == 0 or j == 1):
        k = i * 4 + j
        print("i < 6: i = {}, j = {}, k = {}".format(i, j, k))
    elif i >= 6 and (j == 2 or j == 3):
        i -= 6
        k = i * 4 + j
        print("i >= 6: i = {}, j = {}, k = {}".format(i, j, k))
    else:
        raise ValueError("i must be between 0 and 3, and j must be between 0 and 11")
    return k


def configure_projection_plot(p, field, min_n, max_n, fontsize):
    p.set_axes_unit('pc')
    p.set_font_size(fontsize)
    # Ensure the colorbar limits match for all plots
    p.set_cmap(field, 'viridis')
    p.set_zlim(("gas", field), min_n, max_n)
    p.hide_colorbar()
    return p


def set_ticks_and_labels(grid, k, xticks, plot):
    # ticks + ticklabels
    grid[k].axes.set_xticks([])
    grid[k].axes.set_yticks([])
    grid[k].axes.minorticks_on()
    grid[k].axes.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    grid[k].axes.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    grid[k].axes.set_xticks(xticks, major=True, crs=plot)
    grid[k].axes.set_yticks(xticks, major=True, crs=plot)
    grid[k].axes.set_xticklabels([str(x) for x in xticks], rotation=90)
    grid[k].axes.set_yticklabels([str(x) for x in xticks])


def set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, max_n, fontsize, ticker):
    # y axis time labels
    if j == 0:
        grid[k].axes.set_ylabel("{:.2f} Myr".format((ss_age/1e6)[0]))

    # x axis space labels and colorbar
    if i == (len(dds2)*2-1) or i == (len(DS)-1):
        # ylabel
        grid[k].axes.set_xlabel("(pc)")

        # colorbar
        grid.cbar_axes[k].set_ylabel(r'Number Density \big($\rm \frac{1}{cm^{3}}$\big)')
        grid.cbar_axes[k].minorticks_on()

        # make minorticks
        a1 = np.logspace(min_n_index, max_n_index, np.abs(min_n_index - max_n_index)+1)
        a2 = np.arange(1, 10, 1)
        minorticks = np.outer(a1, a2).flatten()
        atol = 0.9*10**(max_n_index)
        end_i = first_index(minorticks, float(max_n.d), rtol=0.1, atol=atol)
        minorticks = minorticks[:end_i]
        grid.cbar_axes[k].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        grid.cbar_axes[k].set_yticks(minorticks, minor=True)
        grid.cbar_axes[k].tick_params(labelsize=fontsize)


def overlay_quadrants(grid, edgecolor='yellow', linewidth=5):
    # Get the exact coordinates for the lines
    right_edge_of_second_column = grid[1].get_xlim()[1]
    bottom_edge_of_third_row = grid[2 * 4].get_ylim()[0]

    # Draw vertical line
    for i in range(6):  # Six rows in the grid
        grid[i * 4 + 1].axvline(x=right_edge_of_second_column, color=edgecolor, linewidth=linewidth)

    # Draw horizontal line
    for i in range(4):  # Four columns in the grid
        grid[2 * 4 + i].axhline(y=bottom_edge_of_third_row, color=edgecolor, linewidth=linewidth)

    print("Quadrants overlayed")


def main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n_factor=2e5, max_n_factor=0.30):

    # set up figure and axes grid
    configure_plots()
    fig = plt.figure()
    grid = create_axes_grid(fig, dds1, sim)

    # get min and max values for colorbar
    min_n, max_n, min_n_index, max_n_index = get_min_max_values(root_dir, sim, dds2, min_n_factor, max_n_factor)

    # load datasets into list DS and labels into list LABEL
    DS, LABEL = load_datasets(root_dir, sim, dds1, dds2, dds3)

    # loop over datasets (12) and make 2 projection plots per dataset
    for i, ds in enumerate(DS):
        # ... get disk and L, ortho_find ...

        # grab bh properties and make sphere centered on BH
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = ss_pos
        r = 2000*yt.units.pc
        sp = ds.sphere(center, 2 * r)

        # make disk data container and define angular momentum vector L
        disc_r_pc = 2.1
        disc_h_pc = 2.1
        disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

        # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
        # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
        # It's very useful for giving you a vector edge-on to the disk.
        vecs = ortho_find(L)

        # s1 is the first 6 simulations (first 2 columns), s2 is the second 6 simulations (last 2 columns)
        if i < 6:
            js = [0, 1]
        else:
            js = [2, 3]

        # loop over the 2 columns of each seed
        for j in js:
            # set north vector and v
            north = vecs[1] if (j == 0) or (j == 2) else vecs[0]
            v = 0 if (j == 0) or (j == 2) else 1
            if north[-1] < 0:
                north *= [1,1,-1]
            
            # make projection plot and configure basic parameters
            p = make_projection_plot(ds, width_pc, disk, L, field, vecs, v, north, min_n, max_n, fontsize)
            p = configure_projection_plot(p, field, min_n, max_n, fontsize)

            # ... pre-render plot customization ...

            # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[("gas", field)]
            plot.figure = fig
            k = map_to_single_value(i, j)
            plot.axes = grid[k].axes

            # age, mass and simulation label in first and third columns
            if ((j == 0) or (j == 2)):
                p.annotate_text((0.07, 0.86), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), 
                    coord_system="axis", text_args={"color": "white", "fontsize": 8})
                
                # make colorbar
                plot.cax = grid.cbar_axes[k]

            # simulation label in first row and every second column only
            if ((i == 0) or (i == 6) or (i == 3) or (i == 9)) and ((j == 0) or (j == 2)):
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
            p.annotate_marker(center, coord_system="data", color="white")

            # render the plot
            p.render()

            # ... post render customization ...

            # Modify colorbar and axes properties **after** p.render() so that they are not overwritten.

            # ticks + ticklabels
            set_ticks_and_labels(grid, k, xticks, plot)

            # x and y extrenal axis labels and colorbar settngs
            set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, 
                                         max_n, fontsize, ticker)
            
    overlay_quadrants(grid)

    # ... save plot ...

    # Save plot
    num = str(sys.argv[-1])
    plot_name = f"multiplot_axesgrid_time_{width_pc.d}pc_{num}.pdf"
    plt.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/", plot_name)


if __name__ == "__main__":

    root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
                "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/"]
    sim = ["1B.RSm01", "1B.RSb01-2", "2B.RSm01", "2B.RSb01"]
    dds1 = ["DD0130/DD0130", "DD0133/DD0133", "DD0138/DD0138"] 
    dds2 = ["DD0201/DD0201", "DD0204/DD0204", "DD0209/DD0209"]  # 0.19, 0.49, 1.01 Myr for m01, 
    dds3 = ["DD0201/DD0201", "DD0203/DD0203", "DD0208/DD0208"]  # 0.201, 0.501, 1.002 Myr for b01
    fontsize = 14
    width_pc = 1.5 * yt.units.pc # must change xticks if this changes
    xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
    min_n_factor=2e5
    max_n_factor=0.30
    main(root_dir, sim, dds1, dds2, dds3, "number_density", width_pc, xticks, fontsize, min_n_factor, max_n_factor)