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


def load_datasets(root_dir, sim, dds1, dds2, dds3=None):
    DS = []
    LABEL = []
    for l in range(len(sim)):
        for s, _ in enumerate(dds1):
            if "2B.RSb" in sim[l]:
                dd = dds3[s]
            elif ("2B.RSm" in sim[l]) or ("no-SN" in sim[l]):
                dd = dds2[s]
                print("dd2")
                print(root_dir[l], sim[l], dd)
            else:
                dd = dds1[s]
                print("dd1")
                print(root_dir[l], sim[l], dd)
            ds = yt.load(os.path.join(root_dir[l], sim[l], dd))
            label = tidy_data_labels(sim[l])
            DS.append(ds)
            LABEL.append(label)
    return DS, LABEL


def configure_font(fontsize=14):
    pyplot.rcParams['font.size'] = fontsize
    pyplot.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"


def make_projection_plot(ds, width_pc, disk, L, field, min_n, max_n, vecs=None, v=None, north=None, dir="z", fontsize=14, cmap="viridis", center=None):
    #print("vecs: {}, v: {}, north: {}, disk.center: {}, width_pc: {}".format(vecs, v, north, disk.center.d, width_pc))
    center = center if center is not None else disk.center
    if vecs:
        p = yt.ProjectionPlot(ds, vecs[v], ("gas", field), weight_field=("gas", "density"), north_vector=north, center=center, width=width_pc, data_source=disk)
    else:
        p = yt.ProjectionPlot(ds, dir, ("gas", field), weight_field=("gas", "density"), center=center, width=width_pc, data_source=disk)
    p.set_axes_unit('pc')
    p.set_font_size(fontsize)
    p.set_cmap(field, cmap)
    p.set_zlim(("gas", field), min_n, max_n)
    p.hide_colorbar()
    return p


def create_axes_grid(fig, nrows, ncols, dim=(0.01, 0.01, 0.76, 1.152), axes_pad=0, share_all=True, cbar_location="right", cbar_size="2%", cbar_pad="0%"):
    return AxesGrid(
        fig,
        dim,
        nrows_ncols=(nrows,ncols),
        axes_pad=axes_pad,
        label_mode="L",
        aspect=False,
        share_all=share_all,
        cbar_location=cbar_location,
        cbar_mode="single",
        cbar_size=cbar_size,
        cbar_pad=cbar_pad,
    )


def get_min_max_values(root_dir, sim, dds2, field="number_density", min_n_factor=2e5, max_n_factor=0.30):
    ds_hr = yt.load(os.path.join(root_dir[-1], sim[-1], dds2[-1]))

    if field == "number_density":
        min_n = ds_hr.r[("gas", field)].min() * min_n_factor
        max_n = ds_hr.r[("gas", field)].max() * max_n_factor
    elif field == "temperature":
        min_n = ds_hr.r[("gas", field)].min() * min_n_factor
        max_n = ds_hr.r[("gas", field)].max() * max_n_factor
    return min_n, max_n, int(np.log10(min_n)), int(np.log10(max_n))


def map_to_single_value(i, j, i_lim=6):
    i_lim = 6
    if i < i_lim and (j == 0 or j == 1):
        k = i * 4 + j
        print("i < 6/3: i = {}, j = {}, k = {}".format(i, j, k))
    elif i >= i_lim and (j == 2 or j == 3):
        i -= i_lim
        k = i * 4 + j
        print("i >= 6/3: i = {}, j = {}, k = {}".format(i, j, k))
    else:
        raise ValueError(f"j ={j} must be between 0 and 3, and i={i} must be between 0 and 11")
    return k

# def map_to_single_value(i, j, num_rows_per_set=6, num_columns=4):
#     if 0 <= i < 12 and 0 <= j < 4:
#         # Adjust 'i' for the second set of rows
#         if i >= num_rows_per_set:
#             i -= num_rows_per_set
#             # Offset for the second half of the grid
#             offset = num_columns * num_rows_per_set
#         else:
#             offset = 0

#         # Column-major order mapping
#         k = j * num_rows_per_set + i + offset
#         print("i = {}, j = {}, k = {}".format(i, j, k))
#         return k
#     else:
#         raise ValueError(f"j ={j} must be between 0 and 3, and i={i} must be between 0 and 11")


def configure_projection_plot(p, field, cmap, min_n, max_n, fontsize):
    #p.set_axes_unit('pc')
    #p.set_font_size(fontsize)
    # Ensure the colorbar limits match for all plots
    print("min_n: {}, max_n: {}".format(min_n, max_n))
    print("field: {}".format(field))
    p.set_cmap(field, cmap)
    p.set_zlim(("gas", field), min_n, max_n)
    p.hide_colorbar()
    return p


def set_ticks_and_labels(grid, k, xticks, plot, no_xticklabels=True):
    # ticks + ticklabels
    grid[k].axes.set_xticks([])
    grid[k].axes.set_yticks([])
    grid[k].axes.minorticks_on()
    grid[k].axes.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    grid[k].axes.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    grid[k].axes.set_xticks(xticks, major=True, crs=plot)
    grid[k].axes.set_yticks(xticks, major=True, crs=plot)
    if no_xticklabels:
        grid[k].axes.set_xticklabels([])
        grid[k].axes.set_yticklabels([])
        return 0
    grid[k].axes.set_xticklabels([str(x) for x in xticks], rotation=90)
    grid[k].axes.set_yticklabels([str(x) for x in xticks])


def set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, max_n, fontsize, ticker):
    # y axis time labels
    if j == 0:
        grid[k].axes.set_ylabel("(pc)")
    # if (i == 5) or (i == 11):
    #     grid[k].axes.set_ylabel("(pc)")

    # x axis space labels and colorbar
    if (i == 5) or (i == 11): # i == (len(dds2)*2-1) or i == (len(DS)-1)
        print("setting xlabel and colorbar")
        # ylabel
        grid[k].axes.set_xlabel("(pc)")

        # colorbar
        k = i
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


def overlay_quadrants(sims, grid, edgecolor='black', linewidth=5, nrows=6, ncols=4, i_lim=6):
    # Get the exact coordinates for the lines
    right_edge_of_second_column = grid[1].get_xlim()[1]
    bottom_edge_of_third_row = grid[2 * 4].get_ylim()[0]

    # Draw vertical line
    for i in range(nrows):  # Six rows in the grid
        grid[i * ncols + 1].axvline(x=right_edge_of_second_column, color=edgecolor, linewidth=linewidth)

    # Draw horizontal line
    if sims == "baselines_1":
        for i in range(ncols):  # Four columns in the grid
            grid[2 * ncols + i].axhline(y=bottom_edge_of_third_row, color=edgecolor, linewidth=linewidth)

    print("Quadrants overlayed")


def main(sims, root_dir, nrows, ncols, dim, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize=14, min_n_factor=2e5, max_n_factor=0.30):

    # set up figure and axes grid
    configure_font()
    fig = plt.figure()
    grid = create_axes_grid(fig, nrows=nrows, ncols=ncols, dim=dim)

    # get min and max values for colorbar
    min_n, max_n, min_n_index, max_n_index = get_min_max_values(root_dir, sim, dds2, field=field, min_n_factor=min_n_factor, max_n_factor=max_n_factor)

    # load datasets into list DS and labels into list LABEL
    DS, LABEL = load_datasets(root_dir, sim, dds1, dds2, dds3)

    # loop over datasets (12) and make 2 projection plots per dataset
    for i, ds in enumerate(DS):
        # ... get disk and L, ortho_find ...

        # grab bh properties and make sphere centered on BH
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = ss_pos

        # make disk data container and define angular momentum vector L and north vector
        disc_r_pc = disc_h_pc = 0.15
        disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

        # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
        # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
        # It's very useful for giving you a vector edge-on to the disk.
        vecs = ortho_find(L)

        # s1 is the first 6 simulations (first 2 columns), s2 is the second 6 simulations (last 2 columns)
        js = [0, 1] if (i < 6) else [2, 3] # changing this 6 -> 3

        # loop over the 2 columns of each seed
        for j in js:
            
            # make projection plot and configure basic parameters. remake bigger disk for projection
            disc_r_pc = disc_h_pc = 2.1
            disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
            cmap = "RED TEMPERATURE" if field == "temperature" else "viridis"
            if sims == "baselines_1":
                # set north vector and v
                north = vecs[1] if (j == 0) or (j == 2) else vecs[0]
                v = 0 if (j == 0) or (j == 2) else 1
                # if north[-1] < 0:
                #     north *= [1,1,-1]
                p = make_projection_plot(ds, width_pc, disk, L, field, min_n, max_n, vecs=vecs, v=v, north=north, fontsize=fontsize,)

            elif sims == "baselines_2":
                dir = "z" if (j == 0) or (j == 2) else "y"
                p = make_projection_plot(ds, width_pc, disk, L, field, min_n, max_n, dir=dir, fontsize=fontsize, cmap=cmap)
            p = configure_projection_plot(p, field, cmap, min_n, max_n, fontsize)

            # ... pre-render plot customization ...

            # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[("gas", field)]
            plot.figure = fig
            if sims == "baselines_1":
                k = map_to_single_value(i, j)
            elif sims == "baselines_2":
                k = map_to_single_value(i, j, i_lim=3)
            plot.axes = grid[k].axes

            # age, mass and simulation label in first and third columns
            if ((j == 1) or (j == 3)):
                p.annotate_text((0.07, 0.86), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), 
                    coord_system="axis", text_args={"color": "white", "fontsize": 8})
                
                # make colorbar
                plot.cax = grid.cbar_axes[k]
            
            if ((j == 0) or (j == 2)):
                p.annotate_text((0.07, 0.86), "{:.2f} Myr".format((ss_age/1e6)[0]), 
                    coord_system="axis", text_args={"color": "white", "fontsize": 8})

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
            m_color = "black" if sims == "baselines_2" else "black"
            p.annotate_marker(center, coord_system="data", color=m_color, marker="x", plot_args={"s": 15, "alpha": 0.8,'linewidths': 0.2, 'edgecolors': m_color})

            plot.cax = grid.cbar_axes[i]
            # render the plot
            p.render()

            # ... post render customization ...

            # Modify colorbar and axes properties **after** p.render() so that they are not overwritten.

            # ticks + ticklabels
            set_ticks_and_labels(grid, k, xticks, plot, no_xticklabels=False)

            # x and y extrenal axis labels and colorbar settngs
            set_axis_labels_and_colorbar(grid, k, j, i, dds2, DS, ss_age, min_n_index, max_n_index, 
                                         max_n, fontsize, ticker)
            
            # set tick labels on bottom row
            # if ((k == 5) or (k == 11) or (k == 17) or (k == 23)):
            #     grid[k].axes.set_xticklabels([str(x) for x in xticks], rotation=90)
            
    overlay_quadrants(sims, grid, nrows=nrows, ncols=ncols) # skipping this for now

    # ... save plot ...

    # Save plot
    num = str(sys.argv[-1])
    plot_name = f"multiplot_axesgrid_time_{width_pc.d}pc_{num}_{max_n_factor}_factor.pdf"
    plt.savefig('plots/' + plot_name, bbox_inches='tight')
    print("created plots/", plot_name)


if __name__ == "__main__":

    sims = "baselines_1" # or baselines_2

    if sims == "baselines_1":
        # plot 6x4 multipanael of 1B and 2B baselines at 3 times
        # to run: python -i plot_multipanel_time_2.py "1B-2B-baselines-time-fix"
        root_dir = ["/Backup00/sgordon/disk14/cirrus-runs-rsync/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                    "/Backup00/sgordon/disk14/cirrus-runs-rsync/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
                    "/Backup00/sgordon/disk14/cirrus-runs-rsync/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/",
                    "/Backup00/sgordon/disk14/cirrus-runs-rsync/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/"]
        sim = ["1B.RSm01", "1B.RSb01-2", "2B.RSm01", "2B.RSb01"]
        dds1 = ["DD0131/DD0131", "DD0136/DD0136", "DD0138/DD0138"] 
        dds2 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.29, 079, 1.01 Myr for m01, 
        dds3 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.201, 0.801, 1.002 Myr for b01
        fontsize = 14
        width_pc = 1.5 * yt.units.pc # must change xticks if this changes
        xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        min_n_factor=2e5
        max_n_factor=0.2
        nrows=6
        ncols=4
        dim=(0.01, 0.01, 0.76, 1.152)
        field = "number_density"

    elif sims == "baselines_2":
        # plot 3x4 multipanel of 1S.m01 and 1S.m01-no-SN baselines at 3 times
        # to run: python plot_multipanel_time_2.py "1S.m01+1S.m01-no-SN"
        root_dir = ["/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/",
                    "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/",]
        sim = ["1S.RSm01","1S.m01-no-SN",]
        dds1 = ["DD0131/DD0131", "DD0132/DD0132", "DD0136/DD0136"]  # 0.1, 0.2, 0.6 (reaches 29 msun, final mass 30msun)
        dds2 = ["DD0152/DD0152", "DD0162/DD0162", "DD0202/DD0202"]  # 0.1, 0.2, 0.6 (reaches final mass of 60 msun)
        fontsize = 14
        width_pc = 1.5 * yt.units.pc # must change xticks if this changes
        xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        min_n_factor=20 # 2e5 for density baseline_1/2
        max_n_factor=0.5 # 0.001 for density baseline_2 
        nrows=3
        ncols=4
        dim=(0.01, 0.01, 0.73, 0.56)
        dds3=None
        field = "temperature"


    main(sims,root_dir, nrows, ncols, dim, sim, dds1, dds2, dds3,field, width_pc, xticks, fontsize, min_n_factor, max_n_factor)