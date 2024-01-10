##################################  MultiPanel Projections - over time ###################################
# Uses a combination of AxesGrid and yt.ProjectionPlot to produce a multipanel at a fixed time.
# 2x2 plot
# The size of each panel is fixed in spatial dimension.
# to run: python -i plot_disc_projections_fixed_t.py "fixed_t"
##########################################################################################################

import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import pyplot, colors
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.utilities.math_utils import ortho_find
import matplotlib.ticker as ticker
from contextlib import redirect_stdout
import re # complex str searches
from helper_functions import tidy_data_labels, first_index, format_sci_notation, _make_disk_L, ss_properties

# input data - simulations and individual outputs
root_dir = ["/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/",
            "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/"]
sim = ["1B.RSb01-2", "1S.RSb01"]
dds = ["DD0131/DD0131", "DD0133/DD0133"]

# font settings
pyplot.rcParams['font.size'] = 16
pyplot.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)
plt.rcParams["mathtext.default"] = "regular"
fontsize = 16 # for projection annotations

# make AxesGrid figure
fig = plt.figure()
grid = AxesGrid(
    fig,
    (0.01, 0.01, 0.695, 0.70),
    nrows_ncols=(len(dds), len(sim)),
    axes_pad=0,
    label_mode="L",
    aspect=False,
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)

# find min and max field values most refined simulation - for colorbar
field = "number_density"
ds_hr = yt.load(os.path.join(root_dir[-1], sim[-1], dds[-1]))
min_n = ds_hr.r[("gas", field)].min()*1e4 # bring up to ~ 10^2 order
#max_n = ds_hr.r[("gas", field)].max()*1e-3 # bring down to ~ 10^8 order
max_n = 1e8/(yt.units.cm)**3
max_n_index = int(np.log10(max_n))
min_n_index = int(np.log10(min_n))

# iterate over datadumps
for s,_ in enumerate(sim):
    i = s
    # make dataset and simulation label
    ds = yt.load(os.path.join(root_dir[s], sim[s], dds[s]))
    label = tidy_data_labels(sim[s])

    # grab bh properties and make sphere centered on BH
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    r = 2000*yt.units.pc
    sp = ds.sphere(center, 2 * r)

    # make disk data container and define angular momentum vector L
    disc_r_pc = 6.1
    disc_h_pc = 4.1
    disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    vecs = ortho_find(L)
    v = 0

    # set projection width in pc
    width = 8 * ds.units.pc

    # if s == 0:
    #     js = [0, 1]
    # else:
    #     js = [2, 3]

    js = [0, 1]
    for j in js:
        # face-on projection
        if j == 0:
            north = vecs[1]
            p = yt.ProjectionPlot(ds, vecs[v], ("gas", field), weight_field=("gas", "density"), north_vector=north, 
                        center=disk.center, width=width, data_source=disk)
        # edge-on projection
        else:
            north = vecs[0]
            p = yt.ProjectionPlot(ds, vecs[v+2], ("gas", field), weight_field=("gas", "density"), north_vector=north, 
                                center=disk.center, width=2 * disk.radius, data_source=disk)
        p.set_axes_unit('pc')
        p.set_font_size(fontsize)

        # Ensure the colorbar limits match for all plots
        p.set_cmap(field, 'viridis')
        p.set_zlim(("gas", field), min_n, max_n)
        p.hide_colorbar()

        # Annotations

        # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[("gas", field)]
        plot.figure = fig
        k = 2*i + j # grid index
        plot.axes = grid[k].axes

        # first column only
        if j == 0:
            # age, mass and simulation label in first
            #p.set_font_size(10)
            p.annotate_text((0.06, 0.89), r"BH Mass: {} $\rm M_\odot$, Age: {:.2f} Myr".format(int(ss_mass.d), ss_age[0] / 1e6), 
                            coord_system="axis", text_args={"color": "white"})
            
            # make colorbar
            plot.cax = grid.cbar_axes[k]
        else:
            # mark BH position
            p.annotate_marker(center, coord_system="data", color="green")

        p.render()

        # Modify the scalebar, colorbar and axes properties **after** p.render() so that they are not overwritten.

        # ticks + ticklabels
        xticks = grid[k].axes.get_xticks()[1:-1]
        grid[k].axes.set_xticks([])
        grid[k].axes.set_yticks([])
        grid[k].axes.minorticks_on()
        #xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
        grid[k].axes.xaxis.set_minor_locator(ticker.AutoMinorLocator())
        grid[k].axes.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        grid[k].axes.set_xticks(xticks, major=True, crs=plot, rotation=90, fontsize=14)
        grid[k].axes.set_yticks(xticks, major=True, crs=plot, fontsize=14)
        #grid[k].axes.set_xticklabels([str(x) for x in xticks], rotation=90)
        #grid[k].axes.set_yticklabels([str(x) for x in xticks])

        if j == 0:
            grid[k].axes.set_ylabel(str(label), fontsize=fontsize+1)

        # bottom row
        if i == (len(dds)-1):

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

            ax_pos = grid[k].axes.get_position()
            cb_pos = grid.cbar_axes[k].get_position()

num = str(sys.argv[-1])
plot_name = f"multiplot_axesgrid_time_{width.d}pc_{num}.pdf"    
plt.savefig('plots/'+ plot_name, bbox_inches='tight')
print("created plots/", plot_name)
