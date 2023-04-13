##################################  MultiPanel Projections - over time ###################################
# Uses a combination of AxesGrid and yt.ProjectionPlot to produce a multipanel figure accross time.
# The size of each panel is fixed in spatial dimension.
# to run: python plot_multipanel_projections.py 0
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
from smartstar_find import ss_properties
import re # complex str searches
from plot_multi_projections import tidy_data_labels, first_index, format_sci_notation
from plot_disc_projections import _make_disk_L

# input data - simulations and individual outputs
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSb01-2", "1B.RSm01"]
dds = ["DD0128/DD0128", "DD0130/DD0130", 
       "DD0132/DD0132", 
       "DD0134/DD0134", 
       "DD0136/DD0136", 
       "DD0138/DD0138",
       #"DD0132/DD0132", "DD0133/DD0133", 
       #"DD0134/DD0134", "DD0135/DD0135", "DD0136/DD0136", "DD0137/DD0137", "DD0138/DD0138"
       ]

# font settings
pyplot.rcParams['font.size'] = 14
pyplot.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)
plt.rcParams["mathtext.default"] = "regular"
fontsize = 14 # for projection annotations

# make AxesGrid figure
fig = plt.figure()
grid = AxesGrid(
    fig,
    (0.01, 0.01, 0.75, 1.15),
    nrows_ncols=(len(dds), len(sim)*2),
    axes_pad=0,
    label_mode="L",
    aspect=False,
    share_all=True,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3.2%",
    cbar_pad="0.2%",
)

# find min and max field values from last simulation - for colorbar
field = "number_density"
ds_hr = yt.load(os.path.join(root_dir, sim[0], dds[-1]))
min_n = ds_hr.r[("gas", field)].min()*1e4 # bring up to ~ 10^3 order
max_n = ds_hr.r[("gas", field)].max()
max_n_index = int(np.log10(max_n))
min_n_index = int(np.log10(min_n))

# iterate over datadumps
for s,_ in enumerate(sim):
    for i, dd in enumerate(dds):

        # make dataset and simulation label
        ds = yt.load(os.path.join(root_dir, sim[s], dd))
        label = tidy_data_labels(sim[s])

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
        v = 0

        # set projection width in pc
        width = 1.5 * ds.units.pc

        if s == 0:
            js = [0, 1]
        else:
            js = [2, 3]

        for j in js:
            # make projection plot
            if (j == 0) or (j == 2):
                north = vecs[1]
                p = yt.ProjectionPlot(ds, vecs[v], ("gas", field), weight_field=("gas", "density"), north_vector=north, 
                            center=disk.center, width=width, data_source=disk)
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

            # mark BH position
            p.annotate_marker(center, coord_system="data", color="white")

            # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
            plot = p.plots[("gas", field)]
            plot.figure = fig
            k = 4*i + j # grid index
            plot.axes = grid[k].axes

            # first column only
            if (j == 0) or (j == 2):
                # age, mass and simulation label in first
                #p.set_font_size(10)
                p.annotate_text((0.07, 0.86), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), 
                                coord_system="axis", text_args={"color": "white", "fontsize": 8})
                #p.set_font_size(fontsize)
                
                # make colorbar
                plot.cax = grid.cbar_axes[k]
                #grid.cbar_axes[k].set_ymargin(-0.1)
                #grid.cbar_axes[k].set_anchor((0.6, 0.2))

            p.render()

            # Modify the scalebar, colorbar and axes properties **after** p.render() so that they are not overwritten.

            # first row and every second column only
            if (i == 0) and ((j == 0) or (j == 2)):
                #p.annotate_title(str(label), font=16)
                grid[k].axes.set_title(str(label), fontsize=18)

            # ticks + ticklabels
            grid[k].axes.set_xticks([])
            grid[k].axes.set_yticks([])
            grid[k].axes.minorticks_on()
            xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
            grid[k].axes.xaxis.set_minor_locator(ticker.AutoMinorLocator())
            grid[k].axes.yaxis.set_minor_locator(ticker.AutoMinorLocator())
            grid[k].axes.set_xticks(xticks, major=True, crs=plot)
            grid[k].axes.set_yticks(xticks, major=True, crs=plot)
            grid[k].axes.set_xticklabels([str(x) for x in xticks], rotation=90)
            grid[k].axes.set_yticklabels([str(x) for x in xticks])

            if j == 0:
                grid[k].axes.set_ylabel("{:.2f}".format(ds.current_time.to('Myr')))

            if i == (len(dds)-1): #temp 2 -> 0

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
                #cb_pos = grid.cbar_axes[k].set_position([cb_pos.x0 - 0.1, ax_pos.y0, cb_pos.width, cb_pos.height])

plot_name = f"multiplot_axesgrid_time_{width.d}pc.pdf"    
plt.savefig('plots/'+ plot_name, bbox_inches='tight')
print("created plots/", plot_name)
