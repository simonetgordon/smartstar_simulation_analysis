########################################  MultiPanel Projections  ########################################
# Uses a combination of AxesGrid and yt.ProjectionPlot to produce a multipanel figure.
# Can only produce this figure column by column. The full row x col figure needs to be assembled in ppt.
# to run: python plot_multipanel_projections.py [width index k]
##########################################################################################################

import os
import sys
import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import pyplot, colors
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib_scalebar.scalebar import ScaleBar
#from matplotlib_scalebar.dimension import Dimension
import matplotlib.ticker as ticker
from contextlib import redirect_stdout
from smartstar_find import ss_properties
import re # complex str searches


def tidy_data_labels(labels):
    # for lists of labels
    if len(labels) < 5:
        data_labels = [i.replace("-2", "") for i in labels]
        data_labels = [i.replace("RS", "") for i in data_labels]
    # for single label
    else:
        data_labels = labels.replace("-2", "")
        data_labels = data_labels.replace("RS", "")
    return data_labels


def first_index(a, val, rtol=0.1, atol=10):
    return next(j for j, _ in enumerate(a) if np.isclose(_, val, rtol, atol))


def format_sci_notation(x):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'$\rm {} \times 10^{{{}}}$'.format(a, b)


# input data - simulations and individual outputs
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]

# font settings
pyplot.rcParams['font.size'] = 12
pyplot.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)
plt.rcParams["mathtext.default"] = "regular"
fontproperties = {'family': 'serif', 'color':  'black', 'weight': 'normal', 'size': 12}
fontsize = 10 # for projection annotations

# make AxesGrid figure
fig = plt.figure()
grid = AxesGrid(
    fig,
    121,
    nrows_ncols=(3, 1),
    axes_pad=0.3,
    label_mode="L",
    aspect=False,
    share_all=False,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)

# find min and max field values from highest resolution simulation - for colorbar
field = "number_density"
ds_hr = yt.load(os.path.join(root_dir, sim[-1], dds[-1]))
min_n = ds_hr.r[("gas", field)].min()*100 # bring up to ~ 10^-1 order
max_n = ds_hr.r[("gas", field)].max()
max_n_index = int(np.log10(max_n))
min_n_index = int(np.log10(min_n))

# iterate over datadumps
for i, dd in enumerate(dds):

    # make dataset and simulation label
    ds = yt.load(os.path.join(root_dir, sim[i], dd))
    label = tidy_data_labels(sim[i])

    # grab bh properties and make sphere centered on BH
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    r = 2000*yt.units.pc
    sp = ds.sphere(center, 2 * r)

    # set projection width from CLI input
    widths = [4000, 80, 10] * ds.units.pccm
    k = int(sys.argv[-1])
    width = widths[k]

    # make projection plot
    p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center, data_source=sp, weight_field='density')
    p.set_axes_unit('pc')
    p.set_font_size(fontsize)

    # Ensure the colorbar limits match for all plots
    p.set_cmap(field, 'viridis')
    p.set_zlim(("gas", field), min_n, max_n)
    p.hide_colorbar()

    # annotate projection based on its [i, k] ([row, column]) value
    # if (i == 0) and (k == 2):
    #     p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)

    # find smallest cell width
    with open('ds_stats.txt', 'w') as f:
        with redirect_stdout(f):
            ds.print_stats()
    with open('ds_stats.txt', 'r') as f:
        for line in f:
            cell_width = re.search(r'Width: (.{9}) pc', line)
            if cell_width:
                dx = cell_width.group(1)
    f.close()

    if k == 0:
        # age, mass and simulation label in first
        p.annotate_title("dx = {} pc".format(format_sci_notation(float(dx))))
        p.annotate_text((0.24, 0.90), r"BH Mass: {} $\rm M_\odot$".format(int(ss_mass.d)), coord_system="axis",
                        text_args={"color": "white"})
        p.annotate_text([0.07, 0.08], str(label), coord_system="axis", text_args={"color": "black"},
                        inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                        "edgecolor": "white", "alpha": 0.5},
                        )

    if k == 1:
        p.annotate_title("BH Age = {:.2f} Myr".format(ss_age[0] / 1e6))

    # plot scalebars
    if i == 2:
        coeff = [10, 1e4, 1e4]
        unit = ['pc', 'au', 'au']
        p.annotate_scale(corner='lower_right', coeff=coeff[k], unit=unit[k], pos=(0.86, 0.08), max_frac=0.3,
                         scale_text_format="{scale:.2e} {units}")

    if i == 1:
        # Define a custom length unit based on meter conversion
        # pc = {'name': 'pc', 'scale': 3.086e+16}
        # my_dim = 'pixel-length'

        # Create a ScaleBar object using the custom unit
        #scalebar = ScaleBar(0.5, units='pc', dimension=my_dim, length_fraction=0.2, location='lower right')
        #scalebar.length_unit = {'pc': pc}

    # mark BH position
    if k == 2:
        p.annotate_marker(center, coord_system="data", color="white")

    # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", field)]
    plot.figure = fig
    plot.axes = grid[i].axes
    if k == 2: # temp 2 -> 0
        plot.cax = grid.cbar_axes[i]
        grid.cbar_axes[i].set_ymargin(-0.1)
        grid.cbar_axes[i].set_anchor((0.6, 0.2))

    # redraw the plot
    p.render()

    # Modify the scalebar, colorbar and axes properties **after** p.render() so that they are not overwritten.
    if (i == 2) and (k != 0): # i=1 -> i=2 temp
        plot.axes.get_yaxis().set_units('pc')
        #plot.axes.add_artist(scalebar)
    if k == 2:
        grid[i].axes.set_yticklabels("", fontsize=12)
        grid[i].axes.set_ylabel("")
    if i == 2: #temp 2 -> 0

        # ticklabels = grid.cbar_axes[i].get_yticklabels()
        # grid.cbar_axes[i].set_yticklabels(ticklabels, fontsize=12)
        grid.cbar_axes[i].set_ylabel(r'Number Density \big($\rm \frac{1}{cm^{3}}$\big)', fontsize=12)
        grid.cbar_axes[i].minorticks_on()

        # make minorticks
        a1 = np.logspace(min_n_index, max_n_index, np.abs(min_n_index - max_n_index)+1)
        a2 = np.arange(1, 10, 1)
        minorticks = np.outer(a1, a2).flatten()
        atol = 10**(max_n_index-1)
        end_i = first_index(minorticks, float(max_n.d), rtol=0.1, atol=atol)
        minorticks = minorticks[:end_i]
        grid.cbar_axes[i].xaxis.set_minor_locator(ticker.AutoMinorLocator())
        grid.cbar_axes[i].set_yticks(minorticks, minor=True)
        grid.cbar_axes[i].tick_params(labelsize=12)

        ax_pos = grid[i].axes.get_position()
        cb_pos = grid.cbar_axes[i].get_position()
        cb_pos = grid.cbar_axes[i].set_position([cb_pos.x0 - 0.5, ax_pos.y0, cb_pos.width, cb_pos.height])

plt.savefig(f"multiplot_axesgrid_{widths.d[k]}pccm.pdf")
