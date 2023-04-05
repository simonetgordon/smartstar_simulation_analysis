########################################  MultiPanel Projections  ########################################
#
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
from smartstar_find import ss_properties

from matplotlib.ticker import FuncFormatter
import yt.visualization.eps_writer as ytvis
import matplotlib as mpl
from matplotlib.colors import LogNorm


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


pyplot.rcParams['font.size'] = 12
pyplot.rcParams['font.weight'] = 'light'
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]


fig = plt.figure()

fontproperies = {'family': 'serif', 'color':  'black', 'weight': 'normal', 'size': 12}

grid = AxesGrid(
    fig,
    121,
    nrows_ncols=(3, 1),
    axes_pad=0.3,
    label_mode="L",
    aspect=False,
    #add_all=True,
    share_all=False,
    cbar_location="right",
    cbar_mode="single",
    cbar_size="3%",
    cbar_pad="0%",
)

proj = []
labels = []
DS = []
for i, dd in enumerate(dds):

    # data and data labels
    ds = yt.load(os.path.join(root_dir, sim[i], dd))
    label = tidy_data_labels(sim[i])
    DS.append(ds)
    labels.append(label)

    # grab bh properties and set projection dimensions
    field = "number_density"
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    widths = [150, 40, 10] * ds.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    fontsize = 10
    k = sys.argv[-1]
    width = widths[k]

    # projection plot
    p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center,
                          data_source=sp, weight_field='density')
    p.set_axes_unit('pc')

    # Ensure the colorbar limits match for all plots
    p.set_cmap(field, 'viridis')
    p.set_font_size(fontsize)
    p.set_zlim(("gas", field), 1e1, 8e9)
    p.hide_colorbar()

    # annotate projection
    if (i == 0) and (k == 2):
        p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
    if k == 0:
        # transaxes = grid[i].axes.transAxes
        # props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # grid[i].axes.text(0.05, 0.95, label, transform=transaxes, fontsize=14,
        #                   verticalalignment='top', bbox=props)
        p.annotate_text(
            [0.08, 0.08],
            str(label),
            coord_system="axis",
            text_args={"color": "black"},
            inset_box_args={
                "boxstyle": "square,pad=0.3",
                "facecolor": "white",
                "linewidth": 3,
                "edgecolor": "white",
                "alpha": 0.5,
            },
        )
    if i == 1:
        p.annotate_scale(corner='lower_right', coeff=1, unit='pccm')

    p.annotate_text((0.17, 0.90), "BH Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    if k == 2:
        p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
    # p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
    p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3))

    # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", field)]
    plot.figure = fig
    plot.axes = grid[i].axes
    # if j in [1, 2]:
    #     grid[i * 3 + j].set_aspect(2)
    if k == 2:
        plot.cax = grid.cbar_axes[i]

    # redraw the plot
    p.render()

    # Modify the axes properties **after** p.render() so that they
    # are not overwritten.
    if k == 2:
        ticklabels = grid.cbar_axes[i].get_yticklabels()
        grid.cbar_axes[i].set_yticklabels(ticklabels, fontsize=12)
        grid.cbar_axes[i].set_ylabel(r'Number Density \big($\rm cm^{-3}$\big)', fontsize=12)
        grid.cbar_axes[i].minorticks_on()
        # make minorticks
        a1 = np.logspace(1, 9, 9)
        a2 = np.arange(1, 10, 1)
        minorticks = np.outer(a1, a2).flatten()
        grid.cbar_axes[i].set_yticks(minorticks, minor=True)
        grid.cbar_axes[i].tick_params(labelsize=12)
    #plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))

plt.savefig(f"multiplot_axesgrid_{widths.d[k]}pccm.pdf")
