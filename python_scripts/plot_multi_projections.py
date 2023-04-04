from smartstar_find import ss_properties
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
import yt.visualization.eps_writer as ytvis
import matplotlib as mpl
from matplotlib import pyplot, colors
from matplotlib.ticker import FuncFormatter
import numpy as np
import os
from matplotlib.colors import LogNorm
import yt


pyplot.rcParams['font.size'] = 12

from grid_figure import GridFigure
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import G
from yt.visualization.color_maps import yt_colormaps


# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]

# n_col = n_row = 3
# fig, axs = plt.subplots(nrows=n_row, ncols=n_col, figsize=(9, 9))
fig = plt.figure()

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
    cbar_pad="1%",
)

proj = []
labels = []
DS = []
for i, dd in enumerate(dds):

    ds = yt.load(os.path.join(root_dir, sim[i], dd))
    label = "s1_" + str(float(ds.current_time.to('Myr')))[:5] + "_Myr"
    DS.append(ds)
    labels.append(label)

    field = "number_density"
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    widths = [150, 40, 20] * ds.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    fontsize = 10
    width = widths[0]

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
    if i == 0:
        p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
    if i == 2:
        p.annotate_scale(corner='lower_right', coeff=1, unit='pccm')

    p.annotate_text((0.17, 0.90), "BH Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
    # p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
    p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3))

    # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots[("gas", field)]
    plot.figure = fig
    plot.axes = grid[i].axes
    #plot.cax = grid.cbar_axes[i]
    # if j in [1, 2]:
    #     grid[i * 3 + j].set_aspect(2)
    if i == 0:
        plot.cax = grid.cbar_axes[i]

    # redraw the plot
    p.render()

    # Modify the axes properties **after** p.render() so that they
    # are not overwritten.
    if i == 0:
        ticklabels = grid.cbar_axes[i].get_yticklabels()
        grid.cbar_axes[i].set_yticklabels(ticklabels, fontsize=12)
        grid.cbar_axes[i].set_ylabel(r'$\rm Number \, Density \, \left(cm^{-3}\right)$', fontsize=12)
        grid.cbar_axes[i].minorticks_on()

    #plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))

plt.savefig("multiplot_axesgrid.pdf")
