import yt
import sys
import os
from derived_fields import add_fields_ds
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
from scipy.interpolate import interp1d
import sys
import yaml
import yt
import ytree
from ytree.data_structures.tree_container import TreeContainer
from yt.extensions.p2p.tree_analysis_operations import get_progenitor_line

pyplot.rcParams['font.size'] = 16

from grid_figure import GridFigure
from yt.extensions.p2p.stars import get_star_data
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import G
from yt.visualization.color_maps import yt_colormaps


# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]

my_fig = GridFigure(3, 3, figsize=(9, 9),
                    left_buffer=0.15, right_buffer=0.2,
                    bottom_buffer=0.1, top_buffer=0.1,
                    vertical_buffer=0, horizontal_buffer=0.12)


plt.show()
# grid = AxesGrid(
#     fig,
#     121,
#     nrows_ncols=(3, 3),
#     axes_pad=0.3,
#     label_mode="L",
#     #add_all=True,
#     share_all=False,
#     cbar_location="right",
#     cbar_mode="single",
#     cbar_size="3%",
#     cbar_pad="0%",
# )

proj = []
labels = []
DS = []
for i, dd in enumerate(dds):

    ds = yt.load(os.path.join(root_dir, sim[i], dd))
    #add_fields_ds(ds)
    label = "s1_" + str(float(ds.current_time.to('Myr')))[:5] + "_Myr"
    DS.append(ds)
    labels.append(label)

    field = "number_density"
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    widths = [300, 40, 20] * ds.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    fontsize = 8
    for j, width in enumerate(widths):
        p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center,
                              data_source=sp, weight_field='density')
        p.set_axes_unit('pc')

        # Ensure the colorbar limits match for all plots
        p.set_cmap(field, 'viridis')
        p.set_font_size(fontsize)
        p.set_zlim(("gas", field), 1e1, 8e9)
        # if j == 1:
        #     plt.xlim(-1, 1)
        #     plt.ylim(-1, 1)
        # if j == 2:
        #     plt.xlim(-7.5e4, 7.5e4)
        #     plt.ylim(-7.5e4, 7.5e4)g
        p.hide_colorbar()

        # annotate
        p.annotate_scale(corner='lower_left', coeff=10000, unit='au')
        p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
        p.annotate_text((0.62, 0.92), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                        text_args={"color": "white"})
        p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
        # p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
        p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3))

        # this forces the ProjectionPlot to redraw itself on the AxesGrid axes.
        plot = p.plots[("gas", field)]
        plot.figure = fig
        plot.axes = grid[i*3 + j].axes
        plot.cax = grid.cbar_axes[i*3 + j]
        if j in [1, 2]:
            grid[i*3 + j].set_aspect(1)
        # if i == 0:
        #     plot.cax = grid.cbar_axes[i*3 + j]

        # redraw the plot
        p.render()

        # Modify the axes properties **after** p.render() so that they
        # are not overwritten.
        plot.axes.xaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))


plt.savefig("multiplot_axesgrid.pdf")

# mp = ytvis.multiplot(
#     3,
#     3,
#     proj,
#     # xlabels=[False, False, False, True, True, True, False, False, False],
#     # ylabels=[True, False, True, False, True, False, False, True, False],
#     savefig="yt",
#     shrink_cb=0.99,
#     bare_axes=False,
#     yt_nocbar=False,
#     margins=(0.5, 0.5),
#     cb_flags=[False, False, False, False, False, False, False, True, False], # middle, bottom, top
#     xaxis_flags=[False, False, False, True, True, True, False, False, False,
#                  False, False, False, True, True, True, False, False, False],
#     yaxis_flags=[True, False, True, False, True, False, False, True, False],
#)
#mp.scale_line(label="$r_{vir}$", labelloc="top")
#mp.save_fig("multiplot")
