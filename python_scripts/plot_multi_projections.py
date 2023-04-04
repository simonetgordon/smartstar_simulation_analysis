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


pyplot.rcParams['font.size'] = 16

from grid_figure import GridFigure
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import G
from yt.visualization.color_maps import yt_colormaps


# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]

n_col = n_row = 3
my_fig = GridFigure(n_col, n_row, figsize=(9, 9),
                    left_buffer=0.15, right_buffer=0.2,
                    bottom_buffer=0.1, top_buffer=0.1,
                    vertical_buffer=0, horizontal_buffer=0.12)

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
    widths = [300, 40, 20] * ds.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    fontsize = 8
    for j, width in enumerate(widths):
        ax = my_fig[i*n_col + j]
        p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center,
                              data_source=sp, weight_field='density')
        p.set_axes_unit('pc')

        # Ensure the colorbar limits match for all plots
        p.set_cmap(field, 'viridis')
        p.set_font_size(fontsize)
        p.set_zlim(("gas", field), 1e1, 8e9)
        p.hide_colorbar()

        # annotate projection
        p.annotate_scale(corner='lower_left', coeff=10000, unit='au')
        p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
        p.annotate_text((0.62, 0.92), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                        text_args={"color": "white"})
        p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
        # p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
        p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3))

        plot = ax.imshow(p.frb['density'].v, norm=LogNorm())

plt.savefig("multiplot_gridfigure.pdf")
