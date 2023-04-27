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
from grid_figure import GridFigure
from yt.utilities.cosmology import Cosmology
from yt.utilities.physical_constants import G
from yt.visualization.color_maps import yt_colormaps

pyplot.rcParams['font.size'] = 10
n_col = n_row = 2
my_fig = GridFigure(n_col, n_row, figsize=(9, 9),
                    left_buffer=0.15, right_buffer=0.2,
                    bottom_buffer=0.1, top_buffer=0.1,
                    vertical_buffer=0, horizontal_buffer=0.01)

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0148/DD0148", "DD0238/DD0238"]


proj = []
labels = []
DS = []
field = "number_density"
# get max and min density values for cmap
ds_hr = yt.load(os.path.join(root_dir, sim[-1], dds[-1]))
n_max = ds_hr.r[("gas", field)].max().d
n_min = ds_hr.r[("gas", field)].min().d

for i, dd in enumerate(dds[:-1]):
    ds = yt.load(os.path.join(root_dir, sim[i], dd))
    label = "s1_" + str(float(ds.current_time.to('Myr')))[:5] + "_Myr"
    DS.append(ds)
    labels.append(label)

    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    widths = [300, 40] * ds.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    fontsize = 8
    for j, width in enumerate(widths):
        ax = my_fig[i*n_col + j]
        p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center, data_source=sp, weight_field='density')
        p = ax.imshow(p.frb['density'].v, norm=LogNorm())

        # Ensure the colorbar limits match for all plots
        p.set_cmap('viridis')

        # annotate projection
        # ax.annotate_scale(corner='lower_left', coeff=10000, unit='au')
        # ax.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
        # ax.annotate_text((0.62, 0.92), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis", text_args={"color": "white"})
        # ax.annotate_marker(center, coord_system="data", color="white")  # mark ss position

        # plt annotations
        ax.plot(center, 'ro') # bh position
        ax.annotate("BH Mass: {:.2f} Msun".format(ss_mass.d), xy=(0.62, 0.92), xycoords='axes fraction', fontsize=8, color='w')
        ax.annotate("Time: {:.2f} Myr".format(ds.current_time.to('Myr').d), xy=(0.32, 0.12), xycoords='axes fraction', fontsize=8, color='w')
        ax.set_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3), fontsize=10)


for my_axes in my_fig.left_axes:
    my_axes.tick_params(axis="y", left=True, direction="inout", which="both")
    my_axes.tick_params(axis="y", right=True, direction="in", which="both")

# axes labels 1st column
for i in [0, 2]:
    ds = yt.load(os.path.join(root_dir, sim[0], dds[0]))
    w = ds.quan(widths[0]).to('pc')
    my_fig[i].yaxis.set_label_text("$\\rm y [pc]$")
    my_fig[i].yaxis.set_ticks(np.linspace(-w, w, 6))
    #my_fig[i].yaxis.set_ticks(np.logspace(-22, -20, 2), minor=True, labels="")


# colorbar
norm = colors.Normalize(vmin=np.log(n_min), vmax=np.log(n_max))
n_vals = np.logspace(np.log(n_min), np.log(n_max), 5)
cmap = "viridis"
my_cax = my_fig.add_cax(my_fig[1], "right", length=1)
cbar = mpl.colorbar.ColorbarBase(
    my_cax, cmap=cmap, boundaries=n_vals,
    norm=norm, orientation='vertical')

cbar.set_label("Number Density [$\\rm cm^{-3}$]")
#my_ticks = list(np.linspace(-5, -3, 9))
# my_ticks = list(np.linspace(-5, -3, 5))
#cbar.set_ticks()
cbar.solids.set_edgecolor("face")



plt.savefig("multiplot_gridfigure.pdf")