"""
Create inset zoomed projection to illustrate Enzo as as a zoom-in simulation Call like:
python -i plot_zoomed_inset_projection.py
"""

import sys
import yt
import os
yt.enable_parallelism()
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find
from plot_multi_projections import tidy_data_labels
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from grid_figure import GridFigure
from matplotlib import pyplot
from matplotlib.colors import LogNorm

# input
root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16-2"
sims = ["DD0160/DD0160"]
ds = yt.load(os.path.join(root_dir, sims[0]))

sim_str = tidy_data_labels(str(root_dir[77:]))

## find north vector a (near-end) ds
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

# region larger than disc
my_width = (75*yt.units.pc).to("unitary")
region = ds.box(ss_pos - 0.5 * my_width, ss_pos + 0.5 * my_width)
width = (my_width, "unitary")

# set image orientation + north vector
orient = "face-on"

if orient == "face-on":
    orient_str = "face_on_"
    dir = vecs[0]
    north = vecs[1]
else:
    orient_str = "edge_on_"
    dir = vecs[1]
    north = vecs[0]

# initialise plot
pyplot.rcParams['font.size'] = 10
n_col = n_row = 1
my_fig = GridFigure(n_col, n_row, figsize=(9, 9),
                    left_buffer=0.15, right_buffer=0.2,
                    bottom_buffer=0.1, top_buffer=0.1,
                    vertical_buffer=0, horizontal_buffer=0.01)

# h number density
field = "number_density"
ax = my_fig[0]
p1 = yt.ProjectionPlot(ds, dir, ("gas", field), width=(15, 'pc'), north_vector=north, center=center, data_source=region,
                        weight_field=field)
plot = p1.plots[("gas", field)]
p = ax.imshow(p1.frb['number_density'].v, norm=LogNorm())
plot.figure = ax.axes


# p.set_cmap(("gas", field), 'viridis')
# p.set_zlim(("gas", field), 7e3, 1e11)

# inset zoomed proj
pzoom = yt.ProjectionPlot(ds, dir, ("gas", field), width=(1.5, 'pc'), north_vector=north, center=center, 
                          data_source=region, weight_field=field)
axins = zoomed_inset_axes(p, 10, loc=1) # zoom = 10
axins.imshow(pzoom.frb['number_density'].v, norm=LogNorm())

# axins.imshow(p1, interpolation="nearest", origin="lower")
# sub region of the original image
x1, x2, y1, y2 = -0.75, 0.75, -0.75, 0.75
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)

# plt.xticks(visible=False)
# plt.yticks(visible=False)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")


# save
dirname = "plots/inset_plot.png"
p1.save(dirname)

