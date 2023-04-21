"""
Plot projection of simulation in density or temperature. Call like:
python -i projection_plot.py DD0133/DD0133
"""

import yt
import shutil
import sys
import os
from smartstar_find import ss_properties
from find_disc_attributes import _make_disk_L
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time

# set by user
w_pccm = 1.5
field = "density"

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSb01"

# ## First check there is a local /data area                                                                                                                                                                 
# if not os.path.isdir("/data"):
#     sys.exit("Error: no /data")

# ## Second, check we have a directory. If not, then create one.                                                                                                                                             
# UserName=os.getlogin()
# LocalDir=os.path.join("/data",UserName)
# if not os.path.isdir(LocalDir):
#     print("Creating Directory "+LocalDir)
#     os.mkdir(LocalDir)

# ## Third, check if the data is already there, and if not, copy it over.                                                                                                                                    
# DataDumpFull = sys.argv[1]
# DataDump = DataDumpFull.split('/')
# LocalData = os.path.join(LocalDir,DataDump[0])
# if not os.path.isdir(LocalData):
#     print("Copying data to "+LocalData)
#     shutil.copytree(os.path.join(root_dir,DataDump[0]),LocalData)
#     print("Done copying data")
# else:
#     print("Found a local copy in "+LocalData)

input = sys.argv[1]

# Load the data from the local directory                                                                                                                                                                   
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# naming plot
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84


# naming sphere container directory
seed = int(root_dir[43:44])
if seed == 1:
    index = 94
elif seed == 2:
    index = 84
sp_container_dir = root_dir[index:]

# make sphere centred on current ss position
width = (w_pccm, 'pc')
r = 2000 # pc
ss_pos, ss_mass, ss_age = ss_properties(ds)
center = ss_pos
tic = time.perf_counter()
sp = ds.sphere(center, 2*r)
disk, L = _make_disk_L(ds, center, 0.5*yt.units.pc, 0.2*yt.units.pc)
toc = time.perf_counter()
print(f"Made the disc in {toc - tic:0.4f} seconds")

# font settings
fontsize = 26
linewidth = 2
plt.rcParams['font.size'] = fontsize
plt.rcParams['font.weight'] = 'light'
rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
rc('text', usetex=True)

plt.rcParams["mathtext.default"] = "regular"
plt.rcParams['lines.linewidth'] = linewidth

# Gas density
if field == "density":
    tic = time.perf_counter()
    field = "number_density"
    p = yt.ProjectionPlot(ds, 'x', ("gas", field), width=width, center=center, data_source=sp,
                          weight_field='density')
    p.set_cmap(field, 'viridis')
    p.set_font_size(fontsize)
    p.set_background_color(("gas", field))
    p.set_axes_unit('pc')

    # annotate
    #p.annotate_scale(corner='lower_right')
    p.annotate_streamlines(("gas", "velocity_y"), ("gas", "velocity_z"), density = 0.7, linewidth=0.6, color='white')
                           #field_color=("gas", "velocity_x")
                           # )
    p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=False)
    p.annotate_text((0.55, 0.94), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    #p.annotate_grids(min_level=18, cmap='Spectral') # not supported in OffAxisProjection
    p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
    p.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
    #p.annotate_cell_edges(line_width=0.00002, alpha=0.7, color='white')
    #p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
    #p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0]/1e3))

    # save
    plot_name = 'density-' + str(root_dir[81:]) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
    p.save('density_plots/' + plot_name)
    print("created density_plots/" + str(plot_name))

# Temperature
elif field == "temperature":
    p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width=width, center=center, data_source=sp,
                           weight_field='density')
    p1.set_cmap('temperature', 'RED TEMPERATURE')
    plot = p1.plots[list(p1.plots)[0]]
    ax = plot.axes
    # nicen up the plot by setting the background color to the minimum of the colorbar
    p1.set_background_color(("gas", "temperature"))
    # hide the axes, while still keeping the background color correct:
    p1.hide_axes(draw_frame=True)
    p1.set_font_size(28)
    p1.set_axes_unit('pc')
    p1.annotate_scale(corner='lower_left')
    plot_name = 'temperature-' + str(root_dir[70:]) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
    p1.save('temperature_plots/' + plot_name)
