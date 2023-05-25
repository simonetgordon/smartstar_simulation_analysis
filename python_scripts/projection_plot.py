"""
Plot projection of simulation in density, temperature or dark matter. Call like:
python -i projection_plot.py DD0133/DD0133
"""

import yt
import shutil
import sys
import os
import re # complex str searches
from smartstar_find import ss_properties
from find_disc_attributes import _make_disk_L
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time


def extract_ultimate_directory(filepath):
    # split the filepath into directory and file components
    directory, filename = os.path.split(filepath)
    # split the directory component into its path elements
    path_elements = directory.split(os.path.sep)
    # return the penultimate element, or None if not found
    return path_elements[-1] if len(path_elements) > 1 else None


# set by user
w_pccm = 20
field = "density"

# set by user
root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSm01/"
input = sys.argv[1]

## First check there is a local /data area                                                                                                                                                                 
if os.path.isdir("/data"):
    #sys.exit("Error: no /data")

    ## Second, check we have a directory. If not, then create one.                                                                                                                                             
    UserName=os.getlogin()
    LocalDir=os.path.join("/data",UserName)
    if not os.path.isdir(LocalDir):
        print("Creating Directory "+LocalDir)
        os.mkdir(LocalDir)

    ## Third, check if the data is already there, and if not, copy it over.                                                                                                                                    
    DataDumpFull = sys.argv[1]
    DataDump = DataDumpFull.split('/')
    LocalData = os.path.join(LocalDir,DataDump[0])
    if not os.path.isdir(LocalData):
        print("Copying data to "+LocalData)
        shutil.copytree(os.path.join(root_dir,DataDump[0]),LocalData)
        print("Done copying data")
    else:
        print("Found a local copy in "+LocalData)
    root_dir = LocalDir

# Load the data from the local directory                                                                                                                                                                   
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# naming plot
# Extract the string after "seed" and before the next "-" character
seed = re.search(r"seed(\d+)-", root_dir)

# Extract sim name as the penultimate dir in the filepath
sim = extract_ultimate_directory(root_dir)

print(seed)

# make sphere centred on current ss position
width = (w_pccm, 'pc')
r = 2000 # pc
ss_pos, ss_mass, ss_age = ss_properties(ds)
center = ss_pos
sp = ds.sphere(center, 2*r)
disk, L = _make_disk_L(ds, center, 0.5*yt.units.pc, 0.2*yt.units.pc)

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
    field = "number_density"
    p = yt.ProjectionPlot(ds, 'x', ("gas", field), width=width, center=center, data_source=sp,
                            weight_field='density')
    p.set_cmap(field, 'viridis')
    p.set_font_size(fontsize)
    p.set_background_color(("gas", field))
    p.set_axes_unit('pc')

    # annotate
    #p.annotate_scale(corner='lower_right')
    #p.annotate_streamlines(("gas", "velocity_y"), ("gas", "velocity_z"), density = 0.7, linewidth=0.6, color='yellow')
                            #field_color=("gas", "velocity_x")
                            # )
    p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=False)
    p.annotate_text((0.55, 0.94), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    p.annotate_grids(min_level=18, cmap='Spectral') # not supported in OffAxisProjection
    p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
    p.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
    #p.annotate_cell_edges(line_width=0.00002, alpha=0.7, color='white')
    #p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
    p.annotate_title("High Resolution Region")

    # save
    plot_name = 'density-' + str(sim) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
    p.save('density_plots/' + plot_name)
    print("created density_plots/" + str(plot_name))
    field = "dm"
    plt.show()

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

#elif field == "dm":
plt.show()
plt.figure()

#p1 = yt.ProjectionPlot(ds, 'x', "all_cic", width=width, center=center, data_source=sp)
p2 = yt.ParticlePlot(ds, ("all", "particle_position_y"), ("all", "particle_position_z"), ("all", "particle_mass"), 
                     width=width, center=center, data_source=sp)
#p1.set_cmap("all_cic", 'viridis')
p2.set_font_size(fontsize)
p2.set_axes_unit('pc')
p2.set_unit(("all", "particle_mass"), "Msun")

p2.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
p2.annotate_text((0.55, 0.94), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                text_args={"color": "black"})

p2.annotate_marker(center, coord_system="data", color="black")  # mark bh position
p1.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
p2.annotate_title("DM in High Resolution Region")

# save
plot_name = 'dm-particles' + str(sim) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
p2.save('density_plots/' + plot_name)
print("created density_plots/" + str(plot_name))
