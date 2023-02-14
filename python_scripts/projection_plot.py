"""
Plot projection of simulation in density or temperature. Call like:
python projection_plot.py DD0178/DD0178
"""

import yt
import sys
import os
from smartstar_find import ss_properties

# set by user
w_pccm = 8
field = "density"

# set by user
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# make sphere centred on current ss position
width = (w_pccm, 'pccm')
r = 2000 # pc
ss_pos, ss_mass, ss_age = ss_properties(ds, 0.1)
center = ss_pos
sp = ds.sphere(center, 2*r)

fontsize = 18

# Gas density
if field == "density":
    field = "H_nuclei_density"
    p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center, data_source=sp,
                          weight_field='density')
    p.set_cmap(field, 'turbo')
    p.set_font_size(fontsize)
    p.set_background_color(("gas", field))
    p.set_axes_unit('pc')

    # annotate
    p.annotate_scale(corner='lower_left')
    p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
    p.annotate_marker(center, coord_system="data", color="black")  # mark ss position
    p.annotate_text((0.73, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    p.annotate_grids()
    p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
    p.annotate_title("SS Age = {:.2f} kyrs, {} pccm across".format(ss_age[0]/1e3, w_pccm))

    # save
    plot_name = 'density-' + str(root_dir[70:]) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.png'
    p.save('density_plots/' + plot_name)

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
