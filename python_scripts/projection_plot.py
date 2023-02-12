"""
Find and print SmartStar properties. Call like:
python density_plot_zoomed_new.py
"""

import yt
from smartstar_find import ss_pos, ss_mass, ss_age, root_dir, ds, time_array, input

# set by user
w_pccm = 8
field = "density"

# make sphere centred on current ss position
width = (w_pccm, 'pccm')
r = 2000 # pc
center = ss_pos
sp = ds.sphere(center, 2*r)

# format plot
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
    p.annotate_timestamp(corner='lower_right')
    p.annotate_marker(center, coord_system="data", color="black")  # mark ss position
    p.annotate_text((0.73, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"})
    p.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))
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
