import yt
from smartstar_find import ss_pos, ss_mass, ss_age, root_dir, ds, time_array, input

# set by user
w_pccm = 2
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
    p = yt.ProjectionPlot(ds, "x", ("gas", "density"), width=width, center=center, data_source=sp,
                          weight_field='density')
    p.set_cmap('density', 'turbo')
    p.set_font_size(fontsize)
    p.set_background_color(("gas", "density"))
    p.set_axes_unit('pc')

    # annotate
    p.annotate_scale(corner='lower_left')
    p.annotate_marker(center, coord_system="data", color="black")  # mark ss position
    p.annotate_text((0.75, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis", text_args={"color": "black"})
    p.annotate_streamlines(("gas", "velocity_x"), ("gas", "velocity_y"))

    # title
    p.annotate_title("Density plot t = {:.2f}, {} pccm across".format(ds.current_time.to('Myr'), w_pccm))

    # save
    plot_name = 'density-' + str(root_dir[70:]) + '-' + str(input)[10:] + '-' + str(w_pccm) + 'pccm.pdf'
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
    plot_name = 'temperature-' + str(root_dir[70:]) + '-' + str(input[10:]) + '-' + str(w) + 'pccm.pdf'
    p1.save('density_plots/' + plot_name)
