import yt
import ytree
import os
import sys

ds = sys.argv[1]
width = 200
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01"
output = yt.load(os.path.join(root_dir, ds))

# Load merger tree of dataset (up to DD0118 in gas run)                                                                                                                                                  
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
# Load halo                                                                                                                                                                                              
ds_halos = yt.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/halos_36.0.bin')
ds_data = yt.load('/home/sgordon/disk14/pop3/dm-only-L0/DD0036/DD0036')
# Load my_tree and find radius                                                                                                                                                                           
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds.quan(r_halo.d, "pc") # virial radius                                                                                                                                                               

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 5 * virial radius                                                                                               
width = (width, 'pc')
center = 'max'
sp = ds.sphere(center, 5*r)

# Find time since SN                                                                                                                                                                                      
dd_time = ds.current_time.to('Myr').d

# Gas density                                                                                                                                                                                            
p = yt.ProjectionPlot(ds, "x", ("gas", "density"), width=width, center=center, data_source=sp, weight_field='density')
p.set_cmap('density', 'turbo')
p.set_font_size(28)
plot = p.plots[list(p.plots)[0]]
ax = plot.axes

# nicen up the plot by setting the background color to the minimum of the colorbar                                                                                                                      
p.set_background_color(("gas", "density"))
# hide the axes, while still keeping the background color correct:                                                                                                                                       
p.hide_axes(draw_frame=True)
p.set_axes_unit('pc')
p.annotate_scale(corner='lower_left')
plot_name = 'density-' + str(ds) + str(width) + 'pc.pdf'
p.save('plots/' + plot_name)

# Temperature                                                                                                                                                                                            
p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width=width, center=center, data_source=sp,
                       weight_field='density')
p1.set_cmap('temperature', 'RED TEMPERATURE')
plot = p.plots[list(p.plots)[0]]
ax = plot.axes
# nicen up the plot by setting the background color to the minimum of the colorbar                                                                                                                       
p1.set_background_color(("gas", "temperature"))
# hide the axes, while still keeping the background color correct:                                                                                                                                      
p1.hide_axes(draw_frame=True)
p1.set_font_size(28)
p1.set_axes_unit('pc')
p1.annotate_scale(corner='lower_left')
plot_name = 'temperature-' + str(ds) + str(width) + 'pc.pdf'
p1.save('plots/' + plot_name)
