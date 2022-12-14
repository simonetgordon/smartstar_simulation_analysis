import yt
import ytree
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog
import sys

ds = yt.load(sys.argv[1])
ds_sn = yt.load("../DD0122/DD0122") # DD in which star is formed                                                                                                                                             

# Load merger tree of dataset (up to DD0118 in gas run)                                                                                                                                                  
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list') # after eric moved nov 22 /home/sgordon/disk14/pop3                                                                      
# Load halo                                                                                                                                                                                              
ds_halos = yt.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/halos_36.0.bin')
ds_data = yt.load('/home/sgordon/disk14/pop3/dm-only-L0/DD0036/DD0036')
# Load my_tree and find radius                                                                                                                                                                           
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds.quan(r_halo.d, "pc") # virial radius                                                                                                                                                               

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 5 * virial radius                                                                                               
width = (200, 'pc')
center = 'max' # Star pos at time of formation.                                                                                                                                                           
sp = ds.sphere(center, 5*r)

# Find time since SN                                                                                                                                                                                      
dd_time = ds.current_time.to('Myr').d
sn_time = ds_sn.current_time.to('Myr').d
t_sinceSN = dd_time - sn_time

# Gas density                                                                                                                                                                                            
p = yt.ProjectionPlot(ds, "x", ("gas","density"), width=width, center=center, data_source=sp, weight_field='density')
p.set_cmap('density', 'turbo')
#p.annotate_timestamp(corner="upper_right", redshift=True, draw_inset_box=True)                                                                                                                         
p.set_font_size(28)
#p.set_unit('density', 'Msun')                                                                                                                                                                          
# hide the colorbar:                                                                                                                                                                                    
#p.hide_colorbar()                                                                                                                                                                                       
#p.set_colorbar_label('density')                                                                                                                                                                         
plot = p.plots[list(p.plots)[0]]
ax = plot.axes
# nicen up the plot by setting the background color to the minimum of the colorbar                                                                                                                      
p.set_background_color(("gas", "density"))
# hide the axes, while still keeping the background color correct:                                                                                                                                       
p.hide_axes(draw_frame=True)
p.set_axes_unit('pc')
p.annotate_scale(corner='lower_left')
#p.annotate_text(pos=[-1,9.1], text=f"10.8 Msun Black Hole", coord_system="plot")                                                                                                                       
#p.annotate_text(pos=[-1,8.3], text=f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")                                                                                                        
#p.annotate_text(pos=[-1,7.5], text=f"Smallest dx = Bondi Radius", coord_system="plot")                                                                                                                  
p.save('density_halo_200pc-post-SF--nov-22.png')

# Temperature                                                                                                                                                                                            
p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width=width, center=center, data_source=sp, weight_field='density')
p1.set_cmap('temperature', 'RED TEMPERATURE')
# hide the colorbar:                                                                                                                                                                                     
#p1.hide_colorbar()                                                                                                                                                                                      
plot = p.plots[list(p.plots)[0]]
ax = plot.axes
# nicen up the plot by setting the background color to the minimum of the colorbar                                                                                                                       
p1.set_background_color(("gas", "temperature"))
# hide the axes, while still keeping the background color correct:                                                                                                                                      
p1.hide_axes(draw_frame=True)
p1.set_font_size(28)
p1.set_axes_unit('pc')
p1.annotate_scale(corner='lower_left')
#p1.annotate_text(pos=[-1,9.1], text=f"10.8 Msun Black Hole", coord_system="plot")                                                                                                                       
#p1.annotate_text(pos=[-1,8.3], text=f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#p1.annotate_text(pos=[-1,7.5], text=f"Smallest dx = Bondi Radius", coord_system="plot")                                                                                                                 
p1.save('temp_halo_200pc-post-SF-nov-22.png')
