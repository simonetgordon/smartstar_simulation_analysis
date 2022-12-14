import sys
import yt
import ytree
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
from yt.extensions.astro_analysis.halo_analysis import HaloCatalog

ds = yt.load(sys.argv[1])
ds_sn = yt.load("../DD0126/DD0126") # DD in which BH is formed

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
width = (500, 'pc')
center = 'max' # Star pos at time of formation.                                                                                                                                                 
sp = ds.sphere(center, 5*r)

# Find time since SN                                                                                                                                                                              
dd_time = ds.current_time.to('Myr').d
sn_time = ds_sn.current_time.to('Myr').d
t_sinceSN = dd_time - sn_time


# Add Pop III metallicity fields                                                                                                                                                                  
add_p2p_fields(ds)
add_p2p_particle_filters(ds)

p2 = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width=width,center=center, data_source=sp, weight_field='density')
p2.set_cmap('metallicity3', 'kamae')
p2.hide_axes(draw_frame=True)
p2.set_font_size(28)
#p.set_unit('particle_mass', 'Msun')                                                                                                                                                      \       
#p2.annotate_text(pos=[-10,220], text=f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
p2.set_axes_unit('pc')
#p.annotate_timestamp(corner='lower_right')                                                                                                                                                       
p2.annotate_scale(corner='lower_left')
p2.save("metal_halo_500pc-post-SN.png")


# Temperature                                                                                                                                                                                      
p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width=width, center=center, data_source=sp, weight_field='density')
p1.set_cmap('temperature', 'RED TEMPERATURE')                                                                                                                                                             
# nicen up the plot by setting the background color to the minimum of the colorbar                                                                                                                
p1.set_background_color(("gas", "temperature"))
# hide the axes, while still keeping the background color correct:                                                                                                                                
p1.hide_axes(draw_frame=True)
p1.set_font_size(28)
p1.set_axes_unit('pc')
p1.annotate_scale(corner='lower_left')
#p1.annotate_text(pos=[-10,220], text=f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
p1.save('temp_halo_500pc-post-SN.png')
