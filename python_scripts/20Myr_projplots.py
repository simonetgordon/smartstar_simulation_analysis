import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters

ds_10_3 = yt.load("DD0148/DD0148") # 20 Myr after star formation

ds_7_3 = yt.load("../star_raysPerCell_7/DD0172/DD0172") # 20 Myr
ds_5_3 = yt.load("../star_raysPerCell_5.1/DD0172/DD0172")

DD = [ds_10_3, ds_7_3, ds_5_3]

for ds in DD:

    # Add p2p fields and filters
    add_p2p_fields(ds)
    add_p2p_particle_filters(ds)

## Make Sphere of Most Massive Halo ##

# Load merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('../gas+dm-L3/rockstar_halos/out_0.list')

# Load my_tree and find radius
a1 = ytree.load('../gas+dm-L3/tree_810/tree_810.h5')

r_halo = a1[0]["virial_radius"].to('pc')
r = ds_10_3.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 5 * virial radius

star_pos0 = [0.49053251, 0.4946698, 0.50963437]

sp10_c = ds_10_3.sphere(star_pos0, 5*r)
sp7_c = ds_7_3.sphere(star_pos0, 5*r)
sp5_c = ds_5_3.sphere(star_pos0, 5*r)

# Get Pop III star position and make sphere centred on current star

star_pos10_3 = sp10_c['pop3', 'particle_position']
star_pos7_3 = sp7_c['pop3', 'particle_position']
star_pos5_3 = sp5_c['pop3', 'particle_position']

# sphere around current star position
sp10_3 = ds_10_3.sphere(star_pos10_3, 3*r)
sp7_3 = ds_7_3.sphere(star_pos7_3, 3*r)
sp5_3 = ds_5_3.sphere(star_pos5_3, 3*r)

# data_source = sp10_3, (removed)
ds_10_3.all_data()

p = yt.ProjectionPlot(ds_10_3, "x", ("gas", "temperature"), width=(200, "pc"), center= [0.48897534, 0.49257759, 0.50813257], weight_field='density')
p.set_cmap('temperature', 'RED TEMPERATURE')
p.save('20Myr_temp_rays10.png')

p1 = yt.ProjectionPlot(ds_7_3, "x", ("gas", "temperature"), width=(200, "pc"), center= [0.48894602, 0.49261345, 0.50813122], weight_field='density')
p1.set_cmap('temperature', 'RED TEMPERATURE')
p1.save('20Myr_temp_rays7.png')

p2 = yt.ProjectionPlot(ds_5_3, "x", ("gas", "temperature"), width=(200, "pc"), center = [0.48894891, 0.49261072, 0.50812105], weight_field='density')
p2.set_cmap('temperature', 'RED TEMPERATURE')
p2.save('20Myr_temp_rays5.png')
