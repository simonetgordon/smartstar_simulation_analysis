import ytree
import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
from yt.units import kpc
from yt.units import pc

ds_10_1 = yt.load("DD0133/DD0133") # t = 134.3, 5 Myr forward from dd nearest supernova explosion
ds_5_1 = yt.load("../star_raysPerCell_5.1/DD0142/DD0142") # 5 Myr post supernova

all_data = [ds_10_1, ds_5_1]

## Make Sphere of Most Massive Halo ##

# Load merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('../gas+dm-L3/rockstar_halos/out_0.list')

# Load my_tree and find radius
a1 = ytree.load('../gas+dm-L3/tree_810/tree_810.h5')

r_halo = a1[0]["virial_radius"].to('pc')
r = ds_10_1.quan(r_halo.d, "pc") # virial radius

# Make initial sphere centred on the star at the time of formation (DD0122) with radius = 5 * virial radius
star_pos0 = [0.49053251, 0.4946698, 0.50963437]


for ds in all_data:

    # Add p2p fields and filters
    add_p2p_fields(ds)
    add_p2p_particle_filters(ds)

    sp_a = ds.sphere(star_pos0, 5*r)

    star_pos = sp_a['pop3', 'particle_position']

    # sphere around current star position
    sp_1 = ds.sphere(star_pos, 3*r)

    p = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width= 1*kpc ,center=star_pos,  weight_field='density')
    p.set_cmap('metallicity3', 'kamae')
                                                                                                                                                                    
    p.set_axes_unit('kpccm')
    p.annotate_timestamp(corner='lower_right')
    p.annotate_scale(corner='lower_left')
    p.save()
    
    p.zoom(2)
    p.save('zoom')
  

    
