import yt
import ytree
import numpy as np
import matplotlib.pyplot as plt
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters


# Load last star formation snapshot and add pop3/metallicity fields+filters

ds = yt.load("DD0122/DD0122")
add_p2p_fields(ds)
add_p2p_particle_filters(ds)


# Set units of parameters

masses = ds.r['particle_mass'].to('Msun')
positions = ds.r['particle_position'].to('unitary')


# Load merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('../gas+dm-L3/rockstar_halos/out_0.list')


## Make Sphere of Most Massive Halo ##

# Load my_tree
a1 = ytree.load('../gas+dm-L3/tree_810/tree_810.h5')

# Find position and virial radius of the most massive halo
pos_halo = a1[0]["position"].to('unitary')
p = ds.arr(pos_halo.d, "unitary") # position

r_halo = a1[0]["virial_radius"].to('pc')
r = ds.quan(r_halo.d, "pc") # virial radius


# Make sphere containing the most massive halo centred on most dense point, c, with radius = virial radius

v, c = ds.find_max('density') # c is position, v is max density
sp = ds.sphere(c, r)


# Print out Pop III star mass and position and time of formation

star_mass = sp['pop3', 'particle_mass']
star_pos = sp['pop3', 'particle_position']
creation_time = sp['pop3', 'creation_time']

print(star_mass)
print(star_pos)
print(creation_time)
