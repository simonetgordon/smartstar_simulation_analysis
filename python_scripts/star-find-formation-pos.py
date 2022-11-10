import yt
import ytree
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
import os


# Load last star formation snapshot and add pop3/metallicity fields+filters
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/BF-0.0078125"
ds = yt.load(os.path.join(root_dir, "DD0129/DD0129"))
add_p2p_fields(ds)
#add_p2p_particle_filters(ds)

# Set units of parameters
masses = ds.r['particle_mass'].to('Msun')
positions = ds.r['particle_position'].to('unitary')

# Load seed1 merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
r_halo = a1[0]["virial_radius"].to('pc')
r = ds.quan(r_halo.d, "pc") # virial radius


"""Make Sphere of Most Massive Halo"""

# Find position and virial radius of the most massive halo
pos_halo = a1[0]["position"].to('unitary')
p = ds.arr(pos_halo.d, "unitary") # position

r_halo = a1[0]["virial_radius"].to('pc')
r = ds.quan(r_halo.d, "pc") # virial radius

# Make sphere containing the most massive halo centred on most dense point, c, with radius = virial radius
v, c = ds.find_max('density') # c is position, v is max density
sp = ds.sphere(c, r)

# Print out Pop III star mass and position and time of formation
star_mass = sp['SmartStar', 'particle_mass']
star_pos = sp['SmartStar', 'particle_position']
creation_time = sp['SmartStar', 'creation_time']

print("==========================================")
print("BH Mass: ", star_mass)
print("position: ", star_pos)
print("creation time: ", creation_time)
print("==========================================")
