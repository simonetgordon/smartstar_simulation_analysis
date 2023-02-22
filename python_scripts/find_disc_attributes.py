"""
For analysing the nuclear disc feeding a black hole
python find_disc_attributes.py DD0130/DD0130
"""

import yt
import unyt
import sys
import os
from pathlib import Path
import numpy as np
from smartstar_find import ss_properties

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# naming sphere container directory
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84
sp_container_dir = root_dir[index:]
sp_container_path = Path("sphere_containers_{}/{}_sphere.h5".format(sp_container_dir, ds))

# center objects on BH
ss_pos, ss_mass, ss_age = ss_properties(ds)

# make sphere object
width = 0.0001
sp = ds.sphere(ss_pos, width)
# sp.save_as_dataset("sphere_containers_{}/{}_sphere.h5".format(sp_container_dir, ds),
#                    fields=[("gas", "density"), ("SmartStar", "particle_mass"), ("SmartStar", "creation_time"),
#                            ("SmartStar", "particle_position")])
#ad = sp.all_data()

# make disc object
rho_disc = 0 # by density projection plot inspection
print(ss_mass)
disc = ds.cut_region(sp, ["obj['density'].in_units('amu/cm**3') > {0}".format(rho_disc)])
disc_height_sum = disc.sum('dx', axis="x")
print(disc_height_sum)

# Produce a 2D array of uniform pixel sizes of the disc height at the maximum resolution of simulation
sphere_pc = (width * ds.length_unit.in_units("code_length")).in_units("pc")
dx = 7.687095e-04 # pc
print("sphere_pc = ", sphere_pc)
frb_resolution = int(sphere_pc/dx)
print("frb_res = ", frb_resolution)
disc_frb = disc_height_sum.to_frb(width=(2*sphere_pc, 'pc'), resolution=frb_resolution, center=ss_pos)
height_data = disc_frb['dx'].in_units('pc')
print(height_data)


# Used a fixed resolution buffer to grid the height data onto something I could work with. Here “pc” is the total
# size of my sphere in pc and dx is the smallest cell size in the simulation. Sink.pos is the center of my sphere
# frb_resolution=int(pc/dx)
# disc_frb=disc_height_sum.to_frb(width=(2*pc,'pc'),resolution=frb_resolution,center=sink.pos)
# height_data=disc_frb['dz'].in_units('pc')
# If I recall correctly this should give you a 2D array of uniform pixel sizes of the disc height at the maximum
# resolution of your simulation. You could look at it with imshow if you wanted but to make the radial profiles
# I simply binned this data in concentric radial bins centred on the array (making sure I rescaled the x-axis from
# “pixel count” to “pc” using the frb_resolution computed above.