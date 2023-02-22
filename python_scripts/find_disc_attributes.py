import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties

# set by user
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# naming sphere container directory
seed = int(root_dir[43:44])
if seed == 1:
    index = 82
elif seed == 2:
    index = 84
sp_container_dir = root_dir[index:]

# center sphere on BH
ss_pos, ss_mass, ss_age = ss_properties(ds, 0.1, sp_container_dir)
# make sphere
width = (0.2, 'unitary')
width = (width, 'unitary')
sp = ds.sphere(ss_pos, width).save_as_dataset("sphere_containers_{}/{}_sphere.h5".format(ds, sp_container_dir),
                   fields=[("gas", "density"), ("SmartStar", "particle_mass"), ("SmartStar", "creation_time"),
                           ("SmartStar", "particle_position")])
sp = yt.load(sp)
ad = sp.all_data()