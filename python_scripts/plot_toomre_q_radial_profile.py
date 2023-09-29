import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from smartstar_find import ss_properties
import sys
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find
from plot_multi_projections import tidy_data_labels

# Define the orientation direction and load the dataset
ds = yt.load("/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/DD0440/DD0440")
dir == "Lface"

# Find the properties of the BH
ss_pos, ss_mass, ss_age, ss_vel = ss_properties(ds, velocity=True)
center = "max" if ss_pos is None else ss_pos

# Define the radius and height of the disk
disc_r_pc = (0.2, 'pc')
disc_h_pc = (0.05, 'pc')

# Create a disk data container
disk, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)

# Set north vector
vecs = ortho_find(L)
direction = vecs[0] if dir == "Lface" else vecs[1] if dir == "Ledge" else vecs[2]
north_vector = vecs[1] if dir == "Lface" else vecs[2] if dir == "Ledge" else vecs[0]

## TESTING ##
# Project the density along the z-axis (or the axis perpendicular to the disk)
proj = ds.proj('density', 'z', data_source=disk, weight_field='density')

# You can convert this projection to a fixed resolution buffer for plotting or analysis
frb = proj.to_frb((0.2, 'pc'), 1024)

# The 'frb' now contains the surface density (integrated density along the line of sight)
# You can access it as a 2D NumPy array:
surface_density = frb['density']