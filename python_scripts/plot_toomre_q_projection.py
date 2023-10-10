import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.utilities.math_utils import ortho_find
import matplotlib as mpl
import pandas as pd
from matplotlib import rc
from find_disc_attributes import _make_disk_L
from smartstar_find import ss_properties
from plot_radial_profile_from_frb import compute_radial_profile, make_frb, ToomreQ, kappa2D
from matplotlib.colors import LogNorm
from plot_radial_profile_from_frb import extract_dd_segment, extract_simulation_name

fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0167/DD0167"
fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0204/DD0204"
fp = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0250/DD0250"
ds = yt.load(fp)

sim = extract_simulation_name(ds.directory)
dd = extract_dd_segment(ds.directory)


# grab bh particle properties
ss_pos, ss_mass, ss_age = ss_properties(ds)

# constants
G = yt.units.physical_constants.G

# make small disk data container to define L
dx = float(ds.index.get_smallest_dx().in_units('pc'))
disc_r_pc = dx *10*yt.units.pc
disc_h_pc = 0.05*yt.units.pc
print("disk_r_pc = {}, disk_h_pc = {}".format(disc_r_pc, disc_h_pc))
_, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
vecs = ortho_find(L)
north = vecs[1]
dir = vecs[0]

# try density from slice plot x cell height
center = ss_pos 
field = "density"
width_pc = 0.32
pixels = 2048
dx = ds.index.get_smallest_dx().in_units('cm')
p = yt.SlicePlot(ds, dir, ("gas", field), center=center, width=(width_pc, "pc"), north_vector=north)
slc_frb = p.data_source.to_frb((1.0, "pc"), pixels)
slc_dens = slc_frb[("gas", "density")]*dx # surface density
slc_cs = slc_frb[("gas", "sound_speed")].to('cm/s')
slc_kappa = kappa2D(slc_frb)
q = ToomreQ(slc_cs, slc_kappa, G, slc_dens)

# Need to find surface density using FRB and then bin
# npixels = 3040
# center = ss_pos
# frb_height = 0.05*yt.units.pc
# frb_width = 10*yt.units.pc
# cutting = ds.cutting(L, center)
# frb = cutting.to_frb(frb_width, npixels, height=frb_height)
#frb, _ = make_frb(ds, L, ss_pos, width=frb_width, npixels=npixels, height=frb_height)

# Get radius and Toomre Q inputs from frb
# radius = frb['radius'].to('pc')
# surface_density = frb['density'].to('g/cm**3') * frb_height.to('cm') # cm^-2
# cs = frb['sound_speed'].to('cm/s')
# kappa = kappa2D(frb)
# Q = ToomreQ(cs, kappa, G, surface_density)

# Assuming Q is a 2D array with Toomre Q values at each point in the FRB

# Create a figure and axis object
fig, ax = plt.subplots()

# Plot the Toomre Q values
# You can adjust the colormap, color limits, and other parameters as needed
#cax = ax.imshow(Q, cmap='magma', origin="lower", norm=LogNorm())  # adjust extent as necessary
#norm = LogNorm(vmin=Q.min(), vmax=Q.max())
# Convert unyt_quantity to float for the extent
#min_radius = float(-4.98528013 * yt.units.pc)
#max_radius = float(4.98528013 * yt.units.pc)
#ax = ax.imshow(Q, cmap='magma', norm=norm, extent=[np.log10(radius.min()), np.log10(radius.max()), 0, height_value])
#cax = ax.imshow(Q, cmap='magma', norm=norm, extent=[min_radius, max_radius, min_radius, max_radius])
#cax = ax.pcolormesh(radius[0], radius[1], Q, cmap='magma', norm=norm, shading='auto')
cax=ax.imshow(q, cmap='magma', origin="lower", norm=LogNorm())

# Set the scale of the axes to be logarithmic
# ax.set_xscale('log')
# ax.set_yscale('log')

# Optionally, add a colorbar and labels
cbar = fig.colorbar(cax, label='Toomre $Q$', orientation='vertical')  # adjust the label as necessary

# Optionally, set axis limits, #titles, etc.
ax.set_title('Toomre $Q$ at BH Age = {:.2f} Myr'.format(ss_age[0]/1e6))

# set clim
cax.set_clim(1, 1e3)

# Set custom ticks
x_ticks = np.linspace(0, q.shape[1]-1, num=5)  # This creates 6 x ticks evenly spaced across the width of your image
y_ticks = np.linspace(0, q.shape[0]-1, num=5)  # This creates 6 y ticks evenly spaced across the height of your image
ax.set_xticks(x_ticks)
ax.set_yticks(y_ticks)

# Set custom tick labels
x_ticklabels = ['-0.5', '-0.25', '0.', '0.25', '0.5']  # Replace these with your desired x-axis tick labels
y_ticklabels = x_ticklabels
ax.set_xticklabels(x_ticklabels)
ax.set_yticklabels(y_ticklabels)

# Set axis labels
ax.set_xlabel('Radius (pc)')
ax.set_ylabel('Radius (pc)')

# Save the plot
plt.savefig('plots/proj_toomreq_{}_{}.png'.format(sim, dd), dpi=300, bbox_inches='tight')
print("Saved to plots/proj_toomreq_{}_{}.png".format(sim, dd))
