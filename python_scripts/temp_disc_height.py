import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds
from find_disc_attributes import _make_disk_L
from plot_disc_attributes import radial_profile


# data
root_dir =["/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"]
sim = ["1B.RSb01-2"]
dds = ["DD0144/DD0144"]
i = 0
ds = yt.load(os.path.join(root_dir[i], sim[i], dds[i]))
add_fields_ds(ds)
label = str(sim[i]) + "-" + str(float(ds.current_time.to('Myr')))[:5] + "Myr"

# minimum dx
cell_width_pc = [0.012]

# grab bh particle properties
ss_pos, ss_mass, ss_age = ss_properties(ds)

# make disk data container and define L
disc_r_pc = 10*yt.units.pc
disc_h_pc = 1*yt.units.pc
disk, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)

# querying height (abs(cylindrical_z)) on a cut region
n_bins = 64
cuts = [1e6, 4e6, 1e7, 2e7]
hp = plt.figure()
for cut in cuts:
    cr = ds.cut_region(disk, ["obj[('gas', 'number_density')] > {}".format(cut)])
    h_disc = cr[("index", "height")].to('pc')
    r_disc = cr[("index", "radius")].to('pc')
    r_h, h = radial_profile(h_disc, cr, n_bins, cell_width_pc[i])

    hp = plt.loglog(r_h, h, label=r'$> {0:.2E}$'.format(cut))

hp = plt.legend()
hp = plt.title(label + " max disc number density: " + '{0:.2E}'.format(disk[('gas', 'number_density')].d.max()))
hp = plt.ylabel('$H \, (pc)$')
hp = plt.xlabel('$R \, (pc)$')
plot_name = 'disc_height_temp.png'
plt.savefig('plots/' + plot_name, dpi=100)

# Beckmann's method (doesn't work)

# sum dz of cut region (disc height) then make an frb to produce a2D array of uniform pixel sizes of the 
# disc height at the maximum resolution of the simulation.
frb_resolution=int(disk.radius.to('pc')/cell_width_pc[i])
disc_height_sum=cr.sum('dz',axis="z")
disc_frb=disc_height_sum.to_frb(width=2*disk.radius.to('pc'), resolution=frb_resolution, center=disk.center)
height_data=disc_frb['dz'].in_units('pc') # 833x833 ImageArray

# make radius
bds = disc_frb.bounds
shape = disc_frb.buff_size
dx = (bds[1] - bds[0]) / shape[0]
dy = (bds[3] - bds[2]) / shape[1]
px, py = np.meshgrid(np.arange((bds[0] + dx / 2), (bds[1] + dx / 2), dx),
                        np.arange((bds[2] + dy / 2), (bds[3] + dy / 2), dy))
pr = ds.arr(np.sqrt(px ** 2 + py ** 2), "code_length").to('pc')
