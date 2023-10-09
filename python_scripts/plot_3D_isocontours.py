from yt.units import pc
import yt
import numpy as np
import matplotlib.cm as cm
from smartstar_find import ss_properties
import sys
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find  
from plot_toomre_q_radial_profile import extract_simulation_name, extract_dd_segment

"""
make sketchfab 3D rendered plot of the disc 
"""

# Define the orientation direction and load the dataset
fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/DD0440/DD0440"
fp = "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/DD0233/DD0233"
fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0182/DD0182"
ds = yt.load(fp)
sim_name = extract_simulation_name(ds.directory)
dd = extract_dd_segment(ds.directory)

# Define the orientation direction and load the dataset
ss_pos, ss_mass, ss_age = ss_properties(ds)
_, L = _make_disk_L(ds, ss_pos, width_pc=0.1, height_pc=0.01)
disk = ds.disk(ss_pos, L, (1, "pc"), (0.5, "pc"))
rho = 9e5
surface = ds.surface(disk, ("gas", "number_density"), rho)
map = "temperature"
temp_max = surface[("gas", "temperature")].max()
temp_min = surface[("gas", "temperature")].min()
bounds = [[disk.center[i] - 0.1 * pc, disk.center[i] + 0.1 * pc] for i in range(3)]
sketchfab_api_key="11871f901e104be38dcf06b43f5b06ab"
upload_id = surface.export_sketchfab(
    title="{}_{}_BHAge={:.2f}Myr".format(sim_name, dd, ss_age[0]/1e6),
    description="Extraction of density surface of disc at {:.2e} cm^-3, with a {} map of max value={:.2f} and min value={:.2f}.".format(rho, map, temp_max, temp_min),
    api_key=sketchfab_api_key,
    color_field=('gas', map),
    color_map="hot",
    color_log=True,
    bounds=bounds,
    # color_field_max=3000, 
    # color_field_min=30
)