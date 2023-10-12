import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds
from yt.utilities.math_utils import ortho_find
import matplotlib as mpl
from matplotlib.colors import LogNorm
#from plot_radial_profile_from_frb import extract_simulation_name


# make disc data container
def _make_disk_L(ds, center, width_pc, height_pc):
    width = width_pc*yt.units.pc
    height = height_pc*yt.units.pc
    sp = ds.sphere(center, width)
    L = sp.quantities.angular_momentum_vector()
    L /= np.sqrt((L ** 2).sum()) # normal vector is N = L/|L|
    disk = ds.disk(center, L, width, height)
    return disk, L


if __name__ == "__main__":

    root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/"
    input = sys.argv[1]
    ds = yt.load(os.path.join(root_dir, sys.argv[1]))
    add_fields_ds(ds)
    sim_name = extract_simulation_name(ds.directory)

    # grab bh particle properties
    ss_pos, ss_mass, ss_age = ss_properties(ds)

    # make disk data container and define angular momentum vector L
    disc_r_pc = 0.1
    disc_h_pc = 0.01
    _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
    disc_r_pc = 1.0
    disk = ds.disk(ss_pos, L, disc_r_pc*yt.units.pc, disc_h_pc*yt.units.pc)

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    dir, north, _ = ortho_find(L)

    #north = vecs[0] if i > 0 else vecs[1]
    p = yt.ProjectionPlot(ds, dir, ("gas", "number_density"), weight_field=("gas", "density"), north_vector=north,
                            center=disk.center, width=(0.2, "pc"), data_source=disk)
    p.annotate_marker(ss_pos, coord_system='data', plot_args={'color':'white', 's':100})
    p.annotate_grids(line_color='black', alpha=0.5, min_level=18, max_level=18)
    p.set_axes_unit("pc")
    p.set_cmap(("gas", "number_density"), "viridis")
    p.save("disc_gridlines_{}.pdf".format(sim_name))

    # for field in ("height", "cylindrical_radius"):
    #     p = yt.ProfilePlot(disk, ("index", field), ("gas", "density"), weight_field=("gas", "cell_mass"))
    #     p.set_unit(("index", field), "pc")
    #     p.save()