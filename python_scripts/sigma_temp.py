import numpy as np
import yt
from yt.utilities.math_utils import ortho_find
import sys
import yt.utilities
import os
import matplotlib.pyplot as plt

def _height(field, data):
    return np.abs(data["index", "cylindrical_z"])

def ss_properties(ds):
    ad = ds.all_data()
    # find ss properties
    ss_creation = ad['SmartStar', 'creation_time'].to('yr')
    ss_pos = ad['SmartStar', 'particle_position'].to('unitary')[0]
    ss_mass = ad['SmartStar', 'particle_mass'].to('Msun')[0]
    # find ss age
    ss_age = ds.current_time - ss_creation
    return ss_pos, ss_mass, ss_age


def radial_profile(data, center):
    y, x = np.indices((data.shape))
    r = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2)
    r = r.astype(int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile


def radial_profile_app1(data, r):
    mid = data.shape[0]//2
    ids = np.rint((r**2)/r[mid-1,mid]**2).astype(int).ravel()
    print("ids: ", ids)
    count = np.bincount(ids)

    R = data.shape[0]//2 # Radial profile radius
    print("R: ", R)
    R0 = R+1
    dists = np.unique(r[:R0,:R0][np.tril(np.ones((R0,R0),dtype=bool))])
    print("dists: ", dists)

    mean_data = (np.bincount(ids, data.ravel())/count)[count!=0]
    return dists, mean_data


if __name__ == "__main__":
    root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
    input = "DD0148/DD0148"
    ds = yt.load(os.path.join(root_dir, input))
    ds.add_field(("index", "height"), function=_height,
                 sampling_type="local", units="cm")
    pos, mass, age = ss_properties(ds)

    sp = ds.sphere(pos, (1, "pc"))
    L = sp.quantities.angular_momentum_vector()
    L /= np.sqrt((L**2).sum())

    disk = ds.disk(pos, L, (10, "pc"), (0.5, "pc"))
    vecs = ortho_find(L)

    i = 0
    vec = vecs[0]
    north = vecs[0] if i > 0 else vecs[1]
    p = yt.ProjectionPlot(ds, vec, ("gas", "number_density"), weight_field=None,north_vector=north,
                          center=disk.center, width=2*disk.radius, data_source=disk)

    # define sigma_frb and pr
    sigma_frb = p.frb[("gas", "number_density")]
    bds = p.frb.bounds
    shape = p.frb.buff_size
    dx = (bds[1] - bds[0]) / shape[0]
    dy = (bds[3] - bds[2]) / shape[1]
    px, py = np.meshgrid(np.arange((bds[0] + dx / 2), (bds[1] + dx / 2), dx),
                         np.arange((bds[2] + dy / 2), (bds[3] + dy / 2), (dy)))
    pr = ds.arr(np.sqrt(px**2 + py**2), "code_length").to('pc')

    # beckmann's rp method
    surface_density = disk.integrate('density', axis='z')
    dx = 1.229791e-02
    w = 10 * 2 * yt.units.pc
    frb_resolution = int(w / dx)
    surface_density_frb = surface_density.to_frb(width=w, resolution=frb_resolution, center=disk.center)

    xs = np.abs(disk.center[1] - surface_density["px"])
    ys = np.abs(disk.center[2] - surface_density["py"])
    proj_rad = ds.arr(np.sqrt(xs ** 2 + ys ** 2), "code_length").to('pc')

    #rp = radial_profile_app1(sigma_frb, pr)
    #sigma_radii, sigma, frb_radii = radial_profile(surface_density_frb['density'].in_units('amu/cm**2'), bins=list(np.array(pr / dx * frb_resolution / 2)) + [1E10])
    #sigma_radii = sigma_radii * dx / (frb_resolution / 2.0)

    # my method: divide sigma from weighted 1d histogram by r bin counts (works)
    counts_r, r_bin_edges = np.histogram(pr, bins=64)
    sigma, radius = np.histogram(pr, weights=sigma_frb, bins=64)
    sigma = sigma/counts_r

    fig = plt.figure()
    plot_name = "sigma_r.png"
    plt.loglog(radius[:-1], sigma)
    plt.ylabel('$\Sigma \, (cm^{-2})$')
    plt.xlabel('$R \, (pc)$')
    fig.savefig('plots/' + plot_name, dpi=100)