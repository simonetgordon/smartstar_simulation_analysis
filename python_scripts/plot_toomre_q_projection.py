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

def toomre_from_sliceplot(ds, disk, center, width_pc, north, dir, npixels=2048):
    """
    Compute Toomre Q from a slice plot of a dataset.
    Surface Density = slice plot density * cell height
    """
    G = yt.units.physical_constants.G
    dx = ds.index.get_smallest_dx().in_units('cm')
    p = yt.SlicePlot(ds, dir, ("gas", "density"), center=center, width=(width_pc, "pc"), data_source=disk)
    slc_frb = p.data_source.to_frb((width_pc, "pc"), npixels)
    slc_dens = slc_frb[("gas", "density")]*slc_frb['index', 'dy'].to('cm') # replaced dx with array of dy
    slc_cs = slc_frb[("gas", "sound_speed")].to('cm/s')
    slc_kappa = kappa2D(slc_frb)
    q = ToomreQ(slc_cs, slc_kappa, G, slc_dens)

    return q

def field_from_sliceplot(field, ds, disk, center, width_pc, north, dir, npixels=2048, radius=False):
    """
    Compute field from a slice plot of a dataset.
    Surface Density = slice plot density * cell height
    """
    p = yt.SlicePlot(ds, dir, ("gas", field), center=disk.center, width=(width_pc, "pc"), data_source=disk)
    slc_frb = p.data_source.to_frb((1.0, "pc"), npixels)
    slc_field = slc_frb[("gas", field)]
    if radius:
        radius = slc_frb['radius'].to('pc')
        return slc_field, radius
    else:
        return slc_field

def toomre_from_frb(ds, center, L, frb_height_pc, frb_width_pc, npixels=2048):
    """
    Compute Toomre Q from a frb of a dataset.
    Surface Density = disk frb densities * disk height from plane of particle.
    """
    frb_height = frb_height_pc
    frb_width = frb_width_pc
    cutting = ds.cutting(L, center)
    frb = cutting.to_frb(frb_width, npixels, height=frb_height)

    #Get radius and Toomre Q inputs from frb
    radius = frb['radius'].to('pc')
    surface_density = frb['density'].to('g/cm**3') * frb_height.to('cm') # cm^-2
    cs = frb['sound_speed'].to('cm/s')
    kappa = kappa2D(frb)
    q = ToomreQ(cs, kappa, G, surface_density)
    return q, radius

if __name__ == "__main__":
    file_paths = [ 
        # "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0167/DD0167",
        # "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0204/DD0204",
        # "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0250/DD0250",
        "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0240/DD0240",
        "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0279/DD0279",
        "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0197/DD0197"
    ]

    for fp in file_paths:
        # Load dataset
        ds = yt.load(fp)

        # Extract relevant data from dataset
        sim = extract_simulation_name(ds.directory)
        dd = extract_dd_segment(ds.directory)
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        
        # Constants and calculations for disk data container
        G = yt.units.physical_constants.G
        dx = float(ds.index.get_smallest_dx().in_units('pc'))
        disc_r_pc = dx * 100 * yt.units.pc
        disc_h_pc = 0.01 * yt.units.pc
        _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
        print(r"Finding L from disc_r = {:.2f} pc and disc_h = {:.2f} pc".format(disc_r_pc, disc_h_pc))
        vecs = ortho_find(L)

        # Set parameters plot parameters
        center = ss_pos 
        field = "density"
        north = vecs[1]
        dir = vecs[0]
        width_pc = 1
        x_ticklabels = ['-0.5', '-0.25', '0.', '0.25', '0.5']
        width_pc = 0.3
        x_ticklabels = ['-0.15', '-0.08', '0.', '0.08', '0.15']

        # Compute Toomre Q from slice plot
        q = toomre_from_sliceplot(ds, center, width_pc, north, dir, npixels=2048)
        #q, _ = toomre_from_frb(ds, center, L, disc_h_pc, disc_r_pc, npixels=2048)

        
        # Create a figure and axis object
        fig, ax = plt.subplots()
        cax = ax.imshow(q, cmap='magma', origin="lower", norm=LogNorm())
        cbar = fig.colorbar(cax, label='Toomre $Q$', orientation='vertical')
        ax.set_title(f'Toomre $Q$ at BH Age = {ss_age[0]/1e6:.2f} Myr')
        cax.set_clim(1, 1e3)
        
        # Set custom ticks and labels
        x_ticks = np.linspace(0, q.shape[1]-1, num=5)
        y_ticks = np.linspace(0, q.shape[0]-1, num=5)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        y_ticklabels = x_ticklabels
        ax.set_xticklabels(x_ticklabels)
        ax.set_yticklabels(y_ticklabels)
        ax.set_xlabel('Radius (pc)')
        ax.set_ylabel('Radius (pc)')
        
        # Save the plot
        figname = r'plots/proj_toomreq_{}_{}_{:.2f}pc_from_frb.png'.format(sim, dd, disc_r_pc)
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        print(f'Saved plot: {figname}')
        plt.close(fig)  # Close the figure to free up memory