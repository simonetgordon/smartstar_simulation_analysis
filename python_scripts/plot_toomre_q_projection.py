import yt
import numpy as np
import matplotlib.pyplot as plt
from yt.utilities.math_utils import ortho_find
import matplotlib as mpl
import pandas as pd
from matplotlib import rc
from helper_functions import ss_properties
from plot_radial_profile_from_frb import compute_radial_profile, make_frb, ToomreQ, kappa2D
from matplotlib.colors import LogNorm
from plot_radial_profile_from_frb import extract_dd_segment, extract_simulation_name
from helper_functions import _make_disk_L,field_from_sliceplot,toomre_from_sliceplot


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