import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from helper_functions import ss_properties
import sys
from helper_functions import _make_disc_L
from yt.utilities.math_utils import ortho_find
from helper_functions import tidy_data_labels
from yt.units import pc
from matplotlib import rc
import re
import matplotlib.colors as colors
from helper_functions import ToomreQ, extract_dd_segment, extract_simulation_name, kappa2D, kappa2D_vmag, compute_radial_profile


def setup_plot_env(fontsize, linewidth):
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth


if __name__ == "__main__":

    # Define the orientation direction and load the dataset
    fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/DD0440/DD0440"
    #fp = "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/DD0233/DD0233"
    #fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0204/DD0204"
    fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0375/DD0375"
    ds = yt.load(fp)
    dir == "Lface"
    sim_name = extract_simulation_name(ds.directory)
    dd = extract_dd_segment(ds.directory)
    # Constants
    G = yt.units.physical_constants.G # gravitational constant in cgs units    

    # Find the properties of the BH
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    center = "max" if ss_pos is None else ss_pos

    # Define the rxadius and height of the disk for calculating L
    disc_r_pc = 0.2
    disc_h_pc = 0.05

    # Create a disk data container
    disk, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)

    # Cut the disk data container
    frb_width = 20*yt.units.pc
    frb_height = 0.05*yt.units.pc
    npixels = 3040
    frb, _ = make_frb(ds, L, center, width=frb_width, npixels=npixels, height=frb_height)
    dx = frb['dx'].to('pc').min()

    # Get x and y data
    radius = frb['radius'].to('pc')
    density = frb['density'].to('g/cm**3')
    #surface_density = number_density * frb['dx'].to('cm')**2
    surface_density = density * frb_height.to('cm') # cm^-2
    cs = frb['sound_speed'].to('cm/s')
    kappa = kappa2D(frb)
    kappa_vmag = kappa2D_vmag(frb)
    Q = ToomreQ(cs, kappa, G, surface_density)
    Q_vmag = ToomreQ(cs, kappa_vmag, G, surface_density)

    # Plot 2x2 surface density + toomre Q 
    radii_s, sigma = compute_radial_profile(radius, surface_density)
    radii_q, q = compute_radial_profile(radius, Q)
    radii_q2, q2 = compute_radial_profile(radius, Q_vmag)

    setup_plot_env(fontsize=12, linewidth=3)
    plt.figure(figsize=(8, 4))
    plt.subplot(1, 2, 1)
    plt.loglog(radii_s, sigma, linestyle='-', color='turquoise')
    #plt.xlim(0, 1)
    plt.xlabel('Radius (pc)')
    plt.ylabel(r'Surface Density [$\rm g \, cm^{-2}$]')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.axvline(x=dx*4, color='k', linestyle='dotted', linewidth=2, label=r'4 dx')
    # Adding text to the top right corner
    plt.text(0.50, 0.82, 'dx = {:.0e}\n BH age = {:.2f} Myr\n BH mass = {:.0f}'.format(dx, ss_age[0]/1e6, ss_mass), 
                horizontalalignment='left', verticalalignment='bottom', 
                transform=plt.gca().transAxes,  # to make sure the text is positioned relative to the axis
                color='darkblue', bbox=dict(facecolor='white', edgecolor='white', boxstyle='square,pad=0.3', alpha=0.5), 
                fontsize=11)
    plt.title('Surface Density')

    # Log-log plot
    plt.subplot(1, 2, 2)
    plt.loglog(radii_q, q, linestyle='-', color='purple', label=r'$v_{\rm cyl}$')
    plt.loglog(radii_q2, q2, linestyle='-', color='green', label=r'$v_{\rm mag}$')
    plt.xlabel('Radius (pc)')
    plt.ylabel(r'Toomre Q')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.axvline(x=dx*4, color='k', linestyle='dotted', linewidth=2, label=r'4 dx')
    plt.axhline(y=1, color='k', linestyle='dashed', linewidth=1)
    plt.legend(loc='lower right', frameon=True)
    plt.title(r'$Q > 1$ is stable')

    # Save figure
    plt.tight_layout()
    plt.savefig("plots_frb/surface_density_toomreq_vs_radius_{}_{}.png".format(sim_name, dd))
    print("Saved to plots_frb/surface_density_toomreq_vs_radius_{}_{}.png".format(sim_name, dd))
