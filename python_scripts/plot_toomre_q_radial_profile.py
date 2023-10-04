import yt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from smartstar_find import ss_properties
import sys
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find
from plot_multi_projections import tidy_data_labels
from yt.units import pc
from plot_multipanel_velocity import extract_simulation_name
import re

def extract_dd_segment(file_path: str) -> str:
    """
    Extracts the 'DDxxxx' segment from a given file path.

    Parameters:
    file_path (str): The file path from which to extract the 'DDxxxx' segment.

    Returns:
    str: The 'DDxxxx' segment if it exists, otherwise an empty string.
    """
    # Define a regular expression pattern to find 'DDxxxx' where xxxx are numbers
    pattern = re.compile(r'DD[0-9]{4}')
    
    # Search for the pattern in the file path
    match = pattern.search(file_path)
    
    # If a match is found, return it; otherwise return an empty string
    if match:
        return match.group()
    else:
        return ""
    
def bin_data(radial_bins, indices, data):
    """
    Bin the data by radius.
    """
    bins = radial_bins

    # Assume data is a 2D array
    data_flat = data.flatten()

    # Initialize an array to hold the mean data in each bin
    mean_data = np.zeros_like(bins)

    # Compute mean surface density in each bin
    for i in range(1, len(bins)):
        mask = indices == i
        if np.sum(mask) > 0:  # Avoid division by zero
            mean_data[i] = np.mean(data_flat[mask])

    return mean_data

def kappa2D(frb, G):
    """
    Calculate the epicyclic frequency kappa = v/r
    """
    kappa = frb['velocity_magnitude'] / frb['radius'] # try cylindrical theta
    return kappa

def ToomreQ(cs, kappa, G, surface_density):
    """
    Calculate the Toomre Q parameter for linear stability
    """
    Q = cs * kappa / (np.pi * G * surface_density*0.6*1.67e-24)
    return Q

def radial_bins(ds, width, npixels, num_bins=50):
    """
    Create radial bins for the radial profile.
    """
    # Define the physical size of the image
    min_cell_width = ds.index.get_smallest_dx().to('pc')
    max_radius_pc = float(width.d)             # maximum radius in parsecs
    min_radius_pc = float(min_cell_width.d)    # minimum radius in parsecs

    # Create a 2D coordinate grid
    y, x = np.ogrid[-npixels/2:npixels/2, -npixels/2:npixels/2]

    # Normalize the grid to have values ranging from -1 to 1
    x = x / (npixels/2)
    y = y / (npixels/2)

    # Calculate the radius at each point in the normalized grid
    radius_normalized = np.sqrt(x**2 + y**2)

    # Rescale the normalized radius to represent physical distances in parsecs
    radius_pc = radius_normalized * (max_radius_pc - min_radius_pc) + min_radius_pc

    # Now, radius_pc is a 2D array with the radius (in parsecs) at each pixel.
    # Flatten the 2D arrays to 1D arrays
    radius_flat = radius_pc.flatten()

    # Define the radial bins
    bins = np.logspace(np.log10(min_radius_pc), np.log10(max_radius_pc), num_bins)

    # Digitize the data: assign each radius to a bin
    indices = np.digitize(radius_flat, bins)

    return bins, indices


def make_frb(ds, L, center, width= 10*yt.units.pc, npixels=1024, height=0.05*yt.units.pc):
    """
    Make a fixed resolution buffer (frb) of the disk from its center and 
    angular momentum vector L.
    """
    cutting = ds.cutting(L, center)
    frb = cutting.to_frb(width, npixels, height=height)
    return frb, height


def Toomre_Q_profile(ds, ss_pos, L, frb_width=10*yt.units.pc, frb_height=0.05*yt.units.pc, npixels=1024, G=6.6743e-8):
    # Find the properties of the BH
    center = "max" if ss_pos is None else ss_pos

    # Cut the disk data container
    frb, _ = make_frb(ds, L, center, width=frb_width, npixels=1024, height=frb_height)

    # Get surface number density map
    number_density = np.array(frb['number_density'])
    #surface_density = number_density * frb['dx'].to('cm')**2
    surface_density = number_density * frb_height.to('cm') # cm^-2

    # Get sound speed map
    cs = np.array(frb['sound_speed'])

    ### Make radial profile of the surface density ###
    rad_bins, indices = radial_bins(ds, frb_width, npixels, num_bins=64)

    # Get kappa
    kappa = kappa2D(frb, G)

    # weight the data by the density?
    # weighted_hist, bin_edges = np.histogram(data, bins=bins, weights=data*density, density=False)
    # sum_weights, _ = np.histogram(data, bins=bins, weights=density, density=False)
    # weighted_mean = weighted_hist / sum_weights

    # Compute the mean surface density and sound speed in each bin
    surface_density = bin_data(rad_bins, indices, data=surface_density)
    cs = bin_data(rad_bins, indices, data=cs)
    kappa = bin_data(rad_bins, indices, data=kappa)

    # Compute Toomre Q parameter
    toomre_q = ToomreQ(cs, kappa, G, surface_density)
    toomre_q = bin_data(rad_bins, indices, data=toomre_q)

    return rad_bins, toomre_q

if __name__ == "__main__":

    # Define the orientation direction and load the dataset
    fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/DD0440/DD0440"
    fp = "/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/DD0373/DD0373"
    fp = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0204/DD0204"
    ds = yt.load(fp)
    dir == "Lface"
    sim_name = extract_simulation_name(ds.directory)
    dd = extract_dd_segment(ds.directory)
    # Constants
    G = 6.6743e-8  # gravitational constant in cgs units    

    # Find the properties of the BH
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    center = "max" if ss_pos is None else ss_pos

    # Define the rxadius and height of the disk
    disc_r_pc = 0.2
    disc_h_pc = 0.05

    # Create a disk data container
    disk, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)

    # Cut the disk data container
    frb_width = 10*yt.units.pc
    frb_height = 0.05*yt.units.pc
    npixels = 2048
    frb, _ = make_frb(ds, L, center, width=frb_width, npixels=1024, height=frb_height)

    # Get surface number density map
    number_density = np.array(frb['number_density'])
    #surface_density = number_density * frb['dx'].to('cm')**2
    surface_density = number_density * frb_height.to('cm') # cm^-2

    # Get sound speed map
    cs = np.array(frb['sound_speed'])

    ### Make radial profile of the surface density ###
    num_bins = 128
    radial_bins, indices = radial_bins(ds, frb_width, npixels, num_bins=num_bins)

    # Get kappa
    kappa = kappa2D(frb, G)

    # Compute the mean surface density and sound speed in each bin
    surface_density = bin_data(radial_bins, indices, data=surface_density)
    cs = bin_data(radial_bins, indices, data=cs)
    kappa = bin_data(radial_bins, indices, data=kappa)

    # Compute Toomre Q parameter
    tq = ToomreQ(cs, kappa, G, surface_density)

    # Plotting the radial profile
    plt.figure(figsize=(6, 4.5))
    plt.loglog(radial_bins, tq, marker='o', linestyle='-', color='purple')
    plt.xlabel('Radius (pc)')
    #plt.ylabel(r'Mean Surface Density $cm^{-2}$')
    plt.ylabel(r'Toomre Q')
    plt.title('Toomre Q Stability Criterion: Q > 1 is stable to linear perturbations')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.show()
    plt.savefig("radial_profile_toomreQ_{}_{}.png".format(sim_name, num_bins))


    # make sketchfab 3D rendered plot of the disc 
    # dd = disk
    # surface = ds.surface(disk, ("gas", "number_density"), 1e6 )
    # surface = ds.surface(disk, ("gas", "number_density"), 1e7 )
    # bounds = [[dd.center[i] - 1 * pc, dd.center[i] + 1 * pc] for i in range(3)]
    # sketchfab_api_key="11871f901e104be38dcf06b43f5b06ab"
    # upload_id = surface.export_sketchfab(
    #     title="1S.m04-no-SN_DD0440",
    #     description="Extraction of dendity surface of disc at 1e7 cm^-3,cylindrical velocity map.",
    #     color_field=('gas', 'velocity_cylindrical_radius'),
    #     color_map="hot",
    #     color_log=True,
    #     bounds=bounds,
    #     api_key=sketchfab_api_key,
    # )