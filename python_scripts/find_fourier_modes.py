import numpy as np
import yt
from scipy.interpolate import griddata
from plot_toomre_q_projection import field_from_sliceplot
from helper_functions import ss_properties
from helper_functions import _make_disk_L
from helper_functions import tidy_data_labels
from yt.utilities.math_utils import ortho_find
from plot_radial_profile_from_frb import extract_dd_segment, extract_simulation_name


def find_bar_radius(phi_2_values, radii, var_deg=5):
    """
    Find the radius at which the standard deviation std on the phase angle phi_2_values exceeds var_deg.
    phi_2_values is the phase angle for each radius.
    radii is the array of radii values.
    var_deg is the variability threshold in degrees.
    Output: bar_radius, i (bar radius and index of bar radius in radii array)
    """
    # find the incremental standard deviation in phi_2_values
    std = np.zeros(len(phi_2_values))
    # once the standard deviation exceeds var_deg, we have found the bar radius
    for i in range(len(phi_2_values)):
        std[i] = np.std(phi_2_values[:i])
        if std[i] > var_deg:
            bar_radius = radii[i]
            break
    return bar_radius, i


def calculate_theta_array(dimensions):
    # Create an empty array of the same dimensions
    theta_array = np.zeros(dimensions)

    # Calculate the center of the array
    center_x, center_y = dimensions[0] // 2, dimensions[1] // 2

    # Iterate over each pixel
    for x in range(dimensions[0]):
        for y in range(dimensions[1]):
            # Calculate the relative positions
            rel_x, rel_y = x - center_x, y - center_y

            # Calculate the angle and adjust to ensure 0 is 'north'
            theta = np.arctan2(rel_y, rel_x) + np.pi / 2

            # Normalize the angle to be between 0 and 2pi
            theta = theta % (2 * np.pi)

            # Assign the calculated theta to the array
            theta_array[x, y] = theta

    return theta_array


def get_phi_values(surface_density_polar):
    """
    Calculate azimuthal angles for a given polar surface density array.
    
    Parameters:
    - surface_density_polar: 2D array with rows representing unique radii and 
      columns representing data points equally spaced in azimuthal angle.
    
    Returns:
    - phi_values: 1D array of azimuthal angles.
    """
    num_phi = surface_density_polar.shape[1]  # Number of azimuthal data points
    phi_values = np.linspace(0, 2*np.pi, num_phi, endpoint=False)  # exclude 2*pi to avoid redundancy
    
    return phi_values


def cartesian_to_polar_grid(surface_density, radius):
    """
    Convert a 2D surface density map from Cartesian (x,y) to a polar grid (r,phi) where each 
    row corresponds to a constant radius and each column corresponds to a constant azimuthal angle.
    
    Parameters:
    - surface_density: 2D array, representing the surface density at each radius.
    - radius: 2D array, representing the radius at each (x, y) location.

    Returns:
    - surface_density_polar: 2D array, representing the surface density in the polar grid.
    - r_values: Radii values.
    - phi_values: Azimuthal angles.
    """
    
    # Calculate r and phi values for each (x, y) point
    Phi = get_phi_values(surface_density)
    R = radius
    
    # Define the new grid in polar coordinates
    num_r = surface_density.shape[0]
    num_phi = surface_density.shape[1]
    r_values = np.linspace(0, R.max(), num_r)
    phi_values = np.linspace(0, 2*np.pi, num_phi)
    R_polar, Phi_polar = np.meshgrid(r_values, phi_values, indexing="ij")
    
    # Convert data points into 1D array
    points = np.column_stack((R.ravel(), Phi.ravel()))
    values = surface_density.ravel()
    
    # Interpolate onto the new grid
    surface_density_polar = griddata(points, values, (R_polar, Phi_polar), method='cubic')
    
    return surface_density_polar, r_values, phi_values


def compute_m2_mode(surface_density_polar, phi):
    """
    Calculate the m=2 Fourier mode for the given polar surface density array.
    
    Parameters:
    - surface_density_polar: 2D array, where rows correspond to unique radii and columns
      correspond to data points equally spaced in azimuthal angle.
    - phi: 1D array of azimuthal angles.
    
    Returns:
    - A2: 1D array of m=2 Fourier mode amplitudes for each radius.
    """
    
    # Calculate the complex Fourier coefficient
    integrand = surface_density_polar * np.exp(-2j * phi)
    
    # Integrate over azimuthal angle for each radius using the trapezoidal rule
    A2 = np.trapz(integrand, phi, axis=1)
    
    return np.abs(A2)


def compute_bar_strength(densities, theta, dA):
    """
    Compute the bar strength and phase angle variability given the densities 
    of cells and their azimuthal coordinates.
    
    Parameters:
    - densities: 2D array of densities for each cell in the image.
    - theta: 1D array of azimuthal coordinates for each column of the 2D density array (in radians).
    - dA: Area of each cell.
    
    Returns:
    - bar_strength: Computed bar strength.
    - phi_2: Phase angle for each point.
    - phi_2_var: Variability measure of phase angle.
    """
    
    # Compute the mass-equivalent for each cell
    mass_equivalent = densities * dA
    
    # Expand theta into a 2D array matching the shape of densities for multiplication
    theta_2d = np.tile(theta, (densities.shape[0], 1))
    
    # Compute a_2 and b_2 coefficients for m=2 mode
    a_2 = np.sum(mass_equivalent * np.cos(2 * theta_2d))
    b_2 = np.sum(mass_equivalent * np.sin(2 * theta_2d))
    
    # Compute A_0, which is just the sum of all the mass-equivalents
    A_0 = np.sum(mass_equivalent)
    
    # Compute A_2 using the given formula
    A_2 = np.sqrt(a_2**2 + b_2**2)
    
    # Compute the bar strength
    bar_strength = A_2 / A_0

    # Compute the phase angle phi_2
    phi_2 = 0.5 * np.degrees(np.arctan2(b_2, a_2))

    # Variability measure of phase angle, here using standard deviation
    phi_2_var = np.std(phi_2)

    return bar_strength, phi_2, phi_2_var

def compute_bar_strength_in_regions(densities, theta, radius_2d, dA, dr):
    """
    Compute the bar strength and phase angle variability for each annular region.

    Parameters:
    - densities: 2D array of densities for each cell in the image.
    - theta: 1D array of azimuthal coordinates for each column of the 2D density array (in radians).
    - radius_2d: 2D array of radii values for each cell in the image.
    - dA: Area of each cell.
    - dr: width of each annular region.

    Returns:
    - results: A list of computed bar strength, phi_2, and phi_2_var for each annular region.
    """
    
    # Define the maximum radius from the data
    r_max = np.max(radius_2d)
    
    # Number of annular regions
    num_regions = int(r_max / dr)

    # Expand theta into a 2D array matching the shape of densities for multiplication
    theta_2d = np.tile(theta, (densities.shape[0], 1))

    results = []

    for i in range(num_regions):
        r_inner = i * dr
        r_outer = r_inner + dr

        # Mask the data for the current annular region
        mask = (radius_2d >= r_inner) & (radius_2d < r_outer)
        masked_densities = densities[mask]
        masked_theta = theta_2d[mask]

        # Compute the mass-equivalent for each cell
        mass_equivalent = masked_densities * dA
        
        # Compute a_2 and b_2 coefficients for m=2 mode
        a_2 = np.sum(mass_equivalent * np.cos(2 * masked_theta))
        b_2 = np.sum(mass_equivalent * np.sin(2 * masked_theta))
        
        # Compute A_0, which is just the sum of all the mass-equivalents
        A_0 = np.sum(mass_equivalent)
        
        # Compute A_2 using the given formula
        A_2 = np.sqrt(a_2**2 + b_2**2)
        
        # Compute the bar strength
        bar_strength = A_2 / A_0

        # Compute the phase angle phi_2
        phi_2 = 0.5 * np.degrees(np.arctan2(b_2, a_2).value)

        # Variability measure of phase angle, here using standard deviation
        phi_2_var = np.std(phi_2)
        results.append((bar_strength, phi_2, phi_2_var))
        
    return results


if __name__ == "__main__":
    # List of datasets and width_pc values
    ds_list = [ 
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0130/DD0130", # 0.02
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0147/DD0147", # 0.19
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0155/DD0155", # 0.27
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0167/DD0167", # 0.39 (clear arms)
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0175/DD0175", # 0.47 (clear bar + arms)
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0189/DD0189", # 0.6 (max clumpy)
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0197/DD0197", # 0.68 (clear bar)
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0217/DD0217", # 0.87 (clear bar)
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/DD0231/DD0231", # 1.0 (clear bar)

        # "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0250/DD0250",
        # "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0240/DD0240",
        # "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/DD0279/DD0279",
        # "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0197/DD0197"
        ]
    width_pc_list = [0.05, 0.1, 0.2, 0.3] # fill this with corresponding width_pc values
    width_pc_list = [0.1]

    with open("bar_strengths_and_phase_angle.txt", "w") as f:
        for width_pc in width_pc_list:
            for ds_path in ds_list:
                # Load dataset
                ds = yt.load(ds_path)
                sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
                sim_label = sim_label.replace("-2", "")
                sim_label = sim_label.replace("RS", "")
                print("Simulation: " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))

                # Grab bh properties and define center, width, and resolution of sliceplots
                ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
                center = ss_pos
                npixels = 2048

                # Obtain angular momentum vector from a small disk and define a larger disk for plotting
                disc_r_pc = disc_h_pc = 0.01
                dx = ds.index.get_smallest_dx().in_units('cm')
                _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
                vecs = ortho_find(L)
                dir = vecs[0]
                north = vecs[1]
                disc_r_pc_big = disc_h_pc_big = 0.6
                disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)
                density, radius_pc = field_from_sliceplot("density", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=True)
                surface_density = density * dx # g/cm^2

                # Compute bar strength
                theta = get_phi_values(surface_density)
                dA = dx**2
                #bar_strength, phi_2, phi_2_var = compute_bar_strength(density, theta, dA)
                bar_strength, phi_2, phi_2_var = compute_bar_strength_in_regions(density, theta, radius_pc, dA, dr=0.001)
                age = ds.current_time.to('Myr') - 124.76*yt.units.Myr
                print("Simulation: {}, DD: {}, age: {:.2f}, slice width: {} pc, Bar Strength: {:.4f}, Phi: {:.4f} deg, Phi_var: {:.5f} deg".format(sim_label, extract_dd_segment(ds.directory), age, width_pc, bar_strength, phi_2, phi_2_var))
                f.write("Simulation: {}, DD: {}, age: {:.2f}, slice width: {} pc, Bar Strength: {:.4f}, Phi: {:.4f} deg, Phi_var: {:.5f} deg".format(sim_label, extract_dd_segment(ds.directory), age, width_pc, bar_strength, phi_2, phi_2_var))