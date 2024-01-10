import numpy as np
import matplotlib.cm as cm
import yt
from yt.utilities.math_utils import ortho_find
from matplotlib import rc, pyplot
import matplotlib.pyplot as plt


def configure_font(fontsize=14):
    pyplot.rcParams['font.size'] = fontsize
    pyplot.rcParams['font.weight'] = 'light'
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"

def ss_properties(ds):
    print("ds = ", ds)
    ad = ds.all_data()
    try:
        # find ss properties
        ss_creation = ad['SmartStar', 'creation_time'].to('yr')
        ss_pos = ad['SmartStar', 'particle_position'].to('unitary')[0]
        ss_mass = ad['SmartStar', 'particle_mass'].to('Msun')[0]

        # find ss age
        time = ds.current_time.to('yr').d
        creation = ss_creation.d  # strip units off variables
        ss_age = time - creation
    except:
        ss_pos = None
        ss_mass=0*yt.units.msun
        ss_age = [0*yt.units.Myr]

    return ss_pos, ss_mass, ss_age


# make disc data container
def _make_disk_L(ds, center, width, height):
    """
    Make a disk data container and calculate the angular momentum vector
    ds: yt dataset
    center: center of the disk
    width: width of the disk (2*radius)
    height: height of the disk
    """
    if isinstance(width, float):
        width = width*yt.units.pc
    if isinstance(height, float):
        height = height*yt.units.pc
    sp = ds.sphere(center, width)
    L = sp.quantities.angular_momentum_vector()
    L /= np.sqrt((L ** 2).sum())
    disk = ds.disk(center, L, width, height)
    return disk, L


def tidy_data_labels(labels, custom_name=None):
    if custom_name:
        return custom_name
    
    # Function to apply replacements to a single label
    def apply_replacements(label):
        label = label.replace("-2", "")
        label = label.replace("RS", "")
        label = label.replace("-4dx", "")
        label = label.replace("/2B.b08", "") # for plot_disc_attributes.py 2B.b08
        try: 
            parts = labels.split('/')
            if len(parts) == 2 and parts[0] == parts[1]:
                label = parts[0]
        except AttributeError:
            pass
        return label

    # Check if labels is a list or a single label
    if isinstance(labels, list):
        # Process each label in the list
        data_labels = [apply_replacements(label) for label in labels]
    elif isinstance(labels, str):
        # Process a single label
        data_labels = apply_replacements(labels)
    else:
        raise TypeError("labels should be a string or a list of strings")

    return data_labels


def M_sla(dx, dt, v, c, ncells, G=6.67e-20):
    """
    Calculate SLA Mass as in Beckmann et al. 2018b
    Input:
    dx: cell size in pc
    dt: time step in Myr
    v: velocity dispersion in km/s
    c: sound speed in km/s
    ncells: number of cells in the disc
    G: gravitational constant in km^3 kg^-1 s^-2
    e.g. print(M_sla(9.6e10, 1.9e9, 1.71, 0.47, 4))
    """ 
    V=ncells*(dx)**3
    msun_in_kg=1.989e30
    return (((V*(v**2 + c**2)**1.5)/(G**2 * dt))**0.5)/msun_in_kg

def adaptive_moving_average(data, window_size=5):
    """
    Compute an adaptive moving average on the logarithm of the data.
    used in:
    - fourier_modes_2.ipynb
    - 
    
    :param data: The input data (list or array).
    :param window_size: The window size for the moving average.
    :return: An array of the moving average values in the original data space.
    """
    # Take the logarithm of data, excluding non-positive values
    log_data = np.log(data[data > 0])
    data_length = len(log_data)
    log_moving_avg = np.zeros(data_length)

    for i in range(data_length):
        start = max(i - window_size // 2, 0)
        end = min(i + window_size // 2 + 1, data_length)
        log_moving_avg[i] = np.mean(log_data[start:end])

    # Exponentiate to return to original data space
    moving_avg = np.exp(log_moving_avg)

    # Handle edge cases if original data had non-positive values
    moving_avg_full = np.full_like(data, np.nan)
    positive_indices = np.where(data > 0)[0]
    moving_avg_full[positive_indices] = moving_avg

    return moving_avg_full


def extract_colors(cmap_name, n, portion=None, start=None, end=None):
    """
    Extract a list of colors from a matplotlib colormap.
    :param cmap_name: The name of the colormap.
    :param n: The number of colors to extract.
    :param portion: The portion of the colormap to extract from. Can be 'beginning', 'middle', 'end', or None.
    :param start: The starting point of the portion to extract from. Must be between 0 and 1.
    :param end: The ending point of the portion to extract from. Must be between 0 and 1.
    :return: A list of RGB colors.
    """
    cmap = cm.get_cmap(cmap_name)

    if start is not None and end is not None:
        values = np.linspace(start, end, n)
    elif portion == "beginning":
        values = np.linspace(0, 0.3, n)
    elif portion == "middle":
        values = np.linspace(0.3, 0.95, n)
    elif portion == "end":
        values = np.linspace(0.7, 1, n)
    elif portion is None:
        values = np.linspace(0, 1, n)
    else:
        raise ValueError("Invalid portion specified. Use 'beginning', 'middle', 'end', or None.")

    colors = cmap(values)
    return colors


def field_from_sliceplot(field, ds, disk, center, width_pc, north, dir, npixels=2048, radius=False):
    """
    Compute field from a slice plot of a dataset.
    Surface Density = slice plot density * cell height
    """
    p = yt.SlicePlot(ds, dir, ("gas", field), center=disk.center, width=(width_pc, "pc"), data_source=disk)
    slc_frb = p.data_source.to_frb((width_pc, "pc"), npixels)
    slc_field = slc_frb[("gas", field)]
    if radius:
        radius = slc_frb['radius'].to('pc')
        return slc_field, radius
    else:
        return slc_field
    

def ToomreQ(cs, kappa, G, surface_density):
    """
    Calculate the Toomre Q parameter for linear stability
    """
    Q = cs * kappa / (np.pi * G * surface_density*0.6*1.67e-24)
    return Q


def kappa2D(frb):
    """
    Calculate the epicyclic frequency kappa = v/r
    """
    kappa = np.abs(frb['velocity_cylindrical_theta']) / frb['radius']
    return kappa


def kappa2D_vmag(frb):
    """
    Calculate the epicyclic frequency kappa = v/r
    """
    kappa = frb['velocity_magnitude'] / frb['radius']
    return kappa


def orbital_velocity(ds, disk):
    G = 6.67e-8 * (yt.units.cm ** 3)/(yt.units.g*yt.units.s**2) # cgs
    return np.sqrt(G * ds.r['SmartStar', 'particle_mass'].to('g') / disk['index', 'radius'].to('cm'))
    

def critical_density(ds, n_crit=False):
    # Get the cosmological parameters from the dataset parameters
    hubble_constant = yt.YTQuantity(ds.parameters.get("CosmologyHubbleConstantNow")*100, "km/s/Mpc")
    matter_density = yt.YTQuantity(ds.parameters.get("CosmologyOmegaMatterNow"), "")
    cosmological_constant_density = yt.YTQuantity(ds.parameters.get("CosmologyOmegaLambdaNow"), "")

    # Convert the Hubble constant to CGS units (cm/s/Mpc)
    hubble_constant_cgs = (hubble_constant).to("cm/s/Mpc")

    # Conver Hubble constant to 1/s units
    hubble_constant_pers = hubble_constant_cgs / yt.YTQuantity(3.1e19*1e5, "cm/Mpc") 

    # Calculate the critical density
    critical_density = (3 * hubble_constant_pers**2) / (8 * np.pi * yt.physical_constants.G) * (matter_density + cosmological_constant_density)

    # Convert the critical density to the desired units (1/cm^3)
    critical_density = critical_density.to("g/cm**3")

    # Convert to hydrogen nuclei density if desired
    if n_crit:
        critical_density /= yt.physical_constants.mh # /cm**3

    return critical_density


def format_sci_notation(x):
    a, b = '{:.2e}'.format(x).split('e')
    return r'$\rm {} \times 10^{{{}}}$'.format(a, int(b))


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
    
def extract_simulation_name(fp):
    """
    Extract the simulation name from a file path.

    Parameters:
    fp (str): The file path.

    Returns:
    str: The extracted simulation name.
    """
    # Find all substrings that match the pattern
    matches = re.findall(r'/([^/]+)/DD', fp)

    # Return the last match (closest to the end of the string)
    if matches:
        return matches[-1]
    else:
        print("No match found")
        return None

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


# Next few functions pertain to plot_density_slices_toomre.py
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
    """
    Calculate the theta array for a given set of dimensions.
    Input:
        dimensions: tuple of dimensions for the array
    Output:
        theta_array: 2D array of theta values for each pixel in the array
    """
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


def find_fourier_modes_and_phase_angles(radii, radius_pc, densities, theta, dV, dr=0.001, cylindrical_theta_velocity=None, phi2_only=False):
    """
    Compute the m=1, 2 Fourier mode strengths and phase angles for a given annular region of a disk.
    Input:
        radii: list of radii for annular regions
        radius_pc: 2D array of radii values for each cell in the disk
        densities: 2D array of densities values for each cell in the disk
        theta: 2D array of theta values for each cell in the disk
        dV: volume of each cell in the disk
        dr: thickness of annular regions
    Output:
        m1_strengths: list of m=1 Fourier mode strengths for each annular region
        m2_strengths: list of m=2 Fourier mode strengths for each annular region
        phi_1_values: list of m=1 Fourier mode phase angles for each annular region
        phi_2_values: list of m=2 Fourier mode phase angles for each annular region
    """
    m1_strengths = []
    m2_strengths = []
    phi_1_values = []
    phi_2_values = []
    angular_speeds = []
    for r in radii:
        # Define mask for annular region
        mask = (radius_pc >= r) & (radius_pc < r + dr)

        # Get densities and thetas for this annular region, where 0 < theta < 2pi - all same shape (713,)
        masked_densities = densities[mask]
        if not theta.shape == densities.shape:
            theta = np.tile(theta, (densities.shape[0], 1))
        masked_theta = theta[mask]
        if cylindrical_theta_velocity is not None:
            masked_theta_vel = cylindrical_theta_velocity[mask]
            angular_speed = masked_theta_vel / (r*yt.units.pc).to("cm")
            angular_speeds.append(angular_speed) # cm/s / cm = 1/s

        # Compute the mass-equivalent for each cell in this region
        mass_equivalent = masked_densities * dV

        # Compute a and b coefficients for m=2 and m=1 modes
        a_2 = np.sum(mass_equivalent * np.cos(2 * masked_theta))
        b_2 = np.sum(mass_equivalent * np.sin(2 * masked_theta))
        a_1 = np.sum(mass_equivalent * np.cos(masked_theta))
        b_1 = np.sum(mass_equivalent * np.sin(masked_theta))

        if phi2_only:
            phi_2 = 0.5 * np.degrees(np.arctan2(b_2, a_2).value)
            phi_2_values.append(phi_2)
            continue

        # Compute A_0 for this region
        A_0 = np.sum(mass_equivalent)

        # Compute Fourier modes for this region
        A_2 = np.sqrt(a_2**2 + b_2**2)
        A_1 = np.sqrt(a_1**2 + b_1**2)

        # Normalise the Fourier mode strengths for this region
        m2_strength = A_2 / A_0 if A_0 != 0 else 0
        m1_strength = A_1 / A_0 if A_0 != 0 else 0
        m2_strengths.append(m2_strength)
        m1_strengths.append(m1_strength)

        # Compute the phase angle phi_2 for this region
        phi_2 = 0.5 * np.degrees(np.arctan2(b_2, a_2).value) # -90 < phi_2 < 90
        phi_1 = np.degrees(np.arctan2(b_1, a_1).value) # -180 < phi_1 < 180
        phi_2_values.append(phi_2)
        phi_1_values.append(phi_1)
    
    if cylindrical_theta_velocity is not None:
        return m1_strengths, m2_strengths, phi_1_values, phi_2_values, angular_speeds
    
    if phi2_only:
        return phi_2_values
    
    return m1_strengths, m2_strengths, phi_1_values, phi_2_values


def find_previous_ds(ds, root_dir, sim):
    """
    Find the dataset corresponding to the previous snapshot.
    Input:
        ds: current dataset
        root_dir: root directory of simulation
        sim: simulation name
    Output:
        ds_previous: dataset corresponding to previous snapshot
    """
    num_previous = str(int(extract_dd_segment(ds.directory).replace("DD0", "")) - 1)
    dd_previous = "DD0" +  num_previous + "/DD0" + num_previous
    return yt.load(os.path.join(root_dir[0], sim[0], dd_previous))


def find_pattern_speed_rad_per_sec(ds, root_dir, sim, phi_2_values, ss_age, width_pc, disc_h_pc, disc_r_pc, disc_r_pc_big, disc_h_pc_big, npixels, radii, theta, dV):
    """
    Find the pattern speed of the bar in rad/s.
    Method: 
        - compute the difference between phi2 values in current and previous snapshot
        - divide by the difference in time between snapshots
    """
    # Find the previous dataset phi2 values
    ds_previous = find_previous_ds(ds, root_dir, sim)
    center, _, ss_age_previous = ss_properties(ds_previous)
    _, L = _make_disk_L(ds_previous, center, disc_r_pc, disc_h_pc)
    vecs = ortho_find(L)
    dir = vecs[0]
    north = vecs[1]
    disk = ds_previous.disk(center, L, disc_r_pc_big, disc_h_pc_big)
    density, radius_pc = field_from_sliceplot("density", ds_previous, disk, center, width_pc, north, dir, npixels=npixels, radius=True)
    phi_2_values_previous = find_fourier_modes_and_phase_angles(radii, radius_pc, density, theta, dV, phi2_only=True)

    # Compute the difference between phi2 values at each snapshot
    delta_phi = np.array(phi_2_values) - np.array(phi_2_values_previous)
    delta_t = np.array(ss_age[0]) - np.array(ss_age_previous[0])
    pattern_speeds = delta_phi / delta_t

    # Convert from deg/Myr to rad/s
    conversion_factor = np.pi / (180 * 3.154e7) # from deg/yr to rad/s
    pattern_speeds_radians_per_sec = pattern_speeds * conversion_factor

    return pattern_speeds_radians_per_sec


def compute_radial_profile(radius, data, num_bins=128):
    """
    Compute the radial profile of surface density.

    Parameters:
    radius (2D array): frb image map of radius values.
    data (2D array): frb image map of surface density/toomre Q/etc. values.
    num_bins (int): Number of bins to use for the radial profile.

    Returns:
    bin_centers (1D array): Center of each radial bin.
    profile (1D array): Average surface density within each bin.
    """
    # Flatten the arrays
    radius_flat = np.array(radius).flatten()
    data_flat = np.array(data).flatten()
    
    # Define the radial bins in log space
    radial_bins = np.logspace(np.log10(radius_flat.min()+1e-5), np.log10(radius_flat.max()), num_bins+1)
    
    # Compute the bin centers
    bin_centers = np.sqrt(radial_bins[:-1] * radial_bins[1:])  # geometric mean for log bins
    
    # Compute the profile by averaging surface density within each bin
    profile = np.array([data_flat[(radius_flat >= rb_start) & (radius_flat < rb_end)].mean() for 
                        rb_start, rb_end in zip(radial_bins[:-1], radial_bins[1:])])

    return bin_centers, profile



def radial_profile(field, disk, n_bins, cell_width_pc):
    bins = np.logspace(np.log10(cell_width_pc), np.log10(disk["index", "radius"].to('pc').max()), n_bins+1)
    counts_r, r_bin_edges = np.histogram(disk["index", "radius"].to('pc'), bins=bins)
    y, radius = np.histogram(disk["index", "radius"].to('pc'), weights=field, bins=bins)
    y = np.nan_to_num(y)
    y = y/counts_r
    y = np.nan_to_num(y)
    y_no_zeros = y[y > 0]
    r_no_zeros = radius[:-1][y > 0]
    return r_no_zeros, y_no_zeros


def make_frb(ds, L, center, width= 10*yt.units.pc, npixels=1024, height=0.05*yt.units.pc):
    """
    Make a fixed resolution buffer (frb) of the disk from its center and 
    angular momentum vector L.
    """
    cutting = ds.cutting(L, center)
    frb = cutting.to_frb(width, npixels, height=height)
    return frb, height
