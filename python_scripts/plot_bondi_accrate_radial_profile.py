import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as CM
from yt.units import pc, Msun, s, g, K, cm # Import units
from yt.units.physical_constants import G
from csv_data_plotting.helper_functions import extract_simulation_name
from helper_functions import ss_properties, setup_plot_env
from yt.units.yt_array import YTQuantity

# Constants
proton_mass = 1.67e-24*g  # in grams
thomson_cross_section = 6.65e-25*cm**2  # in cm^2
speed_of_light = 3e10*cm/s  # in cm/s
radiative_efficiency = 0.1  # typically assumed for black holes

# Conversion factors
solar_mass_to_g = YTQuantity(1, "Msun").in_units("g").v
yr_to_s = YTQuantity(1, "year").in_units("s").v
mean_molecular_weight = 1.22  # Adjust based on the gas in your simulation
hydrogen_mass = 1.67e-24  # in grams


def bondi_rate(ds, center, radius, thickness, M_BH, edd_frac=False):
    """
    Calculate the Bondi accretion rate through a spherical shell at a given radius.

    Parameters:
    ds : yt dataset
        The simulation dataset.
    center : array-like
        The center of the spherical shell (e.g., black hole position).
    radius : float
        The radius of the spherical shell (in pc).
    thickness : float
        The thickness of the shell to approximate the 2D surface (in pc).
    M_BH : yt.Quantity
        The mass of the black hole in solar masses (Msun).

    Returns:
    bondi_rate : yt.Quantity
        The Bondi accretion rate in solar masses per year (Msun/yr).
    """
    
    # Create a spherical shell with a small thickness
    shell = ds.sphere(center, radius + thickness) - ds.sphere(center, radius)

    # Get the gas  fields
    density = shell[("gas", "density")].in_units("g/cm**3")  # Gas density in g/cm^3
    sound_speed = shell[("gas", "sound_speed")].in_units("cm/s")  # Sound speed in cm/s
    gas_velocity = shell[("gas", "velocity_magnitude")].in_units("cm/s")  # Gas velocity magnitude in cm/s

    # Bondi rate: 4 * pi * G^2 * M_BH^2 * rho / c_s^3
    bondi_rate_per_cell = (4 * np.pi * G **2 * (M_BH**2) * density / (sound_speed**2 + gas_velocity**2)**(3/2))

    # Find the mean Bondi rate over the shell
    mean_bondi_rate = np.mean(bondi_rate_per_cell).to("Msun/yr")

    if edd_frac:
        edd_rate = eddington_accretion_rate(M_BH)
        mean_bondi_rate /= edd_rate
    
    return mean_bondi_rate

def eddington_accretion_rate(bh_mass):
    """
    Compute the Eddington accretion rate for a given black hole mass.
    """
    bh_mass_g = bh_mass.to('g')
    edd_rate_g_s = (4 * np.pi * G * bh_mass_g * proton_mass) / (radiative_efficiency * thomson_cross_section * speed_of_light)
    return edd_rate_g_s.to('Msun/yr')

def bondi_eddington_fraction(bondi_rate, bh_mass):
    """
    Compute the Bondi accretion rate Eddington fraction: bondi_rate / edd_rate.
    """
    edd_rate = eddington_accretion_rate(bh_mass)
    return bondi_rate / edd_rate


def mass_inflow_rate(ds, center, radius, thickness, normalize_by=None):
    """
    Calculate the mass inflow rate through a spherical shell at a given radius.

    Parameters:
    ds : yt dataset
        The simulation dataset.
    center : array-like
        The center of the spherical shell (e.g., black hole position).
    radius : float
        The radius of the spherical shell (in pc).
    thickness : float
        The thickness of the shell to approximate the 2D surface (in pc).

    Returns:
    inflow_rate : yt.Quantity
        The mass inflow rate in solar masses per year (Msun/yr).
    """
    
    # Create a spherical shell with a small thickness
    shell = ds.sphere(center, radius + thickness) - ds.sphere(center, radius)

    # Get the gas density and radial velocity fields
    density = shell[("gas", "density")]  # in units of g/cm^3
    radial_velocity = shell[("gas", "radial_velocity")].in_units("cm/s")  # positive outward, negative inward

    # Select only inward flowing gas (radial_velocity < 0 for inflow)
    inflow_mask = radial_velocity < 0
    inflowing_density = density[inflow_mask]
    inflowing_radial_velocity = radial_velocity[inflow_mask]

    # Get cell volumes and calculate approximate cross-sectional area for the spherical shell
    cell_volumes = shell["index", "cell_volume"].in_units("cm**3")  # Volume of each cell in cm^3
    cell_radii = shell["index", "dx"].in_units("cm")[inflow_mask]  # Cell widths (dx) for area approximation
    cell_areas = (cell_volumes[inflow_mask] / cell_radii)

    # Mass inflow rate is the sum of inflowing mass flux through the shell
    mass_flux = inflowing_density * inflowing_radial_velocity * cell_areas  # in g/s

    # Sum the mass flux to get the total inflow rate
    inflow_rate = np.sum(mass_flux).to("Msun/yr")

    # Normalize
    if normalize_by == 'area':
        normalized_inflow = inflow_rate / shell_area
    elif normalize_by == 'eddington':
        edd_rate = eddington_accretion_rate(bh_mass)
        normalized_inflow = inflow_rate / edd_rate
    elif normalize_by == 'bh_mass':
        normalized_inflow = inflow_rate / bh_mass
    elif normalize_by == 'radius':
        normalized_inflow = inflow_rate / r.in_units("pc").v  # Normalize by the radius in parsecs
    elif normalize_by is None:
        normalized_inflow = inflow_rate
    else:
        raise ValueError("Normalization method not recognized. Choose 'area', 'eddington', 'bh_mass', or 'radius'.")
    
    return normalized_inflow


def compute_accretion_rate_at_radii(ds, bh_pos, bh_mass, rate_type='bondi', r_step=0.02, edd_frac=False, normalization=None):
    """
    Compute accretion rate at each radial distance for a given dataset and black hole properties.
    """

    rates = []
    for r in radial_bins:
        # Bondi rate
        if rate_type == 'bondi':
            rate = bondi_rate(ds, bh_pos, r, r_step, bh_mass, edd_frac=edd_frac)
            
        # Mass flux rate
        elif rate_type == 'inflow':
            rate = mass_inflow_rate(ds, bh_pos, r, r_step, normalize_by=normalization)

        rates.append(rate)
        print(f"Radius = {r:.2f}, Rate = {rate:.2e}")

    return rates

# Set up the figure
plt.figure(figsize=(8, 6))
setup_plot_env(fontsize=14)

# Event 2: initial state at 31.70 Myr, take fb values at 31.90 Myr
simulations = [
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0445/DD0445", "name": "Initial State"},
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2783/DD2783", "name": "Min Feedback"},
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0641/DD0641", "name": "Mid Feedback"},
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695", "name": "Max Feedback"}
]

# Event 1: initial state at 2.50 Myr, take fb values at 2.70 Myr
simulations = [
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0153/DD0153", "name": "Initial State"},
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001/DD0341/DD0341", "name": "Min Feedback"}, # temp
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0352/DD0352", "name": "Mid Feedback"},
    {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0352/DD0352", "name": "Max Feedback"}
]

# Loop over each simulation to compute and plot the accretion rate
rate_type = 'bondi'  # 'bondi' or 'inflow'
edd_frac = False
normal = None  # 'area', 'eddington', 'bh_mass', 'radius'
event = 1

# Parameters for radial distances and increment
r_min = 0.005 * pc
r_max = 5 * pc
r_step = 0.05 * pc
radial_bins = np.arange(r_min.value, r_max.value + r_step.value, r_step.value) * pc

# Generate a list of colors from colormap
n_colors = len(simulations)  # Number of colors you want
cubehelix = CM.get_cmap('cividis', n_colors)
colors = [cubehelix(i) for i in range(n_colors)]
colors = ['grey', 'gold', 'orange', 'crimson']

for i, sim_info in enumerate(simulations):
    print(f"Processing simulation: {sim_info['name']}")
    
    ds = yt.load(sim_info["path"])
    sim_name = extract_simulation_name(sim_info["path"])
    sim_name = sim_info['name'] if 'name' in sim_info else sim_name
    
    bh_pos, bh_mass, bh_age = ss_properties(ds)
    print(f"Black hole mass: {bh_mass:.2e} ")
    print(f"Black hole age: {bh_age[0]/1e6:.2e} Myr")
    
    accretion_rates = compute_accretion_rate_at_radii(ds, bh_pos, bh_mass.to('g'), rate_type=rate_type, r_step=r_step, edd_frac=edd_frac, normalization=normal)
    
    # Plot the accretion rate for this simulation
    plt.plot(radial_bins, accretion_rates, marker='o', color=colors[i], label=sim_name, alpha=0.6)

# Finalize the plot
ylabel = r'$\dot{M}$ ($\mathrm{M_\odot}$/yr) \,' + (normal if normal else '')
ylabel = r'$\dot{M}/\dot{M}_{\rm Edd}$' if edd_frac else ylabel
plt.xlabel('Radius (pc)')
plt.ylabel(ylabel)
plt.axhline(0, color='black', linewidth=0.5)
#plt.yscale('symlog', linthresh=1, linscale=0.5)
plt.yscale('log')
plt.title(f'Event {event}: {rate_type.capitalize()} Rate Comparison')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Save the combined plot
plot_name = f"plots/combined_accretion_profile_{rate_type}_event_{event}_{float(r_max.d):.0f}pc.png"
plt.savefig(plot_name)
print(f"Saved combined plot: {plot_name}")
plt.show()
