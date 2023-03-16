import yt
import numpy as np

G = 6.67e-8 * (yt.units.cm ** 3)/(yt.units.g*yt.units.s**2) # cgs
M = 270*1.989e+33 * yt.units.g # cgs

def _keplerian_frequency_BH(field, data):
    return np.sqrt(G*M / (data[("index", "radius")].to('cm')**3))


# add epicyclic frequency kappa / angular frequency omega derived field
def _angular_frequency_keplerian(field, data):
    return np.abs(data["gas", "radial_velocity"].to('cm/s')) / (2*np.pi*data["index", "radius"].to('cm')) #cm/s / cm = 1/s


def _angular_frequency(field, data):
    omega = np.abs(data["gas", "velocity_cylindrical_theta"].to('cm/s')) / (data["index", "radius"].to('cm'))  # cm/s / cm
    return omega

def _angular_frequency_square(field, data):
    omega = np.abs(data["gas", "velocity_cylindrical_theta"]) / (2 * np.pi * data["index", "radius"])  # cm/s / cm
    return omega**2


# def _angular_frequency_square(field, data):
#     omega = data["gas", "tangential_velocity"].to('cm/s') / data["index", "radius"].to('cm')  # cm/s / cm
#     return omega**2


# For a Keplerian disk, \kappa =\Omega .
def _epicyclic_frequency_ratio(field, data):
    omega = np.abs(data["gas", "velocity_cylindrical_theta"])/(2*np.pi*data["index", "radius"]) #cm/s / cm
    return omega/data["gas", "angular_frequency_keplerian"]


# def _epicycle_frequency_k(field, data):
#     slope = -2.006220329419717 / (yt.units.cm * yt.units.s ** 2)
#     k2 = (4 * data["gas", "angular_frequency_square"]) + (data["index", "radius"].to('cm') * slope)
#     return np.sqrt(np.abs(k2))

def _epicycle_frequency_k(field, data):
    slope = -1.1955583144810638 / (yt.units.cm * yt.units.s)
    Omega = data["gas", "angular_frequency"].to('1/s')
    #Omega = np.abs(data["gas", "velocity_cylindrical_theta"].to('cm/s'))
    R = data["index", "radius"].to('cm')
    k2 = (4 * Omega**2) * (1 + ((0.5 * R * slope)/Omega))
    return np.sqrt(np.abs(k2))


# simulation cooling rate [erg g^-1 s^-1]
def _simulation_cooling_rate(field, data):
    # specific energy / cooling time
    return data['enzo', 'GasEnergy'].to('ergs/g') / data['enzo', 'Cooling_Time'].to('s')


# radiative cooling rate [erg s-1 cm-3]
def _radiative_cooling_rate(field, data):
    # specific energy * mass density / cooling time
    return data["enzo", "simulation_cooling_rate"] * data["gas", "density"]


# advective cooling rate [erg s-1 cm-3]
def _advective_cooling_rate(field, data):
    kb = 1.3807e-16 * ((yt.units.cm ** 2 * yt.units.g) / (yt.units.s ** 2 * yt.units.K))
    xi = -0.6
    m_p = 1.6726e-24*yt.units.g
    num = (data["gas", "density"].to('g/cm**3') * data["gas", "radial_velocity"].to('cm/s') *
           data['gas', 'temperature'] * kb * xi)
    denom = (data["index", "radius"].to('cm') * m_p)
    return (num/denom).to('erg/((cm**3)*s)')


def _height(field, data):
    return np.abs(data["index", "cylindrical_z"])


def add_fields_ds(ds):
    ds.add_field(
        name=("gas", "angular_frequency_keplerian"),
        function=_angular_frequency_keplerian,
        sampling_type="local",
        units="1/s",
    )

    ds.add_field(
        name=("gas", "epicyclic_frequency_ratio"),
        function=_epicyclic_frequency_ratio,
        sampling_type="local",
        units="dimensionless",
    )

    ds.add_field(
        name=("gas", "angular_frequency"),
        function=_angular_frequency,
        sampling_type="local",
        units="1/s",
    )

    ds.add_field(
        name=("gas", "keplerian_frequency_BH"),
        function=_keplerian_frequency_BH,
        sampling_type="local",
        units="1/s",
    )

    ds.add_field(
        name=("enzo", "simulation_cooling_rate"),
        function=_simulation_cooling_rate,
        sampling_type="local",
        units="erg/s/g"
    )

    ds.add_field(
        name=("enzo", "radiative_cooling_rate"),
        function=_radiative_cooling_rate,
        sampling_type="local",
        units="erg/(s*cm**3)"
    )

    ds.add_field(
        name=("enzo", "advective_cooling_rate"),
        function=_advective_cooling_rate,
        sampling_type="local",
        units="erg/(s*cm**3)"
    )

    ds.add_field(
        ("index", "height"),
        function=_height,
        sampling_type="local",
        units="cm"
    )

    ds.add_field(
        ("gas", "angular_frequency_square"),
        function=_angular_frequency_square,
        sampling_type="local",
        units="1/s**2"
    )

    ds.add_field(
        ("gas", "epicycle_frequency_k"),
        function=_epicycle_frequency_k,
        sampling_type="local",
        units="1/s"
    )

