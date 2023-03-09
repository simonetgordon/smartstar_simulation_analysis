"""
For producing radial profiles of advective vs radiative cooling in nuclear disc
python plot_cooling_rates.py DD0130/DD0130
"""

import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
import matplotlib.pyplot as plt
from derived_fields import add_fields_ds

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))
add_fields_ds(ds)
ad = ds.all_data()

q_time = ad['enzo', 'Cooling_Time'] # [s]
print("q_time: ", q_time)
q_energy = ad['enzo', 'GasEnergy'].to('ergs/g') # g cm^2 s^-2 cm^-3 = g cm^-1 s^-2 = ergs g^-1
print("q_energy: ", q_energy)


print("cooling rate: ", ad["enzo", "cooling_rate"])

q_sim = ad["enzo", "cooling_rate"] # [erg g^-1 s^-1] = (g cm^2 s^-2) g^-1 s^-1 = cm^2 s^-3
n = ad["gas", "H_nuclei_density"] # cm^-3
rho = ad["gas", "density"] # g/cm^-3
v_r = ad["gas", "radial_velocity"].to('cm/s')
temp = ad["gas", "temperature"]
r = ad["index", "radius"].to('cm')

# constants
kb = 1.3807e-16*((yt.units.cm**2 * yt.units.g)/(yt.units.s**2 * yt.units.K)) # cm2 g s-2 K-1
m_p = 1.6726e-24*yt.units.g # g
xi = -0.65 # dimensionless parameter following Chen et al. (1995)
# fh2 = ? # molecular fraction - may not need this
print("kb ", kb)

# radiative cooling rate [erg s-1 cm-2]
#q_rad = (q_sim * n * sigma * fh2**2) / m_p

q_rad = q_sim * rho # [erg s^-1 cm^-3]
print("q_rad: ", q_rad)

# advective cooling rate [erg s-1 cm-3]
q_adv = ((rho * v_r * temp * kb * xi) / (r * m_p)).to('erg/((cm**3)*s)')
print("q_adv: ", q_adv)