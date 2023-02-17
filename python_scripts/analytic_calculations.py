# import yt
# import numpy as np

def scale_radius_cgs(m, v, G):
    return (2*G*m)/(v**2)


def interpolated_scale_radius_cgs(m, v, c, G):
    return (G * m) / (v**2 + c**2)


def eddington_rate(mparticle_msun):
    mparticle_g = mparticle_msun*SolarMass
    mdot_edd = 4.0 * PI * GravConst * mparticle_g * mh / (clight * eta_disk * sigma_thompson)
    return mdot_edd*yr_s/SolarMass

# eddington rate constants in cgs units
PI = 3.1415927
GravConst = 6.6740831e-8
mh = 1.67262171e-24
clight = 2.99792458e10
sigma_thompson = 6.65E-25
eta_disk = 0.1 # optically thick and e_radiative = 0.11
SolarMass = 1.9891e33
yr_s = 3.1556952E7

# Avg_Density = 514899 cm^-3,
# AverageTemp = 4.099178e+02 K,
# Average cInfinity = 3.427143e-01 km/s,
#  Average vInfinity = 1.981252e+00 km/s,

m_msun = 270
v_kms = 28.46263486
Temp = 4.099178e+02
c_kms = 3.427143e-01

kms_cms = 100000 # cm per second
msun_g = 1.999999999e33 #grams
G_cgs = 6.67e-8 # cm3 g−1 s−2
cm_pc = 3.24078e-19

print("scale radius (cm) = ", scale_radius_cgs(m_msun*msun_g, v_kms*kms_cms, G_cgs))
print("scale radius (pc) = ", scale_radius_cgs(m_msun*msun_g, v_kms*kms_cms, G_cgs)*cm_pc)
print("interpolated scale radius (pc) = ", interpolated_scale_radius_cgs(
    m_msun*msun_g, v_kms*kms_cms, c_kms*kms_cms, G_cgs)*cm_pc)
