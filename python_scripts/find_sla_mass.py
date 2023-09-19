import numpy as np

def M_sla(dx, dt, v, c, ncells, G=6.67e-20): 
        # dx: cell size in km e.g. dx=9.6e10 for dx=0.003 pc
        # dt: timestep in yr
        # v: relative velocity in km/s
        # c: sound speed in km/s
        # ncells: number of cells in the accretion region
        # G: gravitational constant in km^3 kg^-1 yr^-2
     V=ncells*(dx)**3
     msun_in_kg=1.989e30
     return (((V*(v**2 + c**2)**1.5)/(G**2 * dt))**0.5)/msun_in_kg

print(M_sla(9.6e10, 1.9e9, 1.71, 0.47, 4))