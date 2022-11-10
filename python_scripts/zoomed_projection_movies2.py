import sys
import yt
yt.enable_parallelism()
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
import ytree
import numpy as np
import matplotlib.pyplot as plt

## Make Sphere of Most Massive Halo ##

ds_10_3 = yt.load("DD0148/DD0148") # t = 149.3, 20 Myr forward

# Load merger tree of dataset (up to DD0118 in gas run)
a = ytree.load('../gas+dm-L3/rockstar_halos/out_0.list')

# Load my_tree and find radius
a1 = ytree.load('../gas+dm-L3/tree_810/tree_810.h5')

r_halo = a1[0]["virial_radius"].to('pc')
r = ds_10_3.quan(r_halo.d, "pc") # virial radius in pc

if __name__ == "__main__":
    es = yt.load_simulation(sys.argv[1], "Enzo", find_outputs=True)
    es.get_time_series()

    nested_level = es.parameters["CosmologySimulationNumberOfInitialGrids"] - 1
    if nested_level == 0:
        region = None
        center = [0.49053251, 0.4946698,  0.50963437]
        width = None
    else:
        left_edge = es.parameters["CosmologySimulationGridLeftEdge[%d]" % nested_level]
        right_edge = es.parameters["CosmologySimulationGridRightEdge[%d]" % nested_level]
        my_width =  (right_edge - left_edge).max()
        center = [0.49053251, 0.4946698,  0.50963437]


    for ds in es.piter():

        ds.all_data()

        # Add Pop III metallicity fields
        add_p2p_fields(ds)
        add_p2p_particle_filters(ds)

        sp = ds.sphere(center, 5*r)
        pos_star = sp['pop3', 'particle_position'].to('unitary')
        center_star = pos_star.d


        # Metallicity
        p = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width=(200, 'pc'),center= center_star, weight_field='density')
        p.set_cmap('metallicity3', 'kamae')
        p.set_axes_unit('pc')
        p.annotate_timestamp(corner='lower_right')
        p.annotate_scale(corner='lower_left')
        p.save("frames_zoom/")

