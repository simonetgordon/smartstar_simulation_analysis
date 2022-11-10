import sys
import yt
yt.enable_parallelism()
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
from yt.units import pc

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

        radius = 2172*pc
        r = ds.quan(radius.d, 'unitary')
        sp = ds.sphere(center, r)
        pos_star = sp['pop3', 'particle_position'].to('unitary')
        center_star = pos_star.d
        
      
        # Metallicity
        p = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width=(200, 'pc'),center= center_star, weight_field='density')
        p.set_cmap('metallicity3', 'kamae')
      

        p.set_axes_unit('pc')
        p.annotate_timestamp(corner='lower_right')
        p.annotate_scale(corner='lower_left')
        p.save("frames_zoom/")

        # Temperature
        p1 = yt.ProjectionPlot(ds, "x", ("gas", "temperature"), width= (200, 'pc'),center=center_star, weight_field='density')
        p1.set_cmap('temperature', 'RED TEMPERATURE')
     

        p1.set_axes_unit('pc')
        p1.annotate_timestamp(corner='lower_right')
        p1.annotate_scale(corner='lower_left')
        p1.save("frames_zoom/")

        # Density
        p2 = yt.ProjectionPlot(ds, "x", ("gas", "density"), width=(200, 'pc'),center=center_star, weight_field='density')
        p2.set_cmap('density', 'Rainbow')

        p2.set_axes_unit('pc')
        p2.annotate_timestamp(corner='lower_right')
        p2.annotate_scale(corner='lower_left')
        p2.save("frames_zoom/")

        # H2 fraction
        p3 = yt.ProjectionPlot(ds, "x", ("gas", "H2_p0_fraction"), width=(200,'pc'),center=center_star, weight_field='cell_mass')
        p3.set_cmap("H2_p0_fraction", "kelp")

        
        p3.set_axes_unit('pc')
        p3.annotate_timestamp(corner='lower_right')
        p3.annotate_scale(corner='lower_left')
        p3.save("frames_zoom/")
        
