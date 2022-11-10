import sys
import yt
yt.enable_parallelism()
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters

if __name__ == "__main__":
    es = yt.load_simulation(sys.argv[1], "Enzo", find_outputs=True)
    es.get_time_series()

    nested_level = es.parameters["CosmologySimulationNumberOfInitialGrids"] - 1
    if nested_level == 0:
        region = None
        center0 = [0.5]*3
        width = None
    else:
        left_edge = es.parameters["CosmologySimulationGridLeftEdge[%d]" % nested_level]
        right_edge = es.parameters["CosmologySimulationGridRightEdge[%d]" % nested_level]
        my_width = (right_edge - left_edge).max()
        center0 = 0.5 * (left_edge + right_edge)

    
    for ds in es.piter():

        add_p2p_fields(ds)
        add_p2p_particle_filters(ds)

        
        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

       
        
        #sp = ds.sphere(center0, width)
        #star_pos = sp['pop3', 'particle_position'].to("unitary")
        #center = star_pos[0].d
        center = center0
        
        # Gas density
        p1 = yt.ProjectionPlot(ds, "x", ("gas", "density"), width=(200, 'pc'), center=center, data_source=region, weight_field='density')
        p1.set_cmap('density', 'turbo')
        p1.set_axes_unit('pc')
        p1.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p1.annotate_scale(corner='lower_left')
        time_array = ds.current_time.to('Myr')
        time = time_array.d
        t_sinceSN = time - 128.1
        p1.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
        p1.save("frames_gas/")
        
        
        # Temperature
        p2 = yt.ProjectionPlot(ds, "x", ('gas', 'temperature'), width=(200,'pc'), center=center, data_source=region, weight_field='density')
        p2.set_cmap("temperature", "RED TEMPERATURE")
        p2.set_axes_unit('pc')
        p2.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p2.annotate_scale(corner='lower_left')
        p2.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
        p2.save("frames_gas/")
        
        
        # Metallicity
        p3 = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width=(200, 'pc'),center= center, data_source=region, weight_field='density')
        p3.set_cmap('metallicity3', 'kamae')

        p3.set_axes_unit('pc')
        p3.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
        p3.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p3.annotate_scale(corner='lower_left')
        p3.save("frames_metal/")
        
        
        # H2 fraction
        p4 = yt.ProjectionPlot(ds, "x", ("gas", "H2_p0_fraction"), width=(200,'pc'),center=center, data_source=region, weight_field='cell_mass')
        p4.set_cmap("H2_p0_fraction", "kelp")

        p4.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
        p4.set_axes_unit('pc')
        p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p4.annotate_scale(corner='lower_left')
        p4.save("frames_zoom/")

        

