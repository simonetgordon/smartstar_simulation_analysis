import sys
import yt
yt.enable_parallelism()
from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters
import matplotlib as plt # for colormaps?

# Initial star position at formation = [0.49053251, 0.4946698,  0.50963437]
# SN set at 128.1 Myr post BB
# 100 kpc across (not zoomed in close on star forming region)

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

        # Add Pop III metallicity fields
        add_p2p_fields(ds)
        add_p2p_particle_filters(ds)
        
        # Star location, center the movie on it
#        sp = ds.sphere(center, my_width)
#        star_pos = sp['pop3', 'particle_position'].to("unitary")
#        center = star_pos[0].d
        
        if nested_level > 0:
            region = ds.box(center - 0.5 * my_width,
                            center + 0.5 * my_width)
            width = (my_width, "unitary")
        
        
        # For SN timestamp
        time_array = ds.current_time.to('Myr')
        time = time_array.d
        t_sinceSN = time - 128.1

#        # Hii density - defunct
#        p1 = yt.ProjectionPlot(ds, "x", ('gas', 'hii_density'), width=(200, 'pc'), center=center, data_source=region, weight_field="density")
#        p1.set_cmap("hii_density", 'viridis')
#        p1.set_axes_unit("kpc")
#
#        # annotations
#        p1.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
#        p1.annotate_scale(corner='lower_left')
#
#        p1.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#        p1.save("frames_h2-hii/")
#
#        # H_nuclei_density
#        p2 = yt.ProjectionPlot(ds, "x", ('gas', 'H_nuclei_density'), width=(200, 'pc'), center=center, data_source=region, weight_field="density")
#        p2.set_cmap('H_nuclei_density', 'BrBG')
#        p2.set_axes_unit('pc')
#        p2.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
#        p2.annotate_scale(corner='lower_left')
#        p2.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#        p2.save("frames_h2-hii/")
#
#
#        # H_p1_number_density
#        p3 = yt.ProjectionPlot(ds, "x", ("gas", 'H_p1_number_density'), width=(200, 'pc'), center= center, data_source=region, weight_field='density')
#        p3.set_cmap('H_p1_number_density', 'kelp')
#
#        p3.set_axes_unit('pc')
#        p3.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#        p3.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
#        p3.annotate_scale(corner='lower_left')
#        p3.save("frames_h2-hii/")
        
        
        # H_p1_fraction (HII fraction) - weighted by cell_mass
        p4 = yt.ProjectionPlot(ds, "x", ("gas", 'H_p1_fraction'), width=width, center= center, data_source=region, weight_field='cell_mass')
        p4.set_cmap('H_p1_fraction', 'BrBG')

        p4.set_axes_unit('kpccm')
        p4.annotate_text([10,40], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
        p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
        p4.annotate_scale(corner='lower_left')
        p4.save("frames_h2-hii/")


#        # Fe fraction - defunct
#        p4 = yt.ProjectionPlot(ds, "x", ("gas", 'Fe_fraction'), width=(200, 'pc'), center= center, data_source=region, weight_field='density')
#        p4.set_cmap('Fe_fraction', 'ocean')
#
#        p4.set_axes_unit('kpc')
#        p4.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#        p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
#        p4.annotate_scale(corner='lower_left')
#        p4.save("frames_h2-hii/")
#
        
#        # El number density
#        p5 = yt.ProjectionPlot(ds, "x", ("gas", 'El_number_density'), width=(200, 'pc'), center= center, data_source=region, weight_field='density')
#        p5.set_cmap('El_number_density', 'octarine')
#
#        p5.set_axes_unit('pc')
#        p5.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
#        p5.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
#        p5.annotate_scale(corner='lower_left')
#        p5.save("frames_h2-hii/")
