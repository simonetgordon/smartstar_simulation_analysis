import sys
import yt
import os
yt.enable_parallelism()
from smartstar_find import ss_pos, ss_age, ss_mass

field = "density"
# set by user
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16"
enzo_file = "smartstar-production-runs.enzo"
sim = os.path.join(root_dir, enzo_file)
if __name__ == "__main__":
    es = yt.load_simulation(sim, "Enzo", find_outputs=True)
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
        
        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        center = ss_pos
        w_pccm = 10

        if field == "density":
            # Gas density
            p1 = yt.ProjectionPlot(ds, "x", ("gas", "density"), width=(w_pccm, 'pccm'), center=center, data_source=region,
                                   weight_field='density')
            p1.set_cmap('density', 'turbo')
            p1.set_axes_unit('pc')
            p1.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p1.annotate_scale(corner='lower_left')
            p1.annotate_text((0.73, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                             text_args={"color": "white"})
            p1.annotate_title("SS Age = {:.2f} kyrs, {} pccm across".format(ss_age[0] / 1e3, w_pccm))

            p1.save("frames_gas/")
        
        else:
            # Temperature
            p2 = yt.ProjectionPlot(ds, "x", ('gas', 'temperature'), width=(200,'pc'), center=center, data_source=region, weight_field='density')
            p2.set_cmap("temperature", "RED TEMPERATURE")
            p2.set_axes_unit('pc')
            p2.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p2.annotate_scale(corner='lower_left')
            #p2.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
            p2.save("frames_gas/")


            # Metallicity
            p3 = yt.ProjectionPlot(ds, "x", ("gas", "metallicity3"), width=(200, 'pc'),center= center, data_source=region, weight_field='density')
            p3.set_cmap('metallicity3', 'kamae')

            p3.set_axes_unit('pc')
            #p3.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
            p3.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p3.annotate_scale(corner='lower_left')
            p3.save("frames_metal/")


            # H2 fraction
            p4 = yt.ProjectionPlot(ds, "x", ("gas", "H2_p0_fraction"), width=(200,'pc'),center=center, data_source=region, weight_field='cell_mass')
            p4.set_cmap("H2_p0_fraction", "kelp")

            #p4.annotate_text([20,90], f"Time Since SN = {t_sinceSN:.1f} Myr", coord_system="plot")
            p4.set_axes_unit('pc')
            p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p4.annotate_scale(corner='lower_left')
            p4.save("frames_zoom/")

        

