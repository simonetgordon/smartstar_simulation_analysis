"""
Create projection images from simulation to be made into a movie. Call like:
mpirun -np 16 python projection_movie.py
"""

import sys
import yt
import os
yt.enable_parallelism()
from smartstar_find import ss_properties

map = "density"
# set by user
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01"
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

        ss_pos, ss_mass, ss_age = ss_properties(ds, my_width*0.1)

        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        center = ss_pos
        w_pccm = 10

        if map == "density":
            # Gas density
            field = "H_nuclei_density"
            p1 = yt.ProjectionPlot(ds, "z", ("gas", field), width=(w_pccm, 'pccm'), center=center, data_source=region,
                                   weight_field=field)
            p1.set_cmap(field, 'viridis')

            # format
            p1.set_axes_unit('pc')
            p1.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p1.annotate_scale(corner='lower_left')
            p1.annotate_text((0.73, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                             text_args={"color": "white"})
            p1.annotate_title("SS Age = {:.2f} kyrs, {} pccm across".format(ss_age[0] / 1e3, w_pccm))
            dirname = "frames_edge_on_" + field + "_" + str(root_dir[75:]) + "/"
            p1.save(dirname)
        
        elif map == "temperature":
            # Temperature
            field = "temperature"
            p2 = yt.ProjectionPlot(ds, "x", ('gas', field), width=(w_pccm, 'pccm'), center=center,
                                   data_source=region, weight_field='density')
            p2.set_cmap(field, "RED TEMPERATURE")
            p2.set_axes_unit('pc')
            p2.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p2.annotate_scale(corner='lower_left')
            dirname = "frames_" + field + "_" + str(root_dir[75:])
            p2.save(dirname)

        elif map == "metallicity":

            # Metallicity
            field = "metallicity3"
            p3 = yt.ProjectionPlot(ds, "x", ("gas", field), width=(w_pccm, 'pccm'), center= center, data_source=region,
                                   weight_field='density')
            p3.set_cmap(field, 'kamae')
            p3.set_axes_unit('pc')
            p3.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p3.annotate_scale(corner='lower_left')
            dirname = "frames_" + field + "_" + str(root_dir[75:])
            p3.save(dirname)

        elif map == "h2":
            # H2 fraction
            field = "H2_p0_fraction"
            p4 = yt.ProjectionPlot(ds, "x", ("gas", field ), width=(200,'pc'),center=center,
                                   data_source=region, weight_field='cell_mass')
            p4.set_cmap(field, "kelp")

            p4.set_axes_unit('pc')
            p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p4.annotate_scale(corner='lower_left')
            dirname = "frames_" + field + "_" + str(root_dir[75:])
            p3.save(dirname)

        

