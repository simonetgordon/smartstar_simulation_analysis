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

        center = [0.49053251, 0.4946698,  0.50963437] # centered on initial star position

        
        p1 = yt.ProjectionPlot(ds, "x", ("gas", "hii_density"), width=(200, 'pc'), center=center, data_source=region, weight_field='density')
        p1.set_cmap('density', 'turbo')
        p1.set_axes_unit('pc')
        p1.annotate_timestamp(corner='lower_right')
        p1.annotate_scale(corner='lower_left')
        p1.save("frames_hii/")

