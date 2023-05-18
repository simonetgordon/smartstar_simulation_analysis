"""
Create projection images from simulation to be made into a movie. Call like:
mpirun -np 16 python projection_movie.py
"""

import sys
import yt
import os
yt.enable_parallelism()
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find
from plot_multi_projections import tidy_data_labels
import matplotlib.pyplot as plt
from matplotlib import rc
#from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters


def extract_ultimate_directory(filepath):
    # split the filepath into directory and file components
    directory, filename = os.path.split(filepath)
    # split the directory component into its path elements
    path_elements = directory.split(os.path.sep)
    # return the penultimate element, or None if not found
    return path_elements[-1] if len(path_elements) > 1 else None


##################################  Parameters ###################################

# 1) set map variable
map = "density"

# 2) set image orientation (face-on or edge-on) used to produce north vector
orient = "face-on"

# 3) set width of box
w_pccm = 200
#w_pccm = 30
w_pc = 8 # convert to pc for label
#w_pc = 1 # convert to pc for label

# 4) set colorbar limits
c_min = 3e2
c_max = 8e9

# 5) data
root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb08" 
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSb01"
enzo_file = "smartstar-production-runs.enzo"
sim = os.path.join(root_dir, enzo_file)

#################################################################################

# naming sphere container directory
seed = int(root_dir[38])
if seed == 1:
    index = 82
elif seed == 2:
    index = 84
sp_container_dir = root_dir[index:]

sim_str = tidy_data_labels(extract_ultimate_directory(sim))

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

    ## find north vector a (near-end) ds
    ds_final = yt.load(os.path.join(root_dir, "DD0160/DD0160"))
    ss_pos, ss_mass, ss_age = ss_properties(ds_final)
    center = ss_pos
    r = 2000*yt.units.pc
    sp = ds_final.sphere(center, 2 * r)

    # make disk data container and define angular momentum vector L
    disc_r_pc = 2.1
    disc_h_pc = 2.1
    disk, L = _make_disk_L(ds_final, ss_pos, disc_r_pc, disc_h_pc)
    

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    vecs = ortho_find(L)
    v = 0

    if orient == "face-on":
        orient_str = "face_on_"
        dir = vecs[0]
        north = vecs[1]
    else:
        orient_str = "edge_on_"
        dir = vecs[1]
        north = vecs[0]
    
    for ds in es.piter():

        # add pop3 metallicity fields
        #add_p2p_fields(ds)                                                                                                                                                                                

        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = ss_pos

        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        if map == "density":
            # h number density
            field = "number_density"
            p1 = yt.ProjectionPlot(ds, dir, ("gas", field), width=(w_pccm, 'pccm'), north_vector=north, center=center, data_source=region,
                                   weight_field=("gas", field))
            p1.set_cmap(field, 'viridis')
            p1.set_zlim(("gas", field), c_min, c_max)

            # format
            a = 0.65
            b = 0.95
            b2 = 0.03
            p1.set_axes_unit('pc')

            # annotations
            # if swap_axes:
            #     p1.swap_axes()
            #     a, b = b, a
            #     p1.annotate_timestamp(x_pos=0.10, y_pos=0.87, redshift=True, draw_inset_box=True, coord_system="axis")
            # else:
            #     p1.annotate_timestamp("lower_right", redshift=True, draw_inset_box=True, coord_system="axis")
            #p1.annotate_scale(corner='lower_left')
            p1.set_font({"size": 24})
            p1.annotate_marker(center, coord_system="data", color="white")  # mark ss position
            p1.annotate_text((a, b), r"BH Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                             text_args={"color": "white"})
            p1.annotate_text((a, b-0.05), "BH Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                             text_args={"color": "white"})
            p1.annotate_text((0.82, b2), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                             text_args={"color": "white"})
            #p1.annotate_title("{} pccm across".format(w_pccm))
            p1.annotate_text([0.05, 0.05], sim_str, coord_system="axis", text_args={"color": "black"},
                            inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                            "edgecolor": "white", "alpha": 0.5},
                            )

            dirname = "frames_" + orient_str + field + "_" + str(sim_str) + "_" + str(w_pc) + "pc/"
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
            p2.annotate_text((0.70, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                             text_args={"color": "white"})
            p2.annotate_title("SS Age = {:.2f} kyrs, {} pccm across".format(ss_age[0] / 1e3, w_pccm))
            dirname = "frames_" + field + "_" + str(root_dir[75:]) + "/"
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
            dirname = "frames_" + field + "_" + str(root_dir[75:]) + "/"
            p3.save(dirname)

        elif map == "h2":
            # H2 fraction
            field = "H2_p0_fraction"
            p4 = yt.ProjectionPlot(ds, "x", ("gas", field ), width=(w_pccm,'pccm'), center=center,
                                   data_source=region, weight_field='cell_mass')
            p4.set_cmap(field, "kelp")
            p4.set_axes_unit('pc')
            p4.annotate_timestamp(corner="lower_right", redshift=True, draw_inset_box=True)
            p4.annotate_scale(corner='lower_left')
            dirname = "frames_" + field + "_" + str(root_dir[75:]) + "/"
            p3.save(dirname)

        

