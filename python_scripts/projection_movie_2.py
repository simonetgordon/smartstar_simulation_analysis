import os
import yt
from yt.utilities.math_utils import ortho_find
import sys
yt.enable_parallelism()
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from yt.utilities.math_utils import ortho_find
from plot_multi_projections import tidy_data_labels, find_north_vector
import matplotlib.pyplot as plt
from matplotlib import rc
import re # complex str searches

def apply_annotations_and_save(p, ds, center, ss_mass, ss_age, title=None, orient_str=None, w_pc=None, map=None, sim_str=None, plot_type=None):
    """
    Apply annotations to the projection and save the result.
    """
    # set text coords
    a = 0.03
    b = 0.95
    b2 = 0.03
    p.set_axes_unit('pc')

    # font size
    p.set_font({"size": 24})

    # BH position cross
    p.annotate_marker(center, coord_system="data", color="white")

    # top left text
    p.annotate_text((a, b), r"SS Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                        text_args={"color": "white"}) 
    p.annotate_text((a, b-0.05), "SS Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                        text_args={"color": "white"})
    
    # lower right text
    p.annotate_text((0.82, b2), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                        text_args={"color": "white"})
    p.annotate_text([0.05, 0.05], sim_str, coord_system="axis", text_args={"color": "black"},
                    inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                    "edgecolor": "white", "alpha": 0.5},
                    )
    if title:
        p.annotate_title(str(title))
    
    # save plot
    field = extract_field(map)
    dirname = "frames_" + orient_str + field +  "_" + str(plot_type) + "_" + str(sim_str) + "_" + str(w_pc) + "pc/"
    p.save(dirname)

    return ("saved to frames to {}".format(dirname))

def extract_simulation_name(filepath, custom_name=None):
    # Get the last part of the path
    last_part = os.path.basename(filepath)

    # Use regular expression to extract the full simulation name
    match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', last_part)

    if custom_name:
        return custom_name

    elif match:
        return match.group(0)

    # If the match is not found, try to extract from the parent directories
    path_parts = filepath.split(os.path.sep)
    for i in range(len(path_parts)-1, -1, -1):
        match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', path_parts[i])
        if match:
            return path_parts[i]

    return None


def extract_field(map_tuple):
    """
    Extracts the second element from a tuple and returns it.
    """
    if not isinstance(map_tuple, tuple) or len(map_tuple) < 2:
        raise ValueError("Input must be a tuple with at least two elements")
    
    return map_tuple[1]


def _metal_fraction(field, data):
    """
    Compute the metal fraction from given field and data.
    """
    return (data["enzo", "SN_Colour"] / data["gas", "density"]).to("dimensionless")

def process_data_series(es, map, orient, w_pccm, c_min, c_max, root_dir, use_north_vector=False, use_north_vector_ds=False, cmap="viridis", 
                        weight_field="density", plot_type="projection"):
    """
    Process the data series.
    """
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

    # find north vector a (near-end) ds or use the z-axis 
    if use_north_vector:
        dd = "DD0207/DD0207"
        dir, north, disk, orient_str = find_north_vector(dd, root_dir, orient)
    else:
        if orient == "face-on":
            dir = "z"
            orient_str = "z"
            north = None
        else:
            dir = "y"
            orient_str = "y"
            north = None

    for ds in es.piter():

        # add metal fraction field if map is metal_fraction
        if map == "metal_fraction":
            ds.add_field(("gas", "metal_fraction"), function=_metal_fraction, 
                         units="dimensionless", display_name="Metal Fraction")
            
        # find the location of the smart star and set it as the center (if it exists)
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = center0 if ss_pos is None else ss_pos

        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width, center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        # find north vector for each ds
        if use_north_vector_ds:
            # make disk data container and define angular momentum vector L
            disc_r_pc = disc_h_pc = 0.1
            disk, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
            vecs = ortho_find(L)

            if orient == "face-on":
                orient_str = "face_on_continuous_"
                dir = vecs[0]
                north = vecs[1]
            else:
                orient_str = "edge_on_continuous_"
                dir = vecs[2]
                north = vecs[0]

        # make projection plot
        if plot_type == "projection":
            p = yt.ProjectionPlot(ds, dir, map, center=center, data_source=region, 
                                  weight_field=weight_field, width=(w_pccm, 'pccm'))
        elif plot_type == "slice":
            p = yt.SlicePlot(ds, dir, map, center=center, data_source=region, 
                             width=(w_pccm, 'pccm'), north_vector=north)
        p.set_cmap(map, cmap)
        p.set_zlim(map, c_min, c_max)
        apply_annotations_and_save(p, ds, center, ss_mass, ss_age, title=None, orient_str=orient_str, w_pc=w_pccm, 
                                   map=map, sim_str=extract_simulation_name(root_dir), plot_type=plot_type)

def main(map, orient, w_pccm, c_min, c_max, root_dirs, cmap, weight_field, plot_type="projection", use_north_vector_ds=False):
    """
    Main function to execute the program.
    """
    enzo_file = "smartstar-production-runs.enzo"
    for root_dir in root_dirs:
        sim = os.path.join(root_dir, enzo_file)
        sim_str = extract_simulation_name(sim)
        print("Processing simulation: {}".format(sim_str))
        es = yt.load_simulation(sim, "Enzo", find_outputs=True)
        es.get_time_series()
        process_data_series(es, map, orient, w_pccm, c_min, c_max, root_dir,  use_north_vector_ds=use_north_vector_ds, cmap=cmap, 
                            weight_field=weight_field, plot_type=plot_type)


if __name__ == "__main__":
    ##########################################################
    # Generates a series of projections for a given simulation
    # Call like: mpirun -np 16 python projection_movie_2.py 
    ##########################################################
    map = ("gas", "number_density")
    weight_field = "density"
    orient = "face-on" # "edge-on" or "face-on"
    w_pccm = 30
    w_pc = 1
    c_min = 2e2
    c_max = 3e8
    cmap = "viridis" # e.g. "magma", "viridis", "plasma", "inferno"
    plot_type = "slice" # "projection" or "slice"
    root_dirs = [
        #"/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/",
        "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/",
        "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/",
        "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/"
        ]
    main(map, orient, w_pccm, c_min, c_max, root_dirs, cmap, weight_field, plot_type, use_north_vector_ds=True)
