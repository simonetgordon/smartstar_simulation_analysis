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

def apply_annotations_and_save(p, title=None, orient_str=None, w_pc=None, field=None, sim_str=None):
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

    dirname = "frames_" + orient_str + field + "_" + str(sim_str) + "_" + str(w_pc) + "pc/"
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

def _metal_fraction(field, data):
    """
    Compute the metal fraction from given field and data.
    """
    return (data["enzo", "SN_Colour"] / data["gas", "density"]).to("dimensionless")

def process_data_series(es, map, orient, w_pccm, c_min, c_max, root_dir, use_north_vector=False, use_north_vector_ds=False, cmap="viridis", 
                        weight_field="density"):
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
        p = yt.ProjectionPlot(ds, dir, map, center=center, data_source=region, 
                              weight_field=weight_field, width=(w_pccm, 'pc'))
        p.set_cmap(map, cmap)
        p.set_zlim(map, c_min, c_max)
        apply_annotations_and_save(p)

def main(map, orient, w_pccm, c_min, c_max, root_dir):
    """
    Main function to execute the program.
    """
    enzo_file = "smartstar-production-runs.enzo"
    sim = os.path.join(root_dir, enzo_file)
    sim_str = extract_simulation_name(sim)
    es = yt.load_simulation(sim, "Enzo", find_outputs=True)
    es.get_time_series()
    process_data_series(es, map, orient, w_pccm, c_min, c_max, root_dir)


if __name__ == "__main__":
    map = ("gas", "density")
    weight_field = "density"
    orient = "edge-on"
    w_pccm = 30
    w_pc = 1
    c_min = 2e2
    c_max = 3e8
    cmap = "viridis"
    root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSm01"
    main(map, orient, w_pccm, c_min, c_max, root_dir, cmap, weight_field)
