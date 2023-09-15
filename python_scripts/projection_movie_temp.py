import yt
import os
import numpy as np
yt.enable_parallelism()
from yt.utilities.math_utils import ortho_find
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from plot_multi_projections import tidy_data_labels
import matplotlib.pyplot as plt
from matplotlib import rc
import re # complex str searches


def apply_annotations_and_save(p, title=None, dirname=None, sim_str=None, center=None, ss_mass=0*yt.units.msun, ss_age=0*yt.units.myr):
    # Set text coords
    a = 0.03
    b = 0.95
    b2 = 0.03
    p.set_axes_unit('pc')

    # Font size
    p.set_font({"size": 24})

    # BH position cross
    p.annotate_marker(center, coord_system="data", color="white")

    # Top left text
    p.annotate_text((a, b), r"BH Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                    text_args={"color": "white"}) 
    p.annotate_text((a, b-0.05), "BH Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                    text_args={"color": "white"})
    
    # Lower right text
    p.annotate_text((0.82, b2), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                    text_args={"color": "white"})
    p.annotate_text([0.05, 0.05], sim_str, coord_system="axis", text_args={"color": "black"},
                    inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                    "edgecolor": "white", "alpha": 0.5})
    if title:
        p.annotate_title(str(title))

    p.save(dirname)


def extract_simulation_name(filepath):
    # Get the last part of the path
    last_part = os.path.basename(filepath)

    # Use regular expression to extract the full simulation name
    match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', last_part)

    if match:
        return match.group(0)

    # If the match is not found, try to extract from the parent directories
    path_parts = filepath.split(os.path.sep)
    for i in range(len(path_parts)-1, -1, -1):
        match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', path_parts[i])
        if match:
            return path_parts[i]

    return None


def make_projection_map(ds, map_type, field, w_pccm, w_pc, c_min, c_max, dirname, sim_str, region, center, dir="z", north=None):

    ss_pos, ss_mass, ss_age = ss_properties(ds)
    if map_type == "density":
        # H nuclei number density
        p = yt.ProjectionPlot(ds, dir, ("gas", field), width=(w_pc, 'pc'), north_vector=north, center=center, data_source=region,
                            weight_field=("gas", field))
        p.set_cmap(field, 'viridis')
        p.set_zlim(("gas", field), c_min, c_max)
        apply_annotations_and_save(p, dirname=dirname, sim_str=sim_str, center=center, ss_mass=ss_mass, ss_age=ss_age)
    
    elif map_type == "temperature":
        # Temperature
        p = yt.ProjectionPlot(ds, dir, ("gas", field), width=(w_pccm, 'pccm'), north_vector=north, center=center, data_source=region,
                            weight_field=("gas", "density"))
        p.set_cmap(field, "RED TEMPERATURE")
        p.set_zlim(("gas", field), 50, 1200)
        p.set_axes_unit('pc')
        apply_annotations_and_save(p, dirname=dirname)

    elif map_type == "metallicity":
        # Metallicity
        p = yt.ProjectionPlot(ds, dir, ("enzo", field), width=(w_pccm, 'pccm'), north_vector=north, center=center, data_source=region, 
                            weight_field=("gas", "density"))
        p.set_cmap(field, 'kamae')
        #p.set_zlim(("enzo", field), c_min, c_max)
        p.set_axes_unit('pc')
        apply_annotations_and_save(p, dirname=dirname)

    elif map_type == "h2":
        # H2 fraction
        p = yt.ProjectionPlot(ds, dir, ("gas", field ), width=(w_pccm,'pccm'), north_vector=north, center=center, data_source=region, 
                            weight_field='cell_mass')
        p.set_cmap(field, "kelp")
        p.set_axes_unit('pc')
        apply_annotations_and_save(p, dirname=dirname)


def main(root_dir, map="density", w_pccm=200, w_pc=1.5, c_min=None, c_max=None, orient="face-on", field="H_nuclei_density", dirname=None):
    sim = os.path.join(root_dir, "smartstar-production-runs.enzo")
    sim_str = tidy_data_labels(extract_simulation_name(root_dir))
    es = yt.load_simulation(sim, "Enzo", find_outputs=True)
    es.get_time_series()

    nested_level = es.parameters.get("CosmologySimulationNumberOfInitialGrids", 0) - 1
    if nested_level == 0:
        region = None
        center0 = [0.5]*3
        width = None
    else:
        left_edge = es.parameters.get("CosmologySimulationGridLeftEdge[%d]" % nested_level, np.zeros(3))
        right_edge = es.parameters.get("CosmologySimulationGridRightEdge[%d]" % nested_level, np.zeros(3))
        my_width = (right_edge - left_edge).max()
        center0 = 0.5 * (left_edge + right_edge)

    # Find north vector at (near-end) ds
    ds_final = yt.load(os.path.join(root_dir, "DD0230/DD0230"))
    ss_pos, ss_mass, ss_age = ss_properties(ds_final)
    center = ss_pos
    r = 2000 * yt.units.pc
    sp = ds_final.sphere(center, 2 * r)

    # Make disk data container and define angular momentum vector L
    disc_r_pc = 2.1
    disc_h_pc = 2.1
    disk, L = _make_disk_L(ds_final, ss_pos, disc_r_pc, disc_h_pc)
    
    # Get orthogonal vectors for face-on or edge-on projection
    vecs = ortho_find(L)
    v = 0

    if orient == "face-on":
        dir = vecs[0]
        north = vecs[1]
    else:
        dir = vecs[1]
        north = vecs[0]
    
    for ds in es.piter():
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        if ss_pos is None:
            center = center0
        else:
            center = ss_pos

        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        dirname = "frames_" + orient + field + "_" + str(sim_str) + "_" + str(w_pc) + "pc/"
        os.makedirs(dirname, exist_ok=True)

        make_projection_map(ds, map, field, w_pccm, w_pc, c_min, c_max, dirname, sim_str, region=region, center=center, dir=dir, north=None)

    print("Saved data to", dirname)

if __name__ == "__main__":
    root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb01"
    field = "H_nuclei_density"
    main(root_dir, map="density", w_pccm=200, w_pc=1.5, c_min=None, c_max=None, orient="face-on", field=field)