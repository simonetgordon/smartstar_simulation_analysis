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
from plot_multi_projections import tidy_data_labels, find_north_vector
import matplotlib.pyplot as plt
from matplotlib import rc
import re # complex str searches
#from yt.extensions.p2p import add_p2p_fields, add_p2p_particle_filters

def apply_annotations_and_save(p, title=None, orient_str=None):
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
    return (data["enzo", "SN_Colour"] / data["gas", "density"] ).to("dimensionless")


##################################  Parameters ###################################

# 1) set map variable (density or temperature or metallicity or h2)
map = "density"

# 2) set image orientation (face-on or edge-on) used to produce north vector
orient = "edge-on"

# 3) set width of box
w_pccm = 30
#w_pccm = 200
#w_pccm = 1000
w_pc = 1 # convert to pc for label
#w_pc = 8 # convert to pc for label
#w_pc = 50 # convert to pc for label

# 4) set colorbar limits
c_min = 2e2 # 2e2 for 1pc no-fb
c_max = 3e8 # 3e8 for 1pc no-fb

# 5) data
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb08" 
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm04"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m01-4dx"
#root_dir = "/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm01"
#root_dir = "/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm04"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m08-4dx/2B.m16-4dx-2"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/40msun/replicating-beckmann/1S.RSm04"
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSm04"
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSmf4"
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSb01"
#root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.m01"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSbf16"
#root_dir = "/ceph/cephfs/sgordon/disk14/cirrus-runs-rsync/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSbf16"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.b04-no-SN"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/40msun/replicating-beckmann/1S.m01-no-SN"
root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m01-no-SN"
#root_dir = "/ceph/cephfs/sgordon/disk14/cirrus-runs-rsync/seed2/270msun"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.b04-no-SN"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSmf16-2/2S.RSmf16-2-gap"
#root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/thermal-fb/1B.th.bf128"
enzo_file = "smartstar-production-runs.enzo"
sim = os.path.join(root_dir, enzo_file)

# 6) use north vector
use_north_vector = False # from one ds
use_north_vector_ds = False # for each ds

#################################################################################

sim_str = tidy_data_labels(extract_simulation_name(sim, 
                                                   #custom_name="s2.DCBH.bf128-free"
                                                   ), 
                           #custom_name="s2.DCBH.bf128-free"
                           )

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

    # find north vector a (near-end) ds
    if use_north_vector:
        dd = "DD0207/DD0207"
        dir, north, disk, orient_str = find_north_vector(dd, root_dir, orient)
    else:
        dir = "z"
        orient_str = "z"
        north = None

    for ds in es.piter():

        ds.add_field(
            name=("gas", "metal_fraction"),
            function=_metal_fraction,
            sampling_type="local",
            units="dimensionless",
        )                                                                                                                                                                               

        ss_pos, ss_mass, ss_age = ss_properties(ds)
        center = center0 if ss_pos is None else ss_pos

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
        else:
            if orient == "face-on":
                dir = "z"
                orient_str = "z"
                north = None
            else:
                dir = "y"
                orient_str = "y"
                north = None

        if nested_level > 0:
            region = ds.box(center0 - 0.5 * my_width,
                            center0 + 0.5 * my_width)
            width = (my_width, "unitary")

        if map == "density":
            # H nuclei number density
            field = "number_density"
            p1 = yt.ProjectionPlot(ds, dir, ("gas", field), width=(w_pccm, 'pccm'), 
                                   #north_vector=north, 
                                   center=center, data_source=region,
                                   weight_field=("gas", field))
            p1.set_cmap(field, 'viridis')
            p1.set_zlim(("gas", field), c_min, c_max)

            apply_annotations_and_save(p1, orient_str=orient_str)
        
        elif map == "temperature":
            # Temperature
            field = "temperature"
            p2 = yt.ProjectionPlot(ds, dir, ("gas", field), width=(w_pccm, 'pccm'), 
                                   #north_vector=north, 
                                   center=center, data_source=region,
                                   weight_field=("gas", "density"))
            p2.set_cmap(field, "RED TEMPERATURE")
            p2.set_zlim(("gas", field), 1e2, 1e4) # 1e2, 2e5 for most
            p2.set_axes_unit('pc')

            apply_annotations_and_save(p2, orient_str=orient_str)

        elif map == "metallicity":

            # add pop3 metallicity fields
            # add_p2p_fields(ds) 

            # Metallicity
            # field = "metal_fraction"
            field = "SN_Colour"
            p3 = yt.ProjectionPlot(ds, dir, ("enzo", field), width=(w_pccm, 'pccm'), 
                                   #north_vector=north, 
                                   center=center, data_source=region,
                                   weight_field=("gas", "density"))
            p3.set_cmap(field, 'kamae')
            #p3.set_zlim(("gas", field), 1e-8, 5e-4) # 1e-13, 5e-4 for 8pc, 
            p3.set_axes_unit('pc')

            apply_annotations_and_save(p3, orient_str=orient_str)

        elif map == "h2":
            # H2 fraction
            field = "H2_p0_fraction"
            p4 = yt.ProjectionPlot(ds, dir, ("gas", field ), north_vector=north, width=(w_pccm,'pccm'), center=center, data_source=region, 
                                   weight_field='cell_mass')
            p4.set_cmap(field, "kelp")
            p4.set_axes_unit('pc')

            apply_annotations_and_save(p4, orient_str=orient_str)
