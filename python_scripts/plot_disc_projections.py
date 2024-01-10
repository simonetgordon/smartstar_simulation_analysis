import yt
import sys
import os
import numpy as np
from helper_functions import ss_properties
import matplotlib.pyplot as plt
from matplotlib import pyplot
from matplotlib import rc
from derived_fields import add_fields_ds
from yt.utilities.math_utils import ortho_find
from helper_functions import tidy_data_labels
import re
from helper_functions import extract_dd_segment, extract_simulation_name, _make_disk_L, configure_font


if __name__ == "__main__":

    root_dir = "/Backup00/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/"
    input = sys.argv[1]
    ds = yt.load(os.path.join(root_dir, sys.argv[1]))
    add_fields_ds(ds)
    sim_name = tidy_data_labels(extract_simulation_name(ds.directory))
    dd_name = extract_dd_segment(ds.directory)

    # configure font
    configure_font(fontsize=48)

    # grab bh particle properties
    ss_pos, ss_mass, ss_age = ss_properties(ds)

    # make disk data container and define angular momentum vector L
    disc_r_pc = 0.2
    disc_h_pc = 0.01
    _, L = _make_disk_L(ds, ss_pos, disc_r_pc, disc_h_pc)
    disc_r_pc = 1.0
    disk = ds.disk(ss_pos, L, disc_r_pc*yt.units.pc, disc_h_pc*yt.units.pc)

    # Gives a 3d vector and it will return 3 orthogonal vectors, the first one being the original vector
    # and the 2nd and 3rd being two vectors in the plane orthogonal to the first and orthogonal to each other.
    # It's very useful for giving you a vector edge-on to the disk.
    dir, north, _ = ortho_find(L)

    #north = vecs[0] if i > 0 else vecs[1]
    dir = "x"
    p = yt.ProjectionPlot(ds, dir, ("gas", "number_density"), weight_field=("gas", "density"), center=disk.center, width=(0.25, "pc"), data_source=disk)
    p.annotate_marker(ss_pos, coord_system='data', plot_args={'color':'white', 's':100})
    p.annotate_grids(alpha=0.5, min_level=18, max_level=18)
    p.set_axes_unit("pc")
    #p.set_clim(1e3, 1e9)
    p.set_cmap(("gas", "number_density"), "viridis")

    # annotate information
    a = 0.03
    b = 0.94
    b2 = 0.03
    p.set_axes_unit('pc')
    p.set_font({"size": 24})

    # top left text
    p.annotate_text((a, b), r"BH Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                        text_args={"color": "white"}) 
    p.annotate_text((a, b-0.06), "BH Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                        text_args={"color": "white"})
    
    # lower right text
    p.annotate_text((0.76, b2+0.06), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                        text_args={"color": "white"})
    p.annotate_text((0.76, b2), "dx = {:.0e}".format(ds.index.get_smallest_dx().in_units('pc')), 
                    coord_system="axis", text_args={"color": "white"})
    p.annotate_text([0.05, 0.05], sim_name, coord_system="axis", text_args={"color": "yellow", "fontsize": 48},
                    inset_box_args={"boxstyle": "square,pad=0.5", "facecolor": "black", "linewidth": 5,
                                    "edgecolor": "white", "alpha": 0.5},
    )

    # annotate sphere around BH
    dx = ds.index.get_smallest_dx().in_units('pc')
    r_sphere = dx*4
    #p.annotate_sphere(ss_pos, r_sphere, circle_args={"color": "white"}, coord_system='data', text=None, text_args=None)

    # ticks and tick labels
    print("Ticks: ", p.plots["number_density"].axes.get_xticks())
    # ax = p.plots["number_density"].axes
    # ax.set_xticklabels([])
    # ax.set_yticklabels([])
    # ax.set_xticks([-0.08, -0.04, 0.0, 0.04, 0.08])
    # ax.set_yticks([-0.08, -0.04, 0.0, 0.04, 0.08])
    # ax.set_xticklabels(['-0.08', '-0.04', '0', '0.04', '0.08'], fontsize=24)
    # ax.set_yticklabels(['-0.08', '-0.04', '0', '0.04', '0.08'], fontsize=24)
    # p._setup_plots()
    
    # save figure
    fig_name = "disc_gridlines_{}_{}.pdf".format(sim_name, dd_name)
    p.set_figure_size((8, 8))
    p.save(fig_name, mpl_kwargs={"bbox_inches": "tight"})
    print("Saved " + fig_name)

    # for field in ("height", "cylindrical_radius"):
    #     p = yt.ProfilePlot(disk, ("index", field), ("gas", "density"), weight_field=("gas", "cell_mass"))
    #     p.set_unit(("index", field), "pc")
    #     p.save()