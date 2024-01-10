import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import re
import os
import yt

import sys
from yt.utilities.math_utils import ortho_find
from helper_functions import ss_properties, tidy_data_labels, extract_simulation_name, _make_disk_L


########################################################################################################
# Produces 2x2 sliceplots of the following fields: 
# number_density, cooling_time, relative_velocity_x, velocity_magnitude, relative_velocity_z
########################################################################################################


def main(directory: str, snapshot: str, plot: str, width_pc: float, dir="z", quiver=False, streamlines=False, disc_r=None):

    ds = yt.load(os.path.join(directory, snapshot)) # load data
    fig = plt.figure()

    # Get the simulation name
    sim_name = tidy_data_labels(extract_simulation_name(directory))
    if not isinstance(sim_name, str):
        print("sim_name = ", sim_name)
        raise TypeError("sim_name is not a string - try add / to end of directory path")

    # See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
    # These choices of keyword arguments produce a four panel plot that includes
    # four narrow colorbars, one for each plot.  Axes labels are only drawn on the
    # bottom left hand plot to avoid repeating information and make the plot less
    # cluttered.
    grid = AxesGrid(
        fig,
        (0.075, 0.075, 0.85, 0.78),
        nrows_ncols=(2, 2),
        axes_pad=0.7,
        label_mode="L",
        share_all=True,
        cbar_location="right",
        cbar_mode="each",
        cbar_size="3%",
        cbar_pad="0%",
    )

    fields = [
        ("gas", "number_density"),
        ("gas", "velocity_magnitude"),
        ("gas", "velocity_cylindrical_radius"),
        ("gas", "velocity_cylindrical_theta"),
        ("gas", "cooling_time"),
        ("gas", "relative_velocity_z"),
        ("gas", "relative_velocity_x"),
    ]

    # Get the particle properties
    ss_pos, ss_mass, ss_age, ss_vel = ss_properties(ds, velocity=True)
    center = "max" if ss_pos is None else ss_pos
    print("ss_vel = ", ss_vel)

    # Set dirction with respect to the disc (if it exists)
    if (len(dir) > 1):
        # make disk data container and define angular momentum vector L
        disc_r_pc = disc_h_pc = disc_r
        disk, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
        vecs = ortho_find(L)
        direction = vecs[0] if dir == "Lface" else vecs[1] if dir == "Ledge" else vecs[2]
        north_vector = vecs[1] if dir == "Lface" else vecs[2] if dir == "Ledge" else vecs[0]
    else:
        direction = dir

    # Set bulk velocity
    dd = ds.all_data()
    dd.set_field_parameter("bulk_velocity", ss_vel)

    # Create the plot. Since SlicePlot accepts a list of fields, we need only
    # do this once.
    if plot == "sliceplot":
        p = yt.SlicePlot(ds, direction, fields, center=ss_pos, width=(width_pc, "pc"), north_vector=north_vector)
    elif plot == "projectionplot":
        p = yt.ProjectionPlot(ds, direction, fields, center=ss_pos, width=(width_pc, "pc"))
    else:
        print("Invalid plot type. Choose 'sliceplot' or 'projectionplot'.")

    # Velocity is going to be both positive and negative, so let's make these
    # slices use a linear colorbar scale
    p.set_log(("gas", "relative_velocity_z"), False)
    p.set_log(("gas", "relative_velocity_x"), False)
    p.set_log(("gas", "velocity_cylindrical_theta"), False)
    p.set_log(("gas", "velocity_cylindrical_radius"), False)
    p.set_log(("gas", "velocity_magnitude"), False)

    # set velocity units to km/s
    p.set_unit("cooling_time", "Myr")
    p.set_unit("relative_velocity_z", "km/s")
    p.set_unit("relative_velocity_x", "km/s")
    p.set_unit("velocity_magnitude", "km/s")
    p.set_unit("velocity_cylindrical_radius", "km/s")
    p.set_unit("velocity_cylindrical_theta", "km/s")

    # set cmaps
    p.set_cmap("number_density", "viridis")
    p.set_cmap("relative_velocity_z", "RdBu")
    p.set_cmap("relative_velocity_x", "RdBu")
    p.set_cmap("velocity_cylindrical_radius", "magma")
    p.set_cmap("velocity_cylindrical_theta", "octarine")
    p.set_cmap("velocity_magnitude", "kelp")
    p.set_cmap("cooling_time", "coolwarm")

    # set clim
    p.set_zlim("relative_velocity_z", -5, 5)
    p.set_zlim("relative_velocity_x", -5, 5)
    p.set_zlim("velocity_magnitude", 2, 40)
    p.set_zlim("number_density", 2e3, 6e8)
    p.set_zlim("cooling_time", 7e-4, 6e3)
    p.set_zlim("velocity_cylindrical_theta", -7, 15)
    p.set_zlim("velocity_cylindrical_radius", -20, 20)

    #p.zoom(2)
    p.set_font({"size": 12})

    if quiver:
        p.annotate_quiver(
            ("gas", "relative_velocity_x"),
            ("gas", "relative_velocity_y"),
            #("gas", "vorticity_y"),
            factor=35,
            color="black",
            headwidth=8,
            #scale=1000,
        )

    # BH position cross
    p.annotate_marker(center, coord_system="data", color="red", marker="x")

    # For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
    # axes.
    for j, field in enumerate(fields[:4]):
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[j].axes
        plot.cax = grid.cbar_axes[j]

        if (j == 0) and (field == ("gas", "number_density")):
            print("j = {}, field = {}".format(j, field))
            # set text coords
            a = 0.03
            b = 0.94
            b2 = 0.03
            p.set_axes_unit('pc')
            p.set_font({"size": 14})

            # top left text
            p.annotate_text((a, b), r"BH Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                                text_args={"color": "white"}) 
            p.annotate_text((a, b-0.06), "BH Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                                text_args={"color": "white"})
            
            # lower right text
            p.annotate_text((0.74, b2), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                                text_args={"color": "white"})
            p.annotate_text((0.65, b2+0.06), "dx = {:.0e}".format(ds.index.get_smallest_dx().in_units('pc')), 
                            coord_system="axis", text_args={"color": "white"})
            p.annotate_text([0.05, 0.05], sim_name, coord_system="axis", text_args={"color": "black"},
                            inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                            "edgecolor": "white", "alpha": 0.5},
            )

            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[j].axes
            plot.cax = grid.cbar_axes[j]
        if streamlines:
            p.annotate_streamlines(
                ("gas", "velocity_x"),
                ("gas", "velocity_y"),
                #("gas", "vorticity_y"),
                linewidth=0.2,
                factor=40,
                #color="lightskyblue",
                cmap="bone",
                #length=0.2,
                #headwidth=8,
                #scale=1000,
            )

    # Finally, redraw the plot on the AxesGrid axes.
    p.render()

    sim_name = sim_name + "_" + snapshot[:6] + "_center"
    field_name = sys.argv[-1]
    fig_name = "plots/{}_multiplot_2x2_{}_{}.png".format(sim_name, field_name, dir)
    plt.savefig(fig_name)
    print("Saved " + fig_name)


if __name__ == "__main__":
    ########################################
    # Call like: python plot_multipanel_velocity.py [field_name]
    directory = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/"
    #directory = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.RSm04/1S.RSm04-145+/"
    #directory = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/"
    #directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm04/"
    directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/"
    directory = "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.RSm04-2/"
    directory = "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.RSm01/"
    directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m08-4dx/"
    directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/"
    directory = "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.RSb16/"

    snapshot = "DD0449/DD0449" # 1S.m04-no-SN 1 Myr
    #snapshot = "DD0232/DD0232" # 1S.RSm04-145+ 1 Myr
    #snapshot = "DD0202/DD0202" # 1B.m16-4dx 1 Myr
    #snapshot = "DD0370/DD0370" # 2B.RSm04 1 Myr
    snapshot = "DD0612/DD0612" # 2B.RSb16 0.67 Myr
    snapshot = "DD0138/DD0138" # 1B.RSm04-2 and 1B.RSm01 1 Myr
    snapshot = "DD0298/DD0298" # 2B.m08-4dx 1 Myr
    snapshot = "DD0279/DD0279" # 2B.RSb08-2 1 Myr
    snapshot = "DD0170/DD0170" # 1B.RSb16 1 Myr
    dir = "Lface" # "x", "y", "z", "Lface", "Ledge" (use angular momentum vector L of disc)
    width_pc = 0.35
    disc_r = 0.05 # set to None for no disc
    quiver = False
    streamlines = True
    plot_type = "sliceplot"
    fields = [
        ("gas", "number_density"),
        ("gas", "cooling_time"),
        ("gas", "relative_velocity_x"),
        ("gas", "velocity_magnitude"),
        ("gas", "relative_velocity_z"),
    ]

    main(directory, snapshot, plot_type, width_pc, dir, quiver, streamlines, disc_r=disc_r)
