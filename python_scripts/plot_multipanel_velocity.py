import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import re
import os
import yt
from smartstar_find import ss_properties
#import colorcet as cc
#import typhon - has nice cmaps

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

if __name__ == "__main__":  

    fn = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/"
    snapshot = "DD0187/DD0187"
    ds = yt.load(os.path.join(fn, snapshot)) # load data

    fig = plt.figure()

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
        ("gas", "relative_velocity_z"),
        ("gas", "relative_velocity_x"),
        ("gas", "velocity_magnitude"),
    ]

    # Create the plot.  Since SlicePlot accepts a list of fields, we need only
    # do this once.
    ss_pos, ss_mass, ss_age, ss_vel = ss_properties(ds, velocity=True)
    center = "max" if ss_pos is None else ss_pos
    print("ss_vel = ", ss_vel)
    dd = ds.all_data()
    #bulk_velocity = dd.quantities.bulk_velocity()
    dd.set_field_parameter("bulk_velocity", ss_vel)
    p = yt.SlicePlot(ds, "y", fields, center=ss_pos, width=(1, "pc"))

    # Velocity is going to be both positive and negative, so let's make these
    # slices use a linear colorbar scale
    p.set_log(("gas", "relative_velocity_z"), False)
    p.set_log(("gas", "relative_velocity_x"), False)

    # set velocity units to km/s
    p.set_unit("relative_velocity_z", "km/s")
    p.set_unit("relative_velocity_x", "km/s")
    p.set_unit("velocity_magnitude", "km/s")

    # set cmaps
    p.set_cmap("number_density", "viridis")
    p.set_cmap("relative_velocity_z", "RdBu")
    p.set_cmap("relative_velocity_x", "RdBu")
    p.set_cmap("velocity_magnitude", "magma")

    # set clim
    p.set_zlim("relative_velocity_z", -5, 5)
    p.set_zlim("relative_velocity_x", -5, 5)
    p.set_zlim("velocity_magnitude", 0.3, 12)
    p.set_zlim("number_density", 2e2, 2e7)

    #p.zoom(2)
    p.set_font({"size": 12})

    p.annotate_quiver(
        ("gas", "relative_velocity_z"),
        ("gas", "relative_velocity_x"),
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
    for j, field in enumerate(fields):
        plot = p.plots[field]
        plot.figure = fig
        plot.axes = grid[j].axes
        plot.cax = grid.cbar_axes[j]

        if (j == 0) and (field == ("gas", "number_density")):
            print("j = {}, field = {}".format(j, field))
            # set text coords
            a = 0.03
            b = 0.95
            b2 = 0.03
            p.set_axes_unit('pc')

            p.set_font({"size": 14})

            # top left text
            p.annotate_text((a, b), r"SS Mass: {:.0f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis",
                                text_args={"color": "white"}) 
            p.annotate_text((a, b-0.06), "SS Age = {:.2f} Myr".format(ss_age[0] / 1e6), coord_system="axis",
                                text_args={"color": "white"})
            
            # lower right text
            p.annotate_text((0.75, b2), "z = {:.2f}".format(ds.current_redshift), coord_system="axis",
                                text_args={"color": "white"})
            p.annotate_text([0.05, 0.05], extract_simulation_name(fn), coord_system="axis", text_args={"color": "black"},
                            inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3,
                                            "edgecolor": "white", "alpha": 0.5},
            )

            plot = p.plots[field]
            plot.figure = fig
            plot.axes = grid[j].axes
            plot.cax = grid.cbar_axes[j]

    # Finally, redraw the plot on the AxesGrid axes.
    p.render()

    name = extract_simulation_name(fn) + "_" + snapshot[:6] + "_center"
    print("saving " + name + "_multiplot_2x2_velocity.png")
    plt.savefig("plots/" + name + "_multiplot_2x2_velocity.png")
    print("Saved plots/"  + name + "_multiplot_2x2_velocity.png")
