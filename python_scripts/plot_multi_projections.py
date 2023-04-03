import yt
import sys
import os
from derived_fields import add_fields_ds
from smartstar_find import ss_properties

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/"
sim = ["1B.RSm01", "1B.RSm04-2", "1B.RSm08-2"]
dds = ["DD0148/DD0148", "DD0158/DD0158", "DD0168/DD0168"]

proj = labels = DS = []
for i, dd in enumerate(dds):

    ds = yt.load(os.path.join(root_dir[i], sim[i], dd))
    add_fields_ds(ds)
    label = "s1_" + str(float(ds.current_time.to('Myr')))[:5] + "_Myr"
    DS.append(ds)
    labels.append(label)

    field = "number_density"
    ss_pos, ss_mass, ss_age = ss_properties(ds)
    center = ss_pos
    widths = [40, 8, 2] * yt.units.pccm
    r = 2000  # pc
    sp = ds.sphere(center, 2 * r)
    for width in widths:
        p = yt.ProjectionPlot(ds, "x", ("gas", field), width=width, center=center, data_source=sp,
                              weight_field='density')

        # annotate
        p.annotate_scale(corner='lower_left')
        p.annotate_timestamp(corner='lower_right', redshift=True, draw_inset_box=True)
        p.annotate_text((0.62, 0.95), "Mass: {:.2f} Msun".format(ss_mass.d), coord_system="axis",
                        text_args={"color": "white"})
        p.annotate_grids(min_level=18, cmap='turbo')
        p.annotate_marker(center, coord_system="data", color="white")  # mark ss position
        p.annotate_sphere(ss_pos, radius=(1.23e-2, "pc"), circle_args={"color": "white"})
        # p.annotate_cell_edges(line_width=0.00002, alpha=0.7, color='white')
        # p.annotate_streamlines(("gas", "relative_velocity_x"), ("gas", "relative_velocity_y"))
        p.annotate_title("BH Age = {:.2f} kyrs".format(ss_age[0] / 1e3))

        # add to p list
        proj.append(p)

mp = multiplot_yt(
    3,
    3,
    proj,
    savefig="yt",
    shrink_cb=0.9,
    bare_axes=True,
    yt_nocbar=False,
    margins=(0.5, 0.5),
)
