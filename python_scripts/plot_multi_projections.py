import yt
import sys
import os
from smartstar_find import ss_properties

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16-2"


for i in range(len(sys.argv)):
    input_dd = sys.argv[i]
    ds = yt.load(os.path.join(root_dir, sys.argv[i]))






mp = multiplot_yt(
    2,
    2,
    [p1, p2, p3, p4],
    savefig="yt",
    shrink_cb=0.9,
    bare_axes=True,
    yt_nocbar=False,
    margins=(0.5, 0.5),
)