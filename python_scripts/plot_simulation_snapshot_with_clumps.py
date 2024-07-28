import yt
import matplotlib.pyplot as plt
import os
from matplotlib import rc
from helper_functions import extract_simulation_name
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.transforms import IdentityTransform
import matplotlib.font_manager as fm
from matplotlib.ticker import LogLocator, LogFormatter, LogFormatterExponent
from plot_toomre_q_projection import field_from_sliceplot
from helper_functions import _make_disk_L, ss_properties   # Import the function from the helper_functions.py file
from yt.utilities.math_utils import ortho_find

# Define directories
home_dir_1B_RSb01_2 = "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/"
home_dir_1B_RSb04 = "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/"
clump_dir_1B_RSb01_2 = "_1B.RSb01-2clumps"
clump_dir_1B_RSb04 = "_1B.RSb04clumps"

# Define the snapshots to be plotted
snapshots = ["DD0150_clump_0.h5", "DD0153_clump_0.h5", "DD0156_clump_0.h5", "DD0159_clump_0.h5"]

def setup_plot_env(fontsize, linewidth):
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth

# Function to calculate BH age
def age(ds, seed=1):
    time = ds.current_time.to('Myr')
    birth = 195.59*yt.units.Myr if seed == 2 else 124.76*yt.units.Myr
    return time - birth

def plot_snapshot_with_clumps(home_dir, clump_dir, snapshot, ax, is_first_col, is_bottom_row, plot_clumps=True):
    ds_path = os.path.join(home_dir, snapshot.replace("_clump_0.h5", ""), snapshot.replace("_clump_0.h5", ""))
    clump_file = os.path.join(clump_dir, snapshot)

    # Load dataset
    ds = yt.load(ds_path)

    # Load clumps
    clumps = []
    try:
        clump_ds = yt.load(clump_file)
        clumps = clump_ds.leaves
    except Exception as e:
        print(f"Error loading clumps: {e}")

    # Extract simulation name
    sim_name = extract_simulation_name(home_dir)

    # Obtain angular momentum vector from small disk and define larger disk for plotting
    bh_age = age(ds, seed=1)
    ss_pos, ss_mass, _ = ss_properties(ds)
    center = ss_pos
    disc_r_pc = 0.3 # pc
    npixels = 800
    disc_h_pc = disc_r_pc
    _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
    vecs = ortho_find(L)
    dir = vecs[0]
    north = vecs[1]
    disc_r_pc_big = disc_h_pc_big = 2 # pc
    disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)
    

    # Create the projection plot
    width_pc = 3.0
    p = yt.ProjectionPlot(ds, 'x', ('gas', 'number_density'), width=(width_pc, 'pc'), weight_field=('gas', 'density'), center=center)
    p.set_zlim(('gas', 'number_density'), 5, 5e7)
    p.set_cmap(('gas', 'number_density'), 'viridis')

    if plot_clumps:
        p.annotate_clumps(clumps, cmap='prism')

    p.hide_axes()
    p.hide_colorbar()
    p.save("temp.png")

    # Read the image and apply the same LogNorm
    img = plt.imread("temp.png")
    n = field_from_sliceplot("number_density", ds, disk, center, width_pc, north, dir, npixels=npixels)
    im = ax.imshow(n, norm=LogNorm())
    im_clump = ax.imshow(img)
    im.set_clim(5, 5e7)
    im_clump.set_clim(5, 5e7)
    ax.set_axis_off()

    # Add the scalebar and annotations
    font_text = 18
    scalebar = AnchoredSizeBar(ax.transData, 800/6, '0.5 pc', loc='lower left', pad=0.3, color='white', frameon=False, size_vertical=2, fontproperties=fm.FontProperties(size=font_text))
    ax.add_artist(scalebar)

    ax.text(0.95, 0.05, f"BH Age: {bh_age:.2f}", verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, color='black', fontsize=font_text, bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.5'))
    if is_first_col:
        ax.text(0.05, 0.95, sim_name, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, color='white', fontsize=font_text, bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.5'))

    return im


# Main plotting section
fig, axes = plt.subplots(2, 4, figsize=(20, 9.25))
setup_plot_env(fontsize=16, linewidth=2)

im = None
for i, snapshot in enumerate(snapshots):
    is_first_col = (i == 0)
    im = plot_snapshot_with_clumps(home_dir_1B_RSb01_2, clump_dir_1B_RSb01_2, snapshot, axes[0, i], is_first_col, False, plot_clumps=False)

for i, snapshot in enumerate(snapshots):
    is_first_col = (i == 0)
    im = plot_snapshot_with_clumps(home_dir_1B_RSb04, clump_dir_1B_RSb04, snapshot, axes[1, i], is_first_col, True, plot_clumps=True)

# Add a single colorbar for the entire figure using the last image object
cbar_ax = fig.add_axes([0.90, 0.027, 0.015, 0.94])  # [left, bottom, width, height]
#im.set_clim(1e1, 1e8)

# Create the colorbar with a LogFormatterExponent for the scientific notation
cb = plt.colorbar(im, cax=cbar_ax) #  extend='both', format=LogFormatterExponent(base=10)

# Set ticks explicitly to ensure they are in logarithmic scale
# locator = LogLocator(base=10, subs='auto', numticks=10)
# cb.set_ticks(locator.tick_values(1, 7))
# cb.update_ticks()  # Ensures the colorbar updates with the new tick settings

# Set label for the colorbar
cb.set_label(r'Number Density ($\rm 1/cm^3$)', fontsize=16)

plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("plots/clump_projection_3pc.png")
print("Saved plot to plots/clump_projection_3pc.png")
plt.show()

