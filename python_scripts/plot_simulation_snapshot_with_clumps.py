import yt
import matplotlib.pyplot as plt
import os
import matplotlib.font_manager as fm
from matplotlib import rc
from helper_functions import extract_simulation_name
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from plot_toomre_q_projection import field_from_sliceplot
from helper_functions import _make_disk_L, ss_properties
from yt.utilities.math_utils import ortho_find

# Define directories
home_dir_1 = "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/"
home_dir_2 = "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/"
clump_dir_1 = "clumps/_1B.RSb01-2clumps"
clump_dir_2 = "clumps/_1B.RSb04clumps"

# Define the snapshots to be plotted
snapshots = ["DD0150_clump_0.h5", "DD0153_clump_0.h5", "DD0156_clump_0.h5", "DD0159_clump_0.h5"]


def setup_plot_env(fontsize, linewidth):
    """Set up the plotting environment"""
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth


# Function to calculate BH age
def age(ds, seed=1):
    """Calculate the age of the BH in Myr"""
    time = ds.current_time.to('Myr')
    birth = 195.59*yt.units.Myr if seed == 2 else 124.76*yt.units.Myr
    return time - birth


def plot_snapshot_with_clumps(sim_name, home_dir, clump_dir, snapshot, ax, width_pc, is_first_col, cmin, cmax, plot_clumps=True):
    """
    Plot a snapshot with clumps and other annotatations.
    Returns the image object for the field from which the colorbar is derived.
    """
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
    p = yt.ProjectionPlot(ds, 'x', ('gas', 'number_density'), width=(width_pc, 'pc'), weight_field=('gas', 'density'), center=center)
    p.set_zlim(('gas', 'number_density'), cmin, cmax)
    p.set_cmap(('gas', 'number_density'), 'viridis')

    # Annotate the clumps
    if plot_clumps:
        p.annotate_clumps(clumps, cmap='prism')

    # Save the plot with clumps to a temporary file
    p.hide_axes()
    p.hide_colorbar()
    p.save("temp.png")

    # Make two images: one for the clumps and one for the field from which the cbar is derived
    img = plt.imread("temp.png")
    n = field_from_sliceplot("number_density", ds, disk, center, width_pc, north, dir, npixels=npixels)
    im = ax.imshow(n, norm=LogNorm())
    im_clump = ax.imshow(img) # Overlay the clumps on the image, musy come after the field image
    im.set_clim(cmin, cmax)
    im_clump.set_clim(cmin, cmax)
    ax.set_axis_off()

    # Add the scalebar and annotations
    font_text = 18
    scalebar = AnchoredSizeBar(ax.transData, npixels/(width_pc/0.5), '0.5 pc', loc='lower left', pad=0.3, color='white', 
                               frameon=False, size_vertical=2, fontproperties=fm.FontProperties(size=font_text))
    ax.add_artist(scalebar)

    ax.text(0.95, 0.05, f"BH Age: {bh_age:.2f}", verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, 
            color='black', fontsize=font_text, bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.5'))
    if is_first_col:
        ax.text(0.05, 0.95, sim_name, verticalalignment='top', horizontalalignment='left', transform=ax.transAxes, 
                color='white', fontsize=font_text, bbox=dict(facecolor='black', alpha=0.5, boxstyle='round,pad=0.5'))

    return im


# Main plotting section
fig, axes = plt.subplots(2, 4, figsize=(20, 9.235))
setup_plot_env(fontsize=16, linewidth=2)
sim_name_1 = extract_simulation_name(home_dir_1)
sim_name_2 = extract_simulation_name(home_dir_2)
width_pc = 2.0
cmin = 5
cmax = 5e7

# Loop over the snapshots and and populate the axes
im = None
for i, snapshot in enumerate(snapshots):
    is_first_col = (i == 0)
    im = plot_snapshot_with_clumps(sim_name_1, home_dir_1, clump_dir_1, snapshot, axes[0, i], width_pc, 
                                   is_first_col, cmin, cmax, plot_clumps=False)

for i, snapshot in enumerate(snapshots):
    is_first_col = (i == 0)
    im = plot_snapshot_with_clumps(sim_name_2, home_dir_2, clump_dir_2, snapshot, axes[1, i], width_pc, 
                                   is_first_col, cmin, cmax, plot_clumps=True)

# Create the colorbar for figure from the last im object (created from sliceplot number density field)
cbar_ax = fig.add_axes([0.90, 0.027, 0.015, 0.943])  # [left, bottom, width, height]
cb = fig.colorbar(im, cax=cbar_ax) 
cb.set_label(r'Number Density ($\rm 1/cm^3$)', fontsize=18)
cb.ax.tick_params(labelsize=18) 

# Save the figure
plt.tight_layout(rect=[0, 0, 0.9, 1])
plt.subplots_adjust(wspace=0, hspace=0)
figname = f"plots/clump_projection_2_rows_{sim_name_1}_{sim_name_2}_{width_pc:.1f}pc.png"
plt.savefig(figname)
print("Saved plot to ", figname)
