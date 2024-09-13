import yt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from helper_functions import ss_properties, _make_disk_L, ortho_find, field_from_sliceplot, extract_simulation_name, setup_plot_env
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from matplotlib import rc

def setup_plot_env(fontsize, linewidth):
    """Set up the plotting environment"""
    rc('font', **{'family': 'serif', 'serif': ['Times'], 'weight': 'light'})
    rc('text', usetex=True)
    plt.rcParams["mathtext.default"] = "regular"
    plt.rcParams['font.size'] = fontsize
    plt.rcParams['lines.linewidth'] = linewidth


def check_data_validity(data):
    if not np.isfinite(data).all():
        raise ValueError("Data contains non-finite values (NaN or Inf), please check the input data.")

def field_from_projectionplot(field, ds, disk, width_pc, dir, npixels=2048, radius=False):
    """
    Compute field from a projection plot of a dataset over a disk.
    Surface Density = projection plot density * cell height
    """
    p = yt.ProjectionPlot(ds, dir, ("gas", field), center=disk.center, width=(width_pc, "pc"), data_source=disk)
    proj_frb = p.data_source.to_frb((width_pc, "pc"), npixels)
    proj_field = proj_frb[("gas", field)]
    
    if radius:
        radius = proj_frb['radius'].to('pc')
        return proj_field, radius
    else:
        return proj_field


def make_projection_snapshot(ds, disk, ax, width_pc, dir, field, cmin, cmax, center, L, npixels=800, plot_clumps=False):
    """
    Plot a snapshot with projections of either temperature or density over a disk.
    Returns the image object for the field from which the colorbar is derived.
    """
    # Make the disk projection plot in yt and save to temporary file
    vecs = ortho_find(L)  # Get orthogonal directions based on disk angular momentum
    dir = vecs[0] if dir == "face-on" else vecs[1]  # Choose the correct direction based on input
    north = vecs[1]
    
    # Create the projection plot
    p = yt.ProjectionPlot(ds, dir, ('gas', field), width=(width_pc, 'pc'), weight_field=('gas', 'density'), center=center)
    p.set_zlim(('gas', field), cmin, cmax)
    p.set_cmap(('gas', field), 'magma' if field == 'temperature' else 'viridis')

    # Annotate clumps if needed
    if plot_clumps:
        p.annotate_clumps(clumps, cmap='prism')

    # Save the plot with clumps to a temporary file
    p.hide_axes()
    p.hide_colorbar()
    temp_filename = "plots/temp.png"
    p.save(temp_filename)

    # Load the saved projection image into `img`
    img = plt.imread(temp_filename)

    # Use field_from_sliceplot to create the main field image
    n = field_from_sliceplot(field, ds, disk, center, width_pc, north, dir, npixels=npixels)

    # Plot the data field using imshow
    im_slc = ax.imshow(n, origin='lower', norm=LogNorm(vmin=cmin, vmax=cmax), cmap='magma' if field == 'temperature' else 'viridis')
    
    # Overlay the clump image using imshow
    im_proj = ax.imshow(img)  # Overlay the clumps image

    # Set color limits for both images
    im_slc.set_clim(cmin, cmax)
    im_proj.set_clim(cmin, cmax)
    
    # Hide axis
    ax.set_axis_off()

    return im_slc


def make_split_projection_plot(datasets, box_width_pc=1, npixels=800, nmin=1e4, nmax=1e9 ,tempmin=1e1, tempmax=1e8, dir='face-on', sim_name=None, scale_size=0.5):
    """
    Create a multipanel figure with temperature and number density projections.
    Each row represents a different dataset, and each column represents a different time snapshot.
    """
    fig, axes = plt.subplots(2, len(datasets), figsize=(8.5, 6))  # 2 rows for temp & number density, columns for time snapshots
    setup_plot_env(fontsize=14, linewidth=1.5)

    # Loop through datasets and times to generate the projections
    for i, ds in enumerate(datasets):
        ss_pos, ss_mass, ss_age = ss_properties(ds)
        r = 0.1  # Radius of the disk
        disk, L = _make_disk_L(ds, ss_pos, width=r, height=r)  # Get the disk angular momentum vector

        # Temperature projection (top row)
        im_temp = make_projection_snapshot(ds, disk, axes[0, i], box_width_pc, dir, 'temperature', 
                                           cmin=tempmin, cmax=tempmax, center=ss_pos, L=L, npixels=npixels)

        # Number density projection (bottom row)
        im_density = make_projection_snapshot(ds, disk, axes[1, i], box_width_pc, dir,'number_density',
                                              cmin=nmin, cmax=nmax, center=ss_pos, L=L, npixels=npixels)

        # Add the scalebar and annotations
        font_text = 18
        scalebar = AnchoredSizeBar(axes[0, i].transData, npixels/(box_width_pc/scale_size), f'{scale_size} pc', loc='lower left', pad=0.2, color='white', 
                                frameon=False, size_vertical=2, fontproperties=fm.FontProperties(size=font_text))
        axes[1, i].add_artist(scalebar)

        axes[0, i].text(0.05, 0.95, f"BH Age: {ss_age[0]/1e6:.2f} Myr", verticalalignment='top', horizontalalignment='left', transform=axes[0,i].transAxes, 
                color='black', fontsize=font_text-4, bbox=dict(facecolor='white', alpha=0.5, boxstyle='round,pad=0.5'))

    # Create colorbars
    #fig.subplots_adjust(right=0.85)
    fig.subplots_adjust(left=0.05, right=0.85, top=0.85, bottom=0.1, wspace=0.01, hspace=0.01)
    cbar_ax_temp = fig.add_axes([0.87, 0.48, 0.02, 0.37])  # For temperature (top row) [left, bottom, width, height]
    cbar_ax_density = fig.add_axes([0.87, 0.1, 0.02, 0.37])  # For number density (bottom row)

    # Add colorbars to the side of the plot
    im_temp.set_clim(tempmin, tempmax)
    im_density.set_clim(nmin, nmax)
    fig.colorbar(im_temp, cax=cbar_ax_temp).set_label('Temperature [K]')
    fig.colorbar(im_density, cax=cbar_ax_density).set_label(r'Number Density [cm$^{-3}$]')
    
    # Add simulation name to the top of the plot
    if sim_name:
        fig.suptitle(sim_name, fontsize=16, y=0.90)

    #plt.tight_layout(rect=[0, 0, 0.85, 1])  # Ensure layout fits
    figname = f"{sim_name}_{dir}_{box_width_pc}pc_split_projection.png" if sim_name else "split_projection.png"
    dirname = "plots/projections_temp_density/"
    #plt.tight_layout()
    plt.savefig(dirname + figname, dpi=300)
    print(f"Saved split projection plot to {dirname}{figname}")


# Example usage with datasets and times
filepaths = [
    # 1B.resim.th.b01-3-eta-0.1
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0450/DD0450", # 31.70
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0599/DD0599", # 31.80
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695", # 31.90

    # 1B.resim.th.b01-3-eta-0.01-fb-r=10dx
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0445/DD0445",
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/fb-radius=10dx/DD0603/DD0603",
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/fb-radius=10dx/DD0702/DD0702"

    # 1B.resim.th.b01-3-eta-0.01
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0445/DD0445", # 31.70
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0543/DD0543", # 31.80
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0646/DD0646", # 31.90

    # 1B.resim.th.b04 - end event 2.80
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0153/DD0153", # 2.500
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0353/DD0353", # 2.700
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0552/DD0552", # 2.900

    # 1B.resim.th.b04-eta-0.1
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0153/DD0153", # 2.501
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0352/DD0352", # 2.700
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0421/DD0421", # 2.900

    # 1B.resim.th.b01-eta-0.1-2
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0155/DD0155", # 2.701
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0304/DD0304", # 2.850
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0454/DD0454", # 3.000


    # 1B.resim.th.b04-eta-0.1-fb-r-10dx
    "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0153/DD0153", # 2.500
    "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0352/DD0352", # 2.700
    "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0450/DD0450", # 2.800

    # 1B.resim.th.b04-eta-0.1-fb-r-7dx
    #"/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0153/DD0153", # 2.500
    #"/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0302/DD0302", # 2.650
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0352/DD0352", # 2.700
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0377/DD0377", # 2.725
    # "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0402/DD0402", # 2.750
    #"/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0419/DD0419"
]

datasets = [yt.load(filepath) for filepath in filepaths]

sim_name = extract_simulation_name(datasets[0], custom_name="1B.resim.th.b04-eta-0.1-fb-r-10dx")
print(f"Creating split projection plot for {sim_name}...")
#make_split_projection_plot(datasets, box_width_pc=1, nmin=1e2, nmax=8e7 ,tempmin=2e1, tempmax=8e6, dir='face-on', sim_name=sim_name, scale_size=0.2)
make_split_projection_plot(datasets, box_width_pc=10, nmin=1e0, nmax=6e7 ,tempmin=2e2, tempmax=6e5, dir='edge-on', sim_name=sim_name, scale_size=1)
