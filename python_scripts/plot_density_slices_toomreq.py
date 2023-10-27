import yt
import cmyt
import os
from matplotlib import ticker
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from plot_multi_projections import tidy_data_labels
from yt.utilities.math_utils import ortho_find
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from plot_toomre_q_projection import toomre_from_sliceplot, field_from_sliceplot
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.offsetbox import AnchoredText
from plot_radial_profile_from_frb import extract_simulation_name, extract_dd_segment
import matplotlib.colors as mcolors
from find_fourier_modes import get_theta_values, find_bar_radius

def find_fourier_modes_and_phase_angles(radii, radius_pc, densities, theta, dV, dr=0.001):
    """
    Compute the m=1, 2 Fourier mode strengths and phase angles for a given annular region of a disk.
    Input:
        radii: list of radii for annular regions
        radius_pc: 2D array of radii values for each cell in the disk
        densities: 2D array of densities values for each cell in the disk
        theta: 2D array of theta values for each cell in the disk
        dV: volume of each cell in the disk
        dr: thickness of annular regions
    Output:
        m1_strengths: list of m=1 Fourier mode strengths for each annular region
        m2_strengths: list of m=2 Fourier mode strengths for each annular region
        phi_1_values: list of m=1 Fourier mode phase angles for each annular region
        phi_2_values: list of m=2 Fourier mode phase angles for each annular region
    """
    m1_strengths = []
    m2_strengths = []
    phi_1_values = []
    phi_2_values = []
    for r in radii:
        mask = (radius_pc >= r) & (radius_pc < r + dr)

        # Get densities and thetas for this annular region, where 0 < theta < 2pi
        masked_densities = densities[mask]
        theta_2d = np.tile(theta, (densities.shape[0], 1))
        masked_theta = theta_2d[mask]

        # Compute the mass-equivalent for each cell in this region
        mass_equivalent = masked_densities * dV

        # Compute a_2 and b_2 coefficients for m=2 mode
        a_2 = np.sum(mass_equivalent * np.cos(2 * masked_theta))
        b_2 = np.sum(mass_equivalent * np.sin(2 * masked_theta))

        # for the m = 1 mode
        a_1 = np.sum(mass_equivalent * np.cos(masked_theta))
        b_1 = np.sum(mass_equivalent * np.sin(masked_theta))

        # Compute A_0 for this region
        A_0 = np.sum(mass_equivalent)

        # Compute A_2 for this region
        A_2 = np.sqrt(a_2**2 + b_2**2)
        A_1 = np.sqrt(a_1**2 + b_1**2)

        # Compute the bar strength for this region
        bar_strength = A_2 / A_0 if A_0 != 0 else 0
        m1_strength = A_1 / A_0 if A_0 != 0 else 0
        m2_strengths.append(bar_strength)
        m1_strengths.append(m1_strength)

        # Compute the phase angle phi_2 for this region, where -90 < phi_2 < 90
        phi_2 = 0.5 * np.degrees(np.arctan2(b_2, a_2).value)
        phi_1 = np.degrees(np.arctan2(b_1, a_1).value)
        phi_2_values.append(phi_2)
        phi_1_values.append(phi_1)

    return m1_strengths, m2_strengths, phi_1_values, phi_2_values



root_dir = [#"/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
            #"/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/"
            "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/"
            ]
sim = [#"2B.RSm04", 
       #"2B.m08-4dx"
       "1B.m16-4dx"
       #"2B.RSb08-2"
       ] 

dds = ["DD0167/DD0167", "DD0178/DD0178", "DD0189/DD0189", "DD0231/DD0231"]  # 0.39, 0.50, 0.6, 1 Myr for 1B.m16,
#dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
#dds = ["DD0219/DD0219", "DD0227/DD0227", "DD0236/DD0236", "DD0279/DD0279"]  # 0.2, 0.69, 1 Myr for 2B.b08,

DS = []
for s in range(len(dds)):
    ds = yt.load(os.path.join(root_dir[0], sim[0], dds[s]))
    DS.append(ds)

#### PLOT ####

fig = plt.figure(figsize=(10, 12))

# IMAGE GRID #
nrows = 4
ncols = 3
#grid = ImageGrid(fig, (0.1, 0.12, 0.60, 0.95), nrows_ncols=(nrows, ncols), axes_pad=0.0, share_all=True, cbar_mode=None)
#grid = ImageGrid(fig, (0.25, 0.0, 0.75, 1.0), nrows_ncols=(nrows, ncols), axes_pad=0.0, share_all=True, cbar_mode=None)
plt.rcParams['text.usetex'] = True

# Plot data and annotations
grid_index = 0
size = 0.16 # subplot dims
fourier_width = size*1.5
for row in range(nrows):
    ax_fourier = fig.add_axes([0.02, 0.8 - row*size, size, size], frame_on=True)

    # Load dataset and define axis
    ds = DS[row]
    sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
    sim_label = sim_label.replace("-2", "")
    sim_label = sim_label.replace("RS", "")
    print("Plotting " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))

    # Grab bh properties and define center, width and resolution of sliceplots
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    center = ss_pos
    width_pc = 0.32
    tick_labels = ['', '-0.08', '0.0', '0.08', '']
    npixels = 2048
    dx = ds.index.get_smallest_dx().in_units('cm')

    # Obtain angular momentum vector from small disk and define larger disk for plotting
    disc_r_pc = disc_h_pc = 0.01
    _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
    vecs = ortho_find(L)
    dir = vecs[0]
    north = vecs[1]
    disc_r_pc_big = disc_h_pc_big = 1.0
    disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)

    # fourier modes
    density, radius_pc = field_from_sliceplot("density", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=True)
    surface_density = density * dx # g/cm^2
    
    # List of radii to define annular regions with thickness dr
    dr = 0.001
    r_min = np.min(radius_pc).value
    r_max = np.max(radius_pc).value
    radii = np.arange(r_min, r_max + dr, dr) # 73

    # Compute bar strength and phase angle variability across discrete annular regions
    dV = dx**3
    theta = get_theta_values(surface_density)
    m1_strengths, m2_strengths, phi_1_values, phi_2_values = find_fourier_modes_and_phase_angles(radii, radius_pc, density, theta, dV, dr)

    # Find radius of the bar feature
    var_deg = 3 # degrees
    bar_radius, i = find_bar_radius(phi_2_values, radii, var_deg=var_deg)
    bar_strength = m2_strengths[i]

    # Plot bar strength and phase angle variability across discrete annular regions
    c1 = (0.843, 0.255, 0.655)
    c2 = (0.4157, 0.7608, 0.4118)
    c3 = (0.227, 0.090, 0.447)
    ax_fourier.plot(radii, m1_strengths, color=c1, linestyle='solid', marker='.', label=r'$m=1$', alpha=0.8)
    ax_fourier.plot(radii, m2_strengths, color=c2, linestyle='dashed', marker='x', label=r'$m=2$', alpha=0.8)
    ax_fourier.axvline(x=bar_radius, color=c3, linestyle='dashed', label='Bar Radius = {:.3f} pc'.format(bar_radius), alpha=1)
    ax_fourier.grid(color='grey', linestyle='dotted', alpha=0.5)
    ax_fourier.set_ylim(0.21, 1.03)
    if row == 3:
        ax_fourier.set_xlabel(r'$\rm Radius \, (pc)$', fontsize=14)
    else:
        ax_fourier.set_xticks([])
    ax_fourier.set_ylabel(r'$\rm Amplitude$', fontsize=12)
    if row == 0:
        ax_fourier.legend(loc='upper center', fontsize=11)

    # iterate over columns
    for column in range(ncols):

        grid_index = row * ncols + column
        ax = fig.add_axes([fourier_width + column*size, 0.8 - row*size, size, size])
        #ax = grid[grid_index]

        if column == 0:
            # density
            cmap = 'viridis'
            min_n=9e3
            max_n=4e9
            n = field_from_sliceplot("number_density", ds, disk, center, width_pc, north, dir, npixels=npixels)
            im1 = ax.imshow(n, cmap=cmap, origin="lower", norm=LogNorm())
            im1.set_clim(min_n, max_n)
        elif column == 1:
            # toomre q  
            cmap = 'cubehelix'
            cmap = plt.get_cmap(cmap)
            minval = 0.0
            maxval = 0.8
            n = 100
            cmap = mcolors.LinearSegmentedColormap.from_list('trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
            min_n=1
            max_n=1e3
            q = toomre_from_sliceplot(ds, disk, center, width_pc, north, dir, npixels=npixels)
            im2 = ax.imshow(q, cmap=cmap, origin="lower", norm=LogNorm())
            im2.set_clim(min_n, max_n)
        elif column == 2:
            # cylindrical_radial_velocity (add this part)
            cmap = cmyt.kelp  # divergine, 'coolwarm, 'rainbow'
            cmap = "magma"
            min_v = -13       
            max_v = 4    
            velocity = field_from_sliceplot("velocity_cylindrical_radius", ds, disk, center, width_pc, north, dir, npixels=npixels).to("km/s")
            im3 = ax.imshow(velocity, cmap=cmap, origin="lower")
            im3.set_clim(min_v, max_v)

        # Add annotations
        if grid_index == 0:
            at = AnchoredText(sim_label, loc='upper left', frameon=True, bbox_to_anchor=(0.01, 0.99), bbox_transform=ax.transAxes)
            at.txt._text.set_color("yellow")
            at.txt._text.set_fontsize(12)
            at.patch.set(boxstyle="square,pad=0.05", facecolor="black", linewidth=3, edgecolor="white", alpha=0.5)
            ax.add_artist(at)

        if column == 0:
            # Add BH mass annotation to first column
            at = AnchoredText(r"BH Mass: {:.0f} M$_\odot$".format(ss_mass.d), loc='lower center', frameon=False, bbox_transform=ax.transAxes)
            at.txt._text.set_color("white")
            at.txt._text.set_fontsize(12)
            ax.add_artist(at)

            # Add BH age annotation to leftmost y-labels
            #ax.set_ylabel("{:.2f} Myr".format((ss_age[0]/1e6)), fontsize=14)

        if column == 1:
            # Add BH Age annotation to second column
            at = AnchoredText(r"BH Age: {:.2f} Myr".format((ss_age[0]/1e6)), loc='lower right', frameon=False, bbox_transform=ax.transAxes)
            at.txt._text.set_color("black")
            at.txt._text.set_fontsize(12)
            ax.add_artist(at)

        # Set ticks and ticklabels
        ticks = np.linspace(-width_pc/2, width_pc/2, num=5)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        if row == 3:
            ax.set_xticklabels(tick_labels, fontsize=14)
            ax.set_xlabel("(pc)", fontsize=14)

        # Set tick parameters to create 'inner' ticks
        ax.tick_params(axis='both', direction='in', top=True, right=True, length=2, width=1, colors='black', grid_color='black', grid_alpha=0.5)
        ax.tick_params(which='minor', length=1)
    
        # Set Tick Locations (example with numpy's linspace, adjust as necessary)
        ax.set_xticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks
        ax.set_yticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks

# Adjust colorbars' positions and add a new one
cbar_ax1 = fig.add_axes([0.245, 0.97, 0.15, 0.01])
cbar_ax2 = fig.add_axes([0.40, 0.97, 0.15, 0.01])
cbar_ax3 = fig.add_axes([0.56, 0.97, 0.15, 0.01])

cbar1 = plt.colorbar(im1, cax=cbar_ax1, orientation='horizontal', ticklocation='top')
cbar2 = plt.colorbar(im2, cax=cbar_ax2, orientation='horizontal', ticklocation='top')
cbar3 = plt.colorbar(im3, cax=cbar_ax3, orientation='horizontal', ticklocation='top')

# Adding titles above colorbars
fig.text(0.05, 0.95, 'Fourier Modes', ha='center', va='center')
fig.text(0.32, 0.95, r'Number Density ($\rm cm^{-3}$)', ha='center', va='center')
fig.text(0.485, 0.95, r'Toomre $Q$', ha='center', va='center')
fig.text(0.64, 0.95, r'Radial Velocity ($\rm km/s$)', ha='center', va='center')  # Adjust this title as necessary

# save
plot_name = 'sliceplot-timeseries-' + str(sim_label) + '-' + str(width_pc) + 'pc-fourier_modes-toomreq-radial_vel_r=' + str(disc_r_pc) + 'pc.pdf'
plt.savefig('plots/' + plot_name, bbox_inches='tight')
print("created plots/" + str(plot_name))
plt.close()
