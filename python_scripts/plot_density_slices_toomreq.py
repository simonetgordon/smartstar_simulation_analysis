import yt
import cmyt
import shutil
import sys
import os
import re # complex str searches
from matplotlib import ticker
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from plot_multipanel_time_2 import configure_font, make_projection_plot, create_axes_grid, set_axis_labels_and_colorbar, \
    configure_projection_plot, get_min_max_values, set_ticks_and_labels, make_slice_plot
from plot_multi_projections import tidy_data_labels, find_north_vector
from yt.utilities.math_utils import ortho_find
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import AxesGrid
import time
import numpy as np
import matplotlib as mpl
from matplotlib.colors import LogNorm
from plot_toomre_q_projection import ToomreQ, kappa2D, toomre_from_sliceplot, field_from_sliceplot
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.offsetbox import AnchoredText
from plot_radial_profile_from_frb import extract_simulation_name, extract_dd_segment, setup_plot_env
import matplotlib.colors as mcolors


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
grid = ImageGrid(fig, (0.1, 0.12, 0.480, 0.95), nrows_ncols=(4, 3), axes_pad=0.0, share_all=True, cbar_mode=None)
plt.rcParams['text.usetex'] = True

# Plot data and annotations
grid_index = 0
for row in range(4):
    for column in range(3):

        # Load dataset and define axis
        ds = DS[row]
        sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
        sim_label = sim_label.replace("-2", "")
        sim_label = sim_label.replace("RS", "")
        print("Plotting " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))
        grid_index = row * 3 + column
        ax = grid[grid_index]

        # Grab bh properties and define center, width and resolution of sliceplots
        ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
        center = ss_pos
        width_pc = 0.32
        tick_labels = ['', '-0.08', '0.0', '0.08', '']
        npixels = 2048

        # Obtain angular momentum vector from small disk and define larger disk for plotting
        disc_r_pc = disc_h_pc = 0.01
        _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
        vecs = ortho_find(L)
        dir = vecs[0]
        north = vecs[1]
        disc_r_pc_big = disc_h_pc_big = 1.0
        disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)

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
            #cmap='magma'
            #cmap = cmyt.octarine
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
cbar_ax1 = fig.add_axes([0.102, 0.885, 0.15, 0.01])
cbar_ax2 = fig.add_axes([0.265, 0.885, 0.15, 0.01])
cbar_ax3 = fig.add_axes([0.425, 0.885, 0.15, 0.01])

cbar1 = plt.colorbar(im1, cax=cbar_ax1, orientation='horizontal', ticklocation='top')
cbar2 = plt.colorbar(im2, cax=cbar_ax2, orientation='horizontal', ticklocation='top')
cbar3 = plt.colorbar(im3, cax=cbar_ax3, orientation='horizontal', ticklocation='top')

# Adding titles above colorbars
fig.text(0.18, 0.872, r'Number Density ($\rm cm^{-3}$)', ha='center', va='center')
fig.text(0.34, 0.872, r'Toomre $Q$', ha='center', va='center')
fig.text(0.50, 0.872, r'Radial Velocity ($\rm km/s$)', ha='center', va='center')  # Adjust this title as necessary

# save
plot_name = 'sliceplot-timeseries-' + str(sim_label) + '-' + str(width_pc) + 'pc-toomreq-radial_vel_r=' + str(disc_r_pc) + 'pc.pdf'
plt.savefig('plots/' + plot_name, bbox_inches='tight')
print("created plots/" + str(plot_name))
plt.close()


# if __name__ == "__main__":

#     root_dir = [#"/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
#                 #"/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/"
#                 "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/"

#                 ]
#     sim = [#"2B.RSm04", 
#            #"2B.m08-4dx"
#            "1B.m16-4dx"
#            ] 
#     dds1=None
#     dds2 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.29, 079, 1.01 Myr for m01, - always need to define dds2
#     #xdds3 = ["DD0201/DD0201", "DD0206/DD0206", "DD0208/DD0208"]  # 0.201, 0.801, 1.002 Myr for b01
#     #dds2 = ["DD0229/DD0229", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.8 Myr for 2B.m04, 
#     dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
#     dds2 = ["DD0167/DD0167", "DD0178/DD0178", "DD0188/DD0188"]  # 0.3, 0.69, 0.8 Myr for 1B.m16,
#     fontsize = 14
#     width_pc = 0.18 # must change xticks if this changes
#     xticks = [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
#     xticks = [-0.1, -0.05, 0.0, 0.05, 0.1]
#     xticks = [-0.05, 0.0, 0.05]
#     min_n=2e5
#     max_n=0.1
#     field = "number_density"
#     cmap = 'viridis'
#     slice = True

#     main(root_dir, sim, dds1, dds2, dds3, field, width_pc, xticks, fontsize, min_n, max_n, cmap, slice)