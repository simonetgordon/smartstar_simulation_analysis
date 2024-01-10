import yt
import cmyt
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from yt.utilities.math_utils import ortho_find
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib.offsetbox import AnchoredText
from plot_radial_profile_from_frb import extract_simulation_name, extract_dd_segment
from helper_functions import ss_properties
from helper_functions import _make_disk_L
from helper_functions import tidy_data_labels
from plot_toomre_q_projection import toomre_from_sliceplot, field_from_sliceplot
from helper_functions import calculate_theta_array, find_fourier_modes_and_phase_angles, find_bar_radius, find_pattern_speed_rad_per_sec


def main(root_dir, sim, dds_list):
    print("Loading datasets for simulation " + str(sim[0]))
    for dd_set, dds in enumerate(dds_list):
        print("Plotting dd set " + str(dd_set))
        DS = []
        for s in range(len(dds)):
            ds = yt.load(os.path.join(root_dir[0], sim[0], dds[s]))
            DS.append(ds)
        for disc_r_pc in [0.05]:

            #### PLOT ####

            fig = plt.figure(figsize=(12, 12))

            # IMAGE GRID #
            nrows = 4
            ncols = 3
            #grid = ImageGrid(fig, (0.1, 0.12, 0.60, 0.95), nrows_ncols=(nrows, ncols), axes_pad=0.0, share_all=True, cbar_mode=None)
            #grid = ImageGrid(fig, (0.25, 0.0, 0.75, 1.0), nrows_ncols=(nrows, ncols), axes_pad=0.0, share_all=True, cbar_mode=None)
            plt.rcParams['text.usetex'] = True

            # Plot data and annotations
            grid_index = 0
            size = 0.14 # subplot dims
            fourier_width = size*1.2
            for row in range(nrows):
                
                # Define angular speeds and fourier mode plot axes
                ax_fourier = fig.add_axes([0.02, 0.8 - row*size, fourier_width, size-0.01], frame_on=True)
                ax_speed = fig.add_axes([0.25, 0.8 - row*size, fourier_width, size-0.01], frame_on=True)

                # Load dataset and define axis
                ds = DS[row]
                sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
                sim_label = sim_label.replace("-2", "")
                sim_label = sim_label.replace("RS", "")
                print("Plotting " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))

                # Grab bh properties and define center, width and resolution of sliceplots
                ss_pos, ss_mass, ss_age = ss_properties(ds)
                center = ss_pos
                width_pc = 1
                tick_labels = ['', '-0.05', '0.0', '0.05', '']
                npixels = 2048
                theta = calculate_theta_array((npixels, npixels))
                dx = ds.index.get_smallest_dx().in_units('cm')
                dV = dx**3

                # Obtain angular momentum vector from small disk and define larger disk for plotting
                disc_h_pc = disc_r_pc
                _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
                vecs = ortho_find(L)
                dir = vecs[0]
                north = vecs[1]
                disc_r_pc_big = disc_h_pc_big = 0.8 # pc
                disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)

                # Obtain density and radius values for each cell in the disk
                density, radius_pc = field_from_sliceplot("density", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=True)
                
                # List of radii to define annular regions with thickness dr
                dr = 0.001 # pc
                r_min = 0.004
                r_max = 0.4
                radii = np.arange(r_min, r_max + dr, dr) # 73

                # Compute bar strength and phase angle variability across discrete annular regions
                cylindrical_velocity_theta = field_from_sliceplot("velocity_cylindrical_theta", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=False)
                m1_strengths, m2_strengths, _, phi_2_values, angular_speeds = find_fourier_modes_and_phase_angles(radii, radius_pc, density, theta, dV, dr, cylindrical_velocity_theta)

                ## FOURIER MODES ##

                # Find radius of the bar feature
                var_deg = 7.9 # degrees - 2 deg is too low, tiny bar
                bar_radius, i = find_bar_radius(phi_2_values, radii, var_deg=var_deg)

                # Plot bar strength and phase angle variability across discrete annular regions
                c1 = (1.0,0.2,0.5882)    # m1, rose bonbon
                c2 = (0.1882,0.6706,0.3647) # m2, pigment green
                c3 = (0.9608,0.3922,0.0627) # bar, orange pantone
                if row == 0:
                    ax_fourier.plot(radii, m1_strengths, color=c1, linestyle='solid', marker=None, label=r'$m=1$', alpha=0.8)
                    ax_fourier.plot(radii, m2_strengths, color=c2, linestyle='solid', marker=None, label=r'$m=2$', alpha=0.8)
                    ax_fourier.axvline(x=bar_radius, color=c3, linestyle='dashed', label=r'$R_{{\rm bar}}$ = {:.3f} pc'.format(bar_radius), alpha=1)
                    ax_fourier.legend(loc='upper right', fontsize=11, handlelength=1.4)
                else:
                    ax_fourier.plot(radii, m1_strengths, color=c1, linestyle='solid', marker=None, alpha=0.8)
                    ax_fourier.plot(radii, m2_strengths, color=c2, linestyle='solid', marker=None, alpha=0.8)
                    ax_fourier.axvline(x=bar_radius, color=c3, linestyle='dashed', label=r'$R_{{\rm bar}}$ = {:.3f} pc'.format(bar_radius), alpha=1) 
                    ax_fourier.legend(loc='upper right', fontsize=11, handlelength=1.4) if row != 1 else None
                
                if row == 3:
                    ax_fourier.tick_params(axis='x', direction='out', which='major')
                    ax_fourier.tick_params(axis='x', direction='out', which='minor')
                else:
                    ax_fourier.tick_params(axis='x', direction='in', which='major')
                    ax_fourier.tick_params(axis='x', direction='in', which='minor')
                
                # Set axes limits and ticks
                f_yticks = [0.4, 0.6, 0.8, 1.0]
                ax_fourier.set_ylim(0.23, 1.03)
                ax_fourier.set_yticks(f_yticks)
                ax_fourier.set_xlim(0, 0.11)
                if row == 3:
                    ax_fourier.set_xlabel(r'$\rm Radius \, (pc)$', fontsize=12)
                else:
                    ax_fourier.set_xlabel('')
                    ax_fourier.set_xticklabels([])
                
                #ax_fourier.set_ylabel(r'$\rm Amplitude$', fontsize=12, labelpad=40) if row == 0 else None
                ax_fourier.set_title('Fourier Modes') if row == 0 else None
                ax_fourier.grid(color='grey', linestyle='dotted', alpha=0.5)


                ## COROTATION FREQUENCY ##

                # find mean angular speed per annulus and plot
                angular_speed_means = [np.mean(annulus) for annulus in angular_speeds] # rad/sec
                sec_to_yr = 3.154e7*yt.units.s*yt.units.year # seconds in a year
                angular_freq_means_per_yr = [speed * sec_to_yr *2*np.pi for speed in angular_speed_means] # rad/yr
                label1 = r'Disc Angular Freq.' if row == 0 else None
                ax_speed.plot(radii, angular_freq_means_per_yr, color='teal', label=label1, linestyle='solid', marker=None)

                # find pattern speed and plot
                pattern_speeds_rad_per_sec = find_pattern_speed_rad_per_sec(ds, root_dir, sim, phi_2_values, ss_age, width_pc, disc_h_pc, disc_r_pc, disc_r_pc_big, disc_h_pc_big, npixels, radii, theta, dV)
                lim_ps = (np.abs(radii - bar_radius)).argmin() # limit the number of pattern speed points plotted to almost within the bar radius
                pattern_freq_rad_per_yr = [speed2 * sec_to_yr*2*np.pi for speed2 in pattern_speeds_rad_per_sec]
                ps = np.abs(pattern_freq_rad_per_yr)[0:lim_ps]
                mean_ps = np.mean(ps)
                label2 = r'$m = 2$ Pattern Freq.' if row == 0 else None
                ax_speed.plot(radii[0:lim_ps+5], np.abs(pattern_freq_rad_per_yr)[0:lim_ps+5], label=label2, linestyle=None, color='darkred', marker='x', linewidth=0.6)
                label3 = r'Mean Bar Pattern Freq.' if row == 0 else None
                ax_speed.axhline(y=mean_ps, label=label3, color='darkred', linestyle="-", alpha=0.5)

                # find corotation radius and plot
                corotation_radius = radii[(np.abs(angular_freq_means_per_yr - mean_ps)).argmin()]
                ax_speed.axvline(x=corotation_radius, label=r'$R_{{\rm co-rot}}$ = {:.3f} pc'.format(corotation_radius), color='goldenrod', linestyle="--", alpha=0.7)

                # Add legends (full legend in first row, R_co-rot only in rest)
                ax_speed.legend(fontsize=11, loc='upper center', handlelength=1.4) if row != 1 else None

                # set axes limits and labels
                ax_speed.set_xscale('log')
                ax_speed.set_yscale('log')
                ax_speed.grid(color='grey', linestyle='dotted', alpha=0.5)
                ax_speed.set_ylim(7e-7, 7e-3)
                if row == 3:
                    ax_speed.set_xlabel('Radius (pc)')
                else:
                    ax_speed.set_xlabel('')
                    ax_speed.set_xticklabels([])
                    ax_speed.tick_params(axis='x', which='major', direction='in')
                    ax_speed.tick_params(axis='x', which='minor', direction='in')
                ax_speed.set_title('Co-rotation Radius') if row == 0 else ax_speed.set_title('')
                #ax_speed.set_ylabel('Frequency (1/yr)', labelpad=40)  if row == 0 else None
                    

                ## SLICEPLOTS ##

                for column in range(ncols):
                    # Reset width_pc for the sliceplots
                    width_pc = 0.2

                    # Define axis for sliceplot
                    grid_index = row * ncols + column
                    ax = fig.add_axes([fourier_width*2.85 + column*(size-0.01), 0.8 - row*size, size-0.01, size-0.01])

                    if column == 0:
                        # density
                        cmap = 'viridis'
                        min_n=9e3
                        max_n=2e9
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
                        min_n=4e-1
                        max_n=4e2
                        q = toomre_from_sliceplot(ds, disk, center, width_pc, north, dir, npixels=npixels)
                        im2 = ax.imshow(q, cmap=cmap, origin="lower", norm=LogNorm())
                        im2.set_clim(min_n, max_n)
                    elif column == 2:
                        # cylindrical_radial_velocity (add this part)
                        cmap = cmyt.kelp  # diverging, 'coolwarm, 'rainbow'
                        cmap = "RdBu_r" # magma
                        min_v = -10     
                        max_v = 10   
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
                    pad = 0.02
                    if column == 0:
                        # Add BH Age annotation to second column
                        at = AnchoredText(r"BH Age: {:.2f} Myr".format((ss_age[0]/1e6)), loc='lower center', frameon=False, bbox_transform=ax.transAxes, pad=pad)
                        at.txt._text.set_color("white")
                        at.txt._text.set_fontsize(12)
                        ax.add_artist(at)

                    if column == 1:
                        # Add BH mass annotation to first column
                        at = AnchoredText(r"BH Mass: {:.0f} M$_\odot$".format(ss_mass.d), loc='lower center', frameon=False, bbox_transform=ax.transAxes, pad=pad)
                        at.txt._text.set_color("black")
                        at.txt._text.set_fontsize(12)
                        ax.add_artist(at)

                    if column == 2:
                        size_in_data_units = 205  # Adjust this value based on your data in pixel units
                        label = "0.02 pc"
                        location=4 # lower right
                        scale_color = 'black'
                        bar = AnchoredSizeBar(ax.transData, size_in_data_units, label, location, pad=pad, color=scale_color, frameon=False)
                        ax.add_artist(bar)

                    # Set ticks and ticklabels
                    ax.set_xticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks
                    ax.set_yticks(np.linspace(0, npixels, num=5))  # Adjust num for number of ticks
                    ticks = np.linspace(-width_pc/2, width_pc/2, num=5)
                    ax.set_xticklabels([])
                    ax.set_yticklabels([])
                    if column == 0:
                        ax.set_yticklabels(tick_labels)
                        ax.set_ylabel("(pc)")
                    if row == 3:
                        ax.set_xticklabels(tick_labels)
                        ax.set_xlabel("(pc)")

                    # Set tick parameters to create 'inner' ticks, except on the bottom axis of last row
                    ax.tick_params(which='minor', length=1)
                    ax.tick_params(axis='x', which='both', direction='inout', bottom=True, top=True, length=2, width=1, colors='black')

            # Adjust colorbars' positions and add a new one
            fourier_width *= 1.06
            cbar_ax1 = fig.add_axes([0.301+fourier_width, 0.95, size-0.015, 0.01])
            cbar_ax2 = fig.add_axes([0.432+fourier_width, 0.95, size-0.015, 0.01])
            cbar_ax3 = fig.add_axes([0.563+fourier_width, 0.95, size-0.015, 0.01])

            plt.colorbar(im1, cax=cbar_ax1, orientation='horizontal', ticklocation='top')
            plt.colorbar(im2, cax=cbar_ax2, orientation='horizontal', ticklocation='top')
            plt.colorbar(im3, cax=cbar_ax3, orientation='horizontal', ticklocation='top')

            # Adding titles above colorbars
            title_height = 0.94
            #fig.text(0.125, title_height, r'Fourier Modes', ha='center', va='center')
            fig.text(0.365+fourier_width, title_height, r'Number Density ($\rm cm^{-3}$)', ha='center', va='center')
            fig.text(0.493+fourier_width, title_height, r'Toomre $Q$', ha='center', va='center')
            fig.text(0.625+fourier_width, title_height, r'Radial Velocity ($\rm km/s$)', ha='center', va='center')  # Adjust this title as necessary

            # Add ylabels
            # Assuming your figure is named `fig` and you have your subplots in `ax` array
            fig.text(-0.025, 0.665, 'Amplitude $A_m/A_0$', va='center', rotation='vertical', fontsize=12)
            fig.text(0.02+fourier_width, 0.67, 'Frequency (1/yr)', va='center', rotation='vertical', fontsize=12)

            # save
            plot_name = f'sliceplot-timeseries-{sim_label}-{width_pc}pc-fourier_modes-toomreq-radial_vel_r={disc_r_pc}pc_ageset_low_toomre_{var_deg}deg_dr={dr}_rmin={r_min}.pdf'
            plt.savefig('plots/' + plot_name, bbox_inches='tight')
            print("created plots/" + str(plot_name))
            plt.close()

if __name__ == "__main__":
    """
    Plot a timeseries of sliceplots for a given simulation, showing the density, toomre q and radial velocity fields.
    To run: python plot_density_slices_toomreq.py
    """

    root_dir = [#"/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/",
            #"/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/"
            #"/Backup00/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/",
            "/Backup00/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/"
            ]
    sim = [#"2B.RSm04", 
        #"2B.m08-4dx"
        #"1B.m16-4dx"
        "2B.RSb08-2"
        ] 

    dds = ["DD0167/DD0167", "DD0178/DD0178", "DD0189/DD0189", "DD0231/DD0231"]  # 0.39, 0.50, 0.6, 1 Myr for 1B.m16,
    dds = ["DD0178/DD0178", "DD0189/DD0189", "DD0199/DD0199", "DD0225/DD0225"]  # 0.50, 0.60, 0.69, 0.94 Myr for 1B.m16,
    #dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
    #dds = ["DD0219/DD0219", "DD0227/DD0227", "DD0236/DD0236", "DD0279/DD0279"]  # 0.2, 0.69, 1 Myr for 2B.b08,

    # 1B.m16
    dds_list = [["DD0178/DD0178", "DD0188/DD0188", "DD0199/DD0199", "DD0226/DD0226"],
                #["DD0178/DD0178", "DD0189/DD0189", "DD0199/DD0199", "DD0225/DD0225"],
                # ["DD0179/DD0179", "DD0186/DD0186", "DD0200/DD0200", "DD0228/DD0228"],
                # ["DD0177/DD0177", "DD0188/DD0188", "DD0197/DD0197", "DD0226/DD0226"],
                ]

    #2B.b08
    dds_list = [#["DD0219/DD0219", "DD0227/DD0227", "DD0236/DD0236", "DD0276/DD0276"],
                #["DD0218/DD0218", "DD0226/DD0226", "DD0235/DD0235", "DD0278/DD0278"],
                #["DD0220/DD0220", "DD0228/DD0228", "DD0237/DD0237", "DD0279/DD0279"],
                #["DD0217/DD0217", "DD0225/DD0225", "DD0234/DD0234", "DD0270/DD0270"],
                ["DD0219/DD0219", "DD0246/DD0246", "DD0270/DD0270", "DD0279/DD0279"], # 0.29, 0.50, 0.70, 1 Myr
                ]

    main(root_dir, sim, dds_list)
