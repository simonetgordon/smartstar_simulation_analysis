import os
from plot_toomre_q_projection import field_from_sliceplot
from helper_functions import ss_properties
from helper_functions import _make_disk_L
from helper_functions import tidy_data_labels
from yt.utilities.math_utils import ortho_find
from plot_radial_profile_from_frb import extract_dd_segment, extract_simulation_name
from find_fourier_modes import get_theta_values, find_bar_radius, extract_dd_number
from plot_density_slices_toomreq import find_fourier_modes_and_phase_angles
import numpy as np
import yt
from multiprocessing import Pool
from mpi4py import MPI
import ast
#import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def find_corotation_radius(radii, radius_pc, velocity_cylindrical_theta, omega_pattern):
    """
    Find the corotation radius where the angular speed Ω(r) equals the pattern speed Ω_pattern.
    
    Input:
        radii: 1D array of radii for the disk (annuli)
        radius_pc: 2D array of radii values for each cell in the disk
        velocity_cylindrical_theta: 2D array of circular velocities for each cell in the disk
        omega_pattern: pattern speed
    
    Output:
        corotation_radius: The radius at which Ω(r) ≈ Ω_pattern
    """
    # Calculate the angular speed Ω(r) for each cell in units 1/s
    angular_speeds = velocity_cylindrical_theta / radius_pc.to("cm")

    # Initialize variables to store the minimum difference and corresponding corotation radius
    min_diff = np.inf
    corotation_radius = None

    # Iterate over the array of annuli to find where Ω(r) - Ω_pattern is minimized
    for r in radii:
        # Define mask for cells at this radius
        mask = np.isclose(radius_pc, r, atol=0.5) # Adjust atol as per your discretization

        # Calculate the mean angular speed at this radius
        mean_angular_speed = np.mean(angular_speeds[mask])

        # Calculate the difference between Ω(r) and Ω_pattern
        diff = np.abs(mean_angular_speed - omega_pattern)

        # Update the minimum difference and corotation radius if this is the smallest difference so far
        if diff < min_diff:
            min_diff = diff
            corotation_radius = r

    return corotation_radius


def process_ds(ds, disc_r_pc, find_cylindrical_velocity=False):
    disc_r_pc = 0.1 # pc
    find_cylindrical_velocity = True
    # identify simulation name and label
    sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
    print("Processing " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))

    # [Your code for setting up disk, radii, etc.]
    # Grab bh properties and define center, width and resolution of sliceplots
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    center = ss_pos
    width_pc = 0.2
    npixels = 2048

    # Obtain angular momentum vector from small disk and define larger disk for plotting
    disc_h_pc = disc_r_pc
    _, L = _make_disk_L(ds, center, disc_r_pc, disc_h_pc)
    vecs = ortho_find(L)
    dir = vecs[0]
    north = vecs[1]
    disc_r_pc_big = disc_h_pc_big = 0.6 # pc
    disk = ds.disk(center, L, disc_r_pc_big, disc_h_pc_big)

    # Obtain density and radius values for each cell in the disk
    density, radius_pc = field_from_sliceplot("density", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=True)
    dx = ds.index.get_smallest_dx().in_units('pc')
    surface_density = density * dx # g/cm^2

    # List of radii to define annular regions with thickness dr
    dr = 0.001 # pc
    dx = ds.index.get_smallest_dx().in_units('pc').d
    dx = 0.0008
    r_max_pc = 0.14 # pc 
    r_min = dx
    r_max = r_max_pc
    radii = np.arange(r_min, r_max + dr, dr)  # 141 annuli (fixed at this number)

    # Compute theta values of each annulus
    dV = dx**3
    theta = get_theta_values(surface_density)
    cylindrical_velocity_theta = field_from_sliceplot("velocity_cylindrical_theta", ds, disk, center, width_pc, north, dir, npixels=npixels, radius=False)[0]

    # Compute Fourier modes and phase angles in degrees (not radians), and store them
    _, _, _, phi_2_values, angular_speeds = find_fourier_modes_and_phase_angles(radii, radius_pc, density, theta, dV, dr, cylindrical_velocity_theta)
        
    if find_cylindrical_velocity:
        return phi_2_values, ss_age/1e6, radii, angular_speeds

    return phi_2_values, ss_age/1e6, radii


if __name__ == "__main__":
    import pandas as pd
    # Set paths and parameters
    root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/270msun/replicating-beckmann/"
    sim = "1B.m16-4dx"
    dds = [f"DD0{i:03d}/DD0{i:03d}" for i in range(190, 200)]
    disc_r_pc = 0.1  # pc

    # MPI initialization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size() 

    # Load datasets (this can be optimized if datasets are large)
    DS = []
    for s in range(len(dds)):
        ds = yt.load(os.path.join(root_dir[0], sim[0], dds[s]))
        DS.append(ds)

    # Set disc radius from which to calculate L
    disc_r_pc = 0.1 # pc
    phi_2_values = [] # degrees
    ages = []
    radii = [] # pc
    angular_speeds = [] # 1/s
    for i in range(rank, len(dds), size):
        ds = yt.load(os.path.join(root_dir, sim, dds[i]))
        phi_2_value, current_age, radiuses, angular_speed = process_ds(ds, disc_r_pc, find_cylindrical_velocity=True)
        phi_2_values.append(phi_2_value)
        ages.append(current_age)
        radii.append(radiuses)
        angular_speeds.append(angular_speed) 

    # Gather results at root
    all_phi_2_values = comm.gather(phi_2_values, root=0)
    all_ages = comm.gather(ages, root=0)
    all_angular_speeds = comm.gather(angular_speeds, root=0)
    all_radii = comm.gather(radii, root=0)

    # Only the root process will output the results
    if rank == 0:

        # Define all variables after MPI gather
        phi_2_values_list = [np.array(phi) for phi in all_phi_2_values] # list of arrays
        ages_list = np.array([age[0] for age in all_ages]) # list of floats
        angular_speed_means = [[np.mean(annulus) for annulus in dataset] for dataset in angular_speeds]

        #  Compute the difference between phi2 values at each snapshot
        delta_phi = np.array([phi_2_values_list[i] - phi_2_values_list[i+1] for i in range(len(phi_2_values_list)-1)])
        delta_t = np.array(np.diff(ages_list))
        pattern_speeds = delta_phi / delta_t[:, np.newaxis] # deg/Myr

        # plot bar pattern and disc angular speed vs radius
        fig, ax = plt.subplots(figsize=(8, 5))
        n = len(angular_speed_means)
        # Choose a colormap
        cmap = plt.cm.rainbow
        # Create color list
        colors = [cmap(i) for i in np.linspace(0, 1, n)]
        scalars = np.linspace(all_ages[0][0], all_ages[-1][0], n)
        norm = plt.Normalize(min(scalars), max(scalars))
        # Create a ScalarMappable and set array to the normalized values
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        # Plot angular speed and pattern speed vs radius
        for j in range(len(angular_speed_means)):
            plt.plot(all_radii[j], angular_speed_means[j], label=f"{all_ages[j][0]:.2f} Myr", linestyle="-", color=colors[j])
        conversion_factor = np.pi / (180 * 3.154e13)
        pattern_speeds_radians_per_sec = pattern_speeds * conversion_factor
        for i in range(len(pattern_speeds)):
            plt.plot(all_radii[i+1], abs(pattern_speeds_radians_per_sec)[i], linestyle="--", color=colors[i])
        handles, labels = ax.get_legend_handles_labels()
        # Assuming you want to include only the first two
        legend_handle_1 = handles[0]
        legend_handle_2 = handles[-1]
        plt.legend(handles=[legend_handle_1, legend_handle_2], labels=[labels[0], labels[-1]])
        plt.yscale("log")    
        plt.xlabel("Radius (pc)")
        plt.ylabel("Speed (rad/s)")
        plt.title("Mean Angular Speed and Mean Pattern Speed vs Radius")
        plt.savefig(f"corotation_radius/angular_and_pattern_speed_vs_radius_{all_ages[0][0]:.2f}_{all_ages[-1][0]:.2f}.png", bbox_inches="tight") 
        plt.close()
