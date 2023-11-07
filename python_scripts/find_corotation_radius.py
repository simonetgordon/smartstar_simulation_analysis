import os
from plot_toomre_q_projection import field_from_sliceplot
from smartstar_find import ss_properties
from plot_disc_projections import _make_disk_L
from plot_multi_projections import tidy_data_labels
from yt.utilities.math_utils import ortho_find
from plot_radial_profile_from_frb import extract_dd_segment, extract_simulation_name
from find_fourier_modes import get_theta_values, find_bar_radius, extract_dd_number
from plot_density_slices_toomreq import find_fourier_modes_and_phase_angles
import numpy as np
import yt
from multiprocessing import Pool
from mpi4py import MPI
import pandas as pd
import ast


def process_ds(ds, disc_r_pc, find_cylindrical_velocity=False):
    # identify simulation name and label
    sim_label = tidy_data_labels(extract_simulation_name(ds.directory))
    print("Processing " + str(sim_label) + " " + str(extract_dd_segment(ds.directory)))

    # [Your code for setting up disk, radii, etc.]
    # Grab bh properties and define center, width and resolution of sliceplots
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    center = ss_pos
    width_pc = 0.2
    npixels = 2048
    dx = ds.index.get_smallest_dx().in_units('cm')

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
    surface_density = density * dx # g/cm^2
    
    # List of radii to define annular regions with thickness dr
    dr = 0.001 # pc
    r_min = np.min(radius_pc).value
    r_max = np.max(radius_pc).value
    radii = np.arange(r_min, r_max + dr, dr) 

    # Compute theta values of each annulus
    dV = dx**3
    theta = get_theta_values(surface_density)

    # Compute cylindrical velocity component

    # Compute Fourier modes and phase angles in degrees (not radians), and store them
    m1_strengths, m2_strengths, phi_1_values, phi_2_values = find_fourier_modes_and_phase_angles(radii, radius_pc, density, theta, dV, dr)

    if find_cylindrical_velocity:
        # Compute cylindrical velocity component
        cylindrical_velocity_theta = field_from_sliceplot("velocity_cylindrical_theta", ds, disk, center, width_pc, north, dir, npixels=npixels)

        return phi_2_values, ss_age/1e6, radii, cylindrical_velocity_theta

    return phi_2_values, ss_age/1e6, radii


if __name__ == "__main__":
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
    dds = ["DD0178/DD0178", "DD0189/DD0189", "DD0199/DD0199", "DD0225/DD0225"]  # 0.50, 0.60, 0.69, 0.94 Myr for 1B.m16,
    #dds3 = ["DD0228/DD0228", "DD0268/DD0268", "DD0280/DD0280"]  # 0.3, 0.69, 0.79 Myr for 2B.m08-4dx, 
    #dds = ["DD0219/DD0219", "DD0227/DD0227", "DD0236/DD0236", "DD0279/DD0279"]  # 0.2, 0.69, 1 Myr for 2B.b08,

    # 1B.m16
    dds_list = [["DD0178/DD0178", "DD0188/DD0188", "DD0199/DD0199", "DD0226/DD0226"],
                #["DD0178/DD0178", "DD0189/DD0189", "DD0199/DD0199", "DD0225/DD0225"],
                # ["DD0179/DD0179", "DD0186/DD0186", "DD0200/DD0200", "DD0228/DD0228"],
                # ["DD0177/DD0177", "DD0188/DD0188", "DD0197/DD0197", "DD0226/DD0226"],
                ]
    # 1B.m16 for bar corotation radius
    dds = [f"DD0{i:03d}/DD0{i:03d}" for i in range(190, 220)] # 0.6 - 0.9 Myr

    #2B.b08
    # dds_list = [#["DD0219/DD0219", "DD0227/DD0227", "DD0236/DD0236", "DD0276/DD0276"],
    #             #["DD0218/DD0218", "DD0226/DD0226", "DD0235/DD0235", "DD0278/DD0278"],
    #             #["DD0220/DD0220", "DD0228/DD0228", "DD0237/DD0237", "DD0279/DD0279"],
    #             #["DD0217/DD0217", "DD0225/DD0225", "DD0234/DD0234", "DD0270/DD0270"],
    #             ["DD0219/DD0219", "DD0246/DD0246", "DD0270/DD0270", "DD0279/DD0279"], # 0.29, 0.50, 0.70, 1 Myr
    #             ]

    filename = "phi2_veltheta_over_time.csv"
    if not os.path.exists(filename):
        print(f"Filling in {filename}...")
        # Load datasets
        DS = []
        for s in range(len(dds)):
            ds = yt.load(os.path.join(root_dir[0], sim[0], dds[s]))
            DS.append(ds)

        # Set disc radius from which to calculate L
        disc_r_pc = 0.1 # pc

        # MPI initialization
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        # Distribute datasets among ranks and process them (find the phi2 value array for each dataset)
        phi_2_values = []
        ages = []
        radii = []
        velocity_cylindrical_theta = []
        for i in range(rank, len(DS), size): # starting from 'rank', up to total number of datasets with a step of core number
            ds = DS[i]
            phi_2_value, current_age, radiuses, velocity_theta = process_ds(ds, disc_r_pc, find_cylindrical_velocity=True)
            phi_2_values.append(phi_2_value)
            ages.append(current_age)
            radii.append(radiuses)
            velocity_cylindrical_theta.append(velocity_theta)

        # Gather results at root
        all_phi_2_values = comm.gather(phi_2_values, root=0)
        all_times = comm.gather(ages, root=0)
        all_velocity_cylindrical_theta = comm.gather(velocity_cylindrical_theta, root=0)

        # Only the root process will output the results
        if rank == 0:
            # Flatten lists
            flat_phi_2_values = [item for sublist in all_phi_2_values for item in sublist]
            flat_times = [item[0] for sublist in all_times for item in sublist] # assuming the time arrays always have one element
            flat_velocity_cylindrical_theta = [item for sublist in all_velocity_cylindrical_theta for item in sublist]

            # Pair times with corresponding phi_2_values and sort
            paired_data = sorted(zip(flat_times, flat_phi_2_values, flat_velocity_cylindrical_theta), key=lambda x: x[0])

            # Write sorted results to file
            with open(filename, 'w') as f:
                f.write("Times\t Phi2_Values\t  Velocity_Cylindrical_Theta\n")
                for time, phi_2_value,velocity_cylindrical_theta in paired_data:
                    f.write(f"{time}\t{phi_2_value}\t{velocity_cylindrical_theta}\n")

            print(f"Data written to {filename}")

    print(f"Reading data from {filename}")
    df = pd.read_csv(filename, sep="\t")

    # Extract the columns into variables
    ages = np.array([float(age) for age in df["Times"].tolist()])
    # Convert the string representation of lists to actual lists
    phi2_values_list = df['Phi2_Values'].apply(ast.literal_eval).tolist()

    # Convert the lists to numpy arrays
    phi2_values_array = [np.array(sublist) for sublist in phi2_values_list]

    # Compute pattern speed
    def compute_difference(arr1, arr2):
        if len(arr1) > len(arr2):
            arr2 = np.pad(arr2, (0, len(arr1) - len(arr2)))
        elif len(arr1) < len(arr2):
            arr1 = np.pad(arr1, (0, len(arr2) - len(arr1)))
        return arr2 - arr1

    delta_phi = [compute_difference(phi2_values_array[i], phi2_values_array[i+1]) for i in range(len(phi2_values_array)-1)] # Difference along snapshots
    delta_t = np.diff(ages)
    pattern_speeds = delta_phi / delta_t

    # Use the mean or median pattern speed across all radii and time intervals, or choose an appropriate radii range
    shapes = [arr.shape for arr in pattern_speeds]
    unique_shapes = set(shapes)
    print(unique_shapes)

    # Find the maximum length
    max_length = max(shapes, key=lambda x: x[0])[0]

    # Pad the arrays in pattern_speeds
    padded_pattern_speeds = np.array([np.pad(arr, (0, max_length - len(arr))) if len(arr) < max_length else arr for arr in pattern_speeds])

    # find linear speeds
    # List of radii to define annular regions with thickness dr
    dr = 0.001 # pc
    dx = 0.00077
    r_min = dx
    r_max = 0.14
    radii = np.arange(r_min, r_max + 6*dr, dr) 
    radii_array = np.tile(radii, (29, 1))
    linear_speeds = padded_pattern_speeds*radii_array*np.pi*3.086e16/(180*1e6) # km/s

    # Now compute the mean in a certain time interval
    start_time = 0.7
    end_time = 0.8

    def mean_pattern_speed_within_interval(start, end, ages, pattern_speeds):
        # Filter the data
        
        if len(pattern_speeds) != len(ages):
            # This will truncate or extend the padded_pattern_speeds list to match the length of ages
            pattern_speeds = pattern_speeds[:len(ages)]

        indices_within_range = np.where((ages >= start_time) & (ages <= end_time))[0]
        filtered_pattern_speeds = pattern_speeds[indices_within_range]

        # Compute the mean
        mean_pattern_speed = np.mean(filtered_pattern_speeds)
        return mean_pattern_speed
    
    mean_pattern_speed = mean_pattern_speed_within_interval(start_time, end_time, ages, padded_pattern_speeds)
    mean_linear_speed = mean_pattern_speed_within_interval(start_time, end_time, ages, linear_speeds)
    print(f"Mean pattern speed between time {start_time} and {end_time} Myr: {mean_pattern_speed}")
    print(f"Mean linear speed between time {start_time} and {end_time} Myr: {mean_linear_speed}")
