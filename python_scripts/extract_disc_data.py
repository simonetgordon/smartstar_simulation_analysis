import yt
import logging
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from helper_functions import *
from helper_functions import _make_disk_L, ss_properties

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

# Define the function to process a single dataset
def process_disk_data(filepath, sim_name=None, min_avg_density=1e6, max_disc_r_pc=1.0, step_r_pc=0.1):
    logging.info(f"Started processing {filepath}")

    # Load the dataset
    ds = yt.load(filepath)
    sim = extract_specific_part(filepath) if sim_name is None else sim_name
    dd = extract_dd_segment(filepath)
    logging.info(f"Simulation: {sim} and DD: {dd} at time: {ds.current_time.to('Myr')}")

    ds.add_field(("gas", "thermal_energy"), function=specific_thermal_energy, units="erg/g", sampling_type="cell")
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    logging.info(f"BH properties: position = {ss_pos}, mass = {ss_mass:.2f}, age = {ss_age[0]/1e6:.2f} Myr")
    
    # Initial disk properties
    center = ss_pos
    disc_h_pc = 0.1  # Disk height in parsec
    disc_r_pc = 0.5  # Initial disk radius in parsec
    avg_density = 0  # Placeholder for average density
    
    # Function to determine the disk orientation
    _, L = _make_disk_L(ds, center, disc_r_pc/2, disc_h_pc)

    # Iterate to find the appropriate disk radius with the desired average density
    while avg_density < min_avg_density and disc_r_pc <= max_disc_r_pc:
        # Define the disk region
        disk = ds.disk(center, L, (disc_r_pc, 'pc'), (disc_h_pc, 'pc'))
        
        # Calculate the average density and total mass
        avg_density = disk.mean('number_density', weight='cell_volume').in_units('1/cm**3')
        logging.info(f"At radius {disc_r_pc} pc, avg_density = {avg_density}")
        
        # Increment the disk radius for the next iteration if density threshold is not met
        if avg_density < min_avg_density:
            logging.info(f"Disk radius step is decreased by {step_r_pc} pc")
            disc_r_pc -= step_r_pc
    
    # Extract time information
    time_myr = ds.current_time.in_units('Myr')
    
    # Get the disk radius (r_max) in pc
    disk_radius = disc_r_pc
    total_mass = disk.quantities.total_mass().in_units('Msun')
    
    # Return the computed values
    return {
        'simulation': filepath.split('/')[-3],  # change index based on folder structure
        'age_myr': time_myr.value,
        'datadump': dd,
        'disc_radius': disk_radius,
        'disc_avg_density': avg_density.value,
        'disc_mass': total_mass.value
    }

if __name__ == "__main__":

    filepaths = [
        # 1B.resim.th.b04 - end event 2.80
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0153/DD0153", # 2.500
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0177/DD0177", # 2.525
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0203/DD0203", # 2.550
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0227/DD0227", # 2.575
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0253/DD0253", # 2.600
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0277/DD0277", # 2.625
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0303/DD0303", # 2.651
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0327/DD0327", # 2.675
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0353/DD0353", # 2.700
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0377/DD0377", # 2.725
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0403/DD0403", # 2.751
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0427/DD0427", # 2.775
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0453/DD0453", # 2.800
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0477/DD0477", # 2.825
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0503/DD0503", # 2.851
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0527/DD0527", # 2.875
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0552/DD0552", # 2.900
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0577/DD0577", # 2.925
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0602/DD0602", # 2.950
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0627/DD0627", # 2.975
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0652/DD0652", # 3.000
        # # 1B.resim.th.b04-eta-0.1
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0153/DD0153", # 2.501
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0177/DD0177", # 2.525
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0202/DD0202", # 2.550
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0227/DD0227", # 2.575
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0252/DD0252", # 2.600
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0277/DD0277", # 2.625
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0302/DD0302", # 2.650
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0327/DD0327", # 2.675
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0352/DD0352", # 2.700
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0377/DD0377", # 2.725
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0402/DD0402", # 2.750
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0408/DD0408", # 2.775
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0411/DD0411", # 2.800
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0413/DD0413", # 2.825
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0416/DD0416", # 2.850
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0418/DD0418", # 2.875
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0421/DD0421", # 2.900
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0423/DD0423", # 2.925
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0426/DD0426", # 2.950
        # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0428/DD0428", # 2.975
    ]

    # Initialize an empty list to store the data
    disk_data = []

    # Loop through each filepath and process the data
    for filepath in filepaths:
        data = process_disk_data(filepath)
        disk_data.append(data)

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(disk_data)

    # Save the DataFrame to a CSV file
    filename = "csv_data_plotting/data_misc/disc_data.csv"
    df.to_csv(filename, index=False)

    print(f"Data has been processed and saved to {filename}")