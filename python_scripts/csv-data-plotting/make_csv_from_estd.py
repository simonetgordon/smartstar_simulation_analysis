import os
import re # complex str searches
import numpy as np
import csv
from itertools import zip_longest
import pandas as pd # for sorting chronologically

##########################################################################################################
#                                    Make csv file of BH/gas properties
#
# to run: python make_csv_from_estd.py
##########################################################################################################

def _remove_strings(lst):
    for index, i in enumerate(lst):
        try:
            float(i)
        except ValueError:
            i = i.replace("inf pc \t Jea", "0")
            i = i.replace("inf pc \t Bon", "0")
            i = i.replace("cm^-3", '')
            i = i.replace("cm^-", '')
            i = i.replace("cm^", '')
            i = i.replace("cm", '')
            i = i.replace("c", '')
            i = i.replace(",", '')
            i = i.replace("y", '')
            float(i)
        lst[index] = float(i)
    return lst


def extract_simulation_name(filepath):
    # Get the last part of the path
    last_part = os.path.basename(filepath)

    # Use regular expression to extract the full simulation name
    match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', last_part)

    if match:
        return match.group(0)

    # If the match is not found, try to extract from the parent directories
    path_parts = filepath.split(os.path.sep)
    for i in range(len(path_parts)-1, -1, -1):
        match = re.search(r'\b(?:\d+[A-Za-z]+\d+|[A-Za-z]+\d+)\b', path_parts[i])
        if match:
            return path_parts[i]
    return None


def get_data_arrays(root_dir):
    # Get list of estd*.out files in root directory
    estds = [f for f in os.listdir(root_dir) if f.endswith(".out") and f.startswith("estd")]

    # Collect all estd*.out files
    out_files = []
    misc_files = []
    for file in os.listdir(root_dir):
        if file.endswith(".out"):
            if file.startswith("estd_"):
                out_files.append(file)
            else:
                misc_files.append(file)

    # Sort .out files by the 1 or 2 digit number in their filename
    estds = sorted(out_files, key=lambda x: int(x.split("_")[1].split(".")[0]) 
                if len(x.split("_")[1].split(".")[0]) <= 2 else float('inf'))
    # Append estd.out file to end
    estds.extend(misc_files)

    # Extract the string after "seed" and before the next "-" character
    match = re.search(r"seed(\d+)-", root_dir)

    # Extract sim name as the penultimate dir in the filepath
    sim = extract_simulation_name(root_dir)
    print("Creating combined .csv for simulation ", sim)

    ##########################################################################################################
    #                                  Read from simulation output txt file(s)
    ##########################################################################################################

    output_combined = 'estd_files/' + sim

    with open(output_combined, "w") as outfile:
        # Append contents of each .out file to output file
        for out_file in estds:# writing data arrays to this file
            with open(root_dir + out_file, "r") as infile:
                outfile.write(infile.read())

    print("Created combined .csv for simulation ", sim)

    ##########################################################################################################
    #                                  Accretion Arrays: accrates, ages
    ##########################################################################################################
    accrates = []
    ages = []
    for line in open(output_combined):
        accrate = re.search(r'accrate = (.{8})', line)
        age = re.search(r'Age = (.{12})', line)
        if accrate:
            accrates.append(accrate.group(1))
            ages.append(age.group(1))

    accrates = np.array([float(i) for i in accrates])
    ages = np.array([float(i) for i in ages])


    ##########################################################################################################
    #                        Average Velocities, Densities + Total Gas Mass
    ##########################################################################################################

    average_density = []
    total_gas_mass = []
    average_vinfinity = []
    average_cinfinity = []
    for line in open(output_combined):
        avg_dens = re.search(r'Avg_rho = (.{12})', line)
        avg_cinf = re.search(r'Avg_cinf = (.{12})', line)
        avg_vinf = re.search(r'Avg_vinf = (.{12})', line)
        gas_mass = re.search(r'TotalGasMass within r_k = (.{12})', line)
        if avg_dens:
            average_density.append(avg_dens.group(1))
            total_gas_mass.append(gas_mass.group(1))
            average_cinfinity.append(avg_cinf.group(1))
            average_vinfinity.append(avg_vinf.group(1))

    average_density = np.array(_remove_strings(average_density))
    total_gas_mass = np.array([float(i) for i in total_gas_mass])
    average_vinfinity = np.array([float(i) for i in average_vinfinity])
    average_cinfinity = np.array([float(i) for i in average_cinfinity])


    ##########################################################################################################
    #                                HL radius, Bondi radius and JeansLength
    ##########################################################################################################

    hl_radius = []
    bondi_radius = []
    jeans_length = []
    for line in open(output_combined):
        hl = re.search(r'HLRadius = (.{12})', line)
        bondi = re.search(r'BondiRadius = (.{12})', line)
        jeans = re.search(r'JeansLengthOfRegion = (.{12})', line)
        if hl:
            hl_radius.append(hl.group(1))
            bondi_radius.append(bondi.group(1))
            try:
                jeans_length.append(jeans.group(1))
            except AttributeError:
                continue

    hl_radius = np.array(_remove_strings(hl_radius))
    bondi_radius = np.array(_remove_strings(bondi_radius))
    jeans_length = np.array(_remove_strings(jeans_length))

    ##########################################################################################################
    #                                              BH Mass
    ##########################################################################################################

    mass = []
    for line in open(output_combined):
        bh_mass = re.search(r'cmass = (.{12})', line)
        if bh_mass:
            mass.append(bh_mass.group(1))

    mass = np.array([float(i) for i in mass])

    ##########################################################################################################
    #                                            Temperature
    ##########################################################################################################

    # Assume RegionTemperature is collected here
    region_temperatures = []
    # Temp storage for RegionTemperature values, assuming they might appear more frequently
    temp_region_temps = []

    for line in open(output_combined):
        region_temp = re.search(r'RegionTemperature = (.{12})', line)
        age = re.search(r'Age = (.{12})', line)
        if region_temp:
            # If a RegionTemperature is found, store it temporarily
            temp_region_temps.append(region_temp.group(1))
        elif age:
            # If an 'age' is found, process the most recent RegionTemperature (or however you decide to handle multiple temps)
            if temp_region_temps:
                # For example, taking the last RegionTemperature value before this 'age'
                region_temperatures.append(temp_region_temps[-1])
                temp_region_temps = []  # Reset for the next batch

    ##########################################################################################################
    #                                           List of Arrays
    ##########################################################################################################

    all_data = [ages, accrates, average_density, average_vinfinity, average_cinfinity,
                total_gas_mass, hl_radius, bondi_radius, jeans_length, mass, region_temperatures]

    return all_data


if __name__ == "__main__":

    # to run: python make_csv_from_estd.py

    # Set root directory(ies) of simulation output files
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m08-4dx/2B.m16-4dx-2/"
    #root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/"
    #root_dir = "/Backup01/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSm01-2/"
    #root_dir_2 = "/Backup01/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSm01-2/2S.m01-386+/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSbf4/estd-297+/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSbf4/estd-earlier/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.m16-4dx/"
    #root_dir = "/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/estd-data/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/estd_DD0233+/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-end/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m08-4dx/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-start/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/2B.RSb08-2/"
    #root_dir = "/disk14/sgordon/cirrus-runs-rsync/seed1-bh-only/40msun/replicating-beckmann/1S.RSbf4/estd-data/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/40msun/replicating-beckmann/1S.RSb01/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/estd-early/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/2S.b01-234+/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m04-no-SN/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.b04-no-SN/1S.b04-no-SN-2/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed1-bh-only/seed1-bh-only/40msun/replicating-beckmann/1S.m01-no-SN/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.mf4-no-SN/"
    # root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSmf8-2/"
    # root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSmf4-2/"
    # root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSm01-2/2S.m01-386+/"
    # root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/estd-early/"
    # root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/2S.b01-234+/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.m01-no-SN/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.b01-no-SN/"
    #root_dir_2 = "/Backup01/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/"
    #root_dir = "/Backup01/sgordon/pleiades/seed1-bh-only/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/estd-start/"
    #root_dir = "/Backup01/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb01/"
    #root_dir = "/Backup01/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.m01-no-SN/"
    #root_dir_2 = "/disk14/sgordon/pleiades-11-12-23/seed2-bh-only/40msun/replicating-beckmann-2/2S.m01-no-SN/"
    #root_dir = "/Backup01/sgordon/pleiades/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/estd-early/"
    #root_dir_2 = "/disk14/sgordon/pleiades-11-12-23/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/2S.b01-234+/estds-0.18-1Myr/"
    #root_dir = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/40msun/replicating-beckmann/1S.m01-no-SN/"
    #root_dir = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01/"
    root_dir = "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/"
    #root_dir = "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/"
    #root_dir = "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/"
    #root_dir = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/thermal-fb/1B.th.bf128-eps-0.0001/"
    #root_dir = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/thermal-fb/1B.th.bf128/"
    #root_dir = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/thermal-fb/1B.th.bf128-eps-0.01/"

    # Extract simulation name
    sim = extract_simulation_name(root_dir)

    # Set output file name
    write_to = "data_files/data-" + str(sim + ".csv")

    # Extract data arrays
    all_data = get_data_arrays(root_dir)

    # Assign header columns corresponding to data arrays
    headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
                'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass', 'temperature']
    

    missing_temps_count = 0  # Initialize a counter for missing RegionTemperature values
    temp_index = headerList.index('age') # to identify which column to check for missing values

    # Write to csv file (overwrites existing file)
    with open(write_to, "w+") as f:
        dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
        dw.writeheader()
        writer = csv.writer(f, delimiter=',')
        for values in zip_longest(*all_data):
            # Check if the RegionTemperature value is present
            if values[temp_index] is not None:
                writer.writerow(values)
            else:
                missing_temps_count += 1
                if missing_temps_count > 100:
                    print("More than 40 entries missing for temperature. Current count:", missing_temps_count)
                    break  # Or continue, based on how you want to handle this scenario

    print("Saved data to {}".format(write_to))
    
    # Extract end data from other directory (if it exists) and append to existing csv file
    try:
        all_data_2 = get_data_arrays(root_dir_2)
        with open(write_to, "a+") as f:
            dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
            writer = csv.writer(f, delimiter=',')
            for values in zip_longest(*all_data_2):
                # Check if the RegionTemperature value is present
                if values[temp_index] is not None:
                    writer.writerow(values)
                else:
                    missing_temps_count += 1
                    if missing_temps_count > 100:
                        print("More than 40 entries missing for temperature. Current count:", missing_temps_count)
                        break  # Or continue, based on how you want to handle this scenario
        print("Appended end data to {}".format(write_to))
    except:
        print("No end data found for simulation {}, done.".format(sim))

    # Sort the csv file by age
    print("Sorting the csv file by age...")
    with open(write_to, 'r') as f:
        # Load the CSV data into a pandas DataFrame
        df = pd.read_csv(write_to)

        # Delete rows where 'age' is less than 1000
        df_filtered = df.loc[df['age'] >= 1000]

        # Sorting the DataFrame by the 'age' column in ascending order
        df_sorted = df_filtered.sort_values(by='age', ascending=True)

        # Saving the sorted DataFrame back to the same file, effectively replacing it
        df_sorted.to_csv(write_to, index=False)

        print(f"File '{write_to}' has been sorted and saved successfully.")
