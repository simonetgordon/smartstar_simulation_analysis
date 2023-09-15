import os
import re
import numpy as np
import csv
from itertools import zip_longest
from pathlib import Path

def _remove_strings(lst):
    # Use list comprehension for string replacements
    replacements = {
        "inf pc \t Jea": "0",
        "inf pc \t Bon": "0",
        "cm^-3": '',
        "cm^-": '',
        "cm^": '',
        "cm": '',
        "c": '',
        ",": '',
        "y": ''
    }
    lst = [replacements[i] if i in replacements else i for i in lst]
    return [float(i) for i in lst]

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


def read_data_from_file(filename):
    # Initialize lists for data
    accrates, ages, average_density, total_gas_mass, average_vinfinity = [], [], [], [], []
    average_cinfinity, hl_radius, bondi_radius, jeans_length, mass = [], [], [], [], []

    with open(filename) as infile:
        for line in infile:
            # Use regex groups to directly extract values
            accrate = re.search(r'accrate = (\S+)', line)
            age = re.search(r'Age = (\S+)', line)
            avg_dens = re.search(r'Avg_rho = (\S+)', line)
            avg_cinf = re.search(r'Avg_cinf = (\S+)', line)
            avg_vinf = re.search(r'Avg_vinf = (\S+)', line)
            gas_mass = re.search(r'TotalGasMass within r_k = (\S+)', line)
            hl = re.search(r'HLRadius = (\S+)', line)
            bondi = re.search(r'BondiRadius = (\S+)', line)
            jeans = re.search(r'JeansLengthOfRegion = (\S+)', line)
            bh_mass = re.search(r'cmass = (\S+)', line)

            if accrate:
                accrates.append(float(accrate.group(1)))
                ages.append(float(age.group(1)))
            if avg_dens:
                average_density.append(float(avg_dens.group(1)))
                average_cinfinity.append(float(avg_cinf.group(1)))
                average_vinfinity.append(float(avg_vinf.group(1)))
                total_gas_mass.append(float(gas_mass.group(1)))
            if hl:
                hl_radius.append(float(hl.group(1)))
                bondi_radius.append(float(bondi.group(1)))
                try:
                    jeans_length.append(jeans.group(1))
                except AttributeError:
                    continue
            if bh_mass:
                mass.append(float(bh_mass.group(1)))

    return ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass

def write_csv_header(filename, headerList):
    if not os.path.isfile(filename):
        with open(filename, "w", newline='') as f:
            dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
            dw.writeheader()

def append_data_to_csv(filename, data):
    with open(filename, 'a', newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for values in zip_longest(*data):
            writer.writerow(values)

def get_sorted_estd_files(directory):
    estd_files = [f for f in os.listdir(directory) if f.endswith(".out") and f.startswith("estd")]
    estd_files.sort(key=lambda x: int(re.search(r'estd_(\d+)', x).group(1)) if re.search(r'estd_(\d+)', x) else float('inf'))
    return estd_files

def get_combined_estds_file(directory, sim):
    # Get list of estd*.out files in directory
    estds = [f for f in os.listdir(directory) if f.endswith(".out") and f.startswith("estd")]

    # Collect all estd*.out files in the directory
    out_files = []
    misc_files = []
    for file in os.listdir(directory):
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

    output_combined = 'estd_files/' + sim

    with open(output_combined, "w") as outfile:
        # Append contents of each .out file to output file
        for out_file in estds:# writing data arrays to this file
            with open(directory + out_file, "r") as infile:
                outfile.write(infile.read())

    return output_combined

def get_data_from_estd_file(output_combined):
     # Initialize lists for data columns
    accrates, ages, average_density, total_gas_mass, average_vinfinity = [], [], [], [], []
    average_cinfinity, hl_radius, bondi_radius, jeans_length, mass = [], [], [], [], []

    for line in open(output_combined):
        accrate = re.search(r'accrate = (.{8})', line)
        age = re.search(r'Age = (.{12})', line)
        if accrate:
            accrates.append(accrate.group(1))
            ages.append(age.group(1))

    accrates = np.array([float(i) for i in accrates])
    ages = np.array([float(i) for i in ages])

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

    for line in open(output_combined):
        bh_mass = re.search(r'cmass = (.{12})', line)
        if bh_mass:
            mass.append(bh_mass.group(1))

    mass = np.array([float(i) for i in mass])

    return accrates, ages, average_density, total_gas_mass, average_vinfinity, average_cinfinity, hl_radius, bondi_radius, jeans_length, mass

def combine_data_from_estd_files(directory):
    estd_files = get_sorted_estd_files(directory)
    #estd_files = collect_and_sort_files(directory)

    # Initialize lists for data columns
    accrates, ages, average_density, total_gas_mass, average_vinfinity = [], [], [], [], []
    average_cinfinity, hl_radius, bondi_radius, jeans_length, mass = [], [], [], [], []

    for file in estd_files:
        with open(os.path.join(directory, file), "r") as infile:
            for line in infile:
                # Use regex groups to directly extract values
                accrate = re.search(r'accrate =(.{8})', line)
                age = re.search(r'Age = (.{12})', line)
                avg_dens = re.search(r'Avg_rho = (\S+)', line)
                avg_cinf = re.search(r'Avg_cinf = (\S+)', line)
                avg_vinf = re.search(r'Avg_vinf = (\S+)', line)
                gas_mass = re.search(r'TotalGasMass within r_k = (\S+)', line)
                hl = re.search(r'HLRadius = (\S+)', line)
                bondi = re.search(r'BondiRadius = (\S+)', line)
                jeans = re.search(r'JeansLengthOfRegion = (\S+)', line)
                bh_mass = re.search(r'cmass = (.{12})', line)

                if accrate:
                    accrates.append(float(accrate.group(1)))
                    ages.append(float(age.group(1)))
                if avg_dens:
                    average_density.append(float(avg_dens.group(1)))
                    average_cinfinity.append(float(avg_cinf.group(1)))
                    average_vinfinity.append(float(avg_vinf.group(1)))
                    total_gas_mass.append(float(gas_mass.group(1)))
                if hl:
                    hl_radius.append(float(hl.group(1)))
                    bondi_radius.append(float(bondi.group(1)))
                    try:
                        jeans_length.append(jeans.group(1))
                    except AttributeError:
                        continue
                if bh_mass:
                    mass.append(float(bh_mass.group(1)))

    return ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass

def main(root_dir):
    sim = extract_simulation_name(root_dir)
    print("Creating combined .csv for simulation", sim)
    output_combined = get_combined_estds_file(root_dir, sim)

    # Combine data from the root directory
    #ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass = combine_data_from_estd_files(root_dir)
    ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass = get_data_from_estd_file(output_combined)
    

    # Define header columns
    headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
                  'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass']

    # Write data to CSV file
    write_to = "data_files/data-" + str(sim) + ".csv"
    # Check if file exists
    if os.path.isfile(write_to):
        # If file exists, remove it
        os.remove(write_to)
        print(f"The current {write_to} file has been removed and will be replaced.")
    else:
        print(f"No action taken. File {write_to} will be created.")
    # Write to clear file
    write_csv_header(write_to, headerList)

    # Append data from the root directory to the CSV file
    append_data_to_csv(write_to, [ages, accrates, average_density, average_vinfinity, average_cinfinity,
                                  total_gas_mass, hl_radius, bondi_radius, jeans_length, mass])
    
    # Combine data from the second directory (if it exists)
    second_directory = ""
    second_directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-end/"  # Replace with the correct path to the second directory

    # Check if the second directory exists and append its data to the CSV file
    if os.path.isdir(second_directory):
        sim = extract_simulation_name(root_dir) + "_2"
        print("Creating combined .csv for secondary part of simulation", sim)
        output_combined = get_combined_estds_file(root_dir, sim)
        second_ages, second_accrates, second_average_density, second_average_vinfinity, second_average_cinfinity, second_total_gas_mass, second_hl_radius, second_bondi_radius, second_jeans_length, second_mass = get_data_from_estd_file(output_combined)

        # Append data from the second directory to the CSV file
        append_data_to_csv(write_to, [second_ages, second_accrates, second_average_density, second_average_vinfinity,
                                      second_average_cinfinity, second_total_gas_mass, second_hl_radius, second_bondi_radius,
                                      second_jeans_length, second_mass])

    print("Saved data to", write_to)


if __name__ == "__main__":
    # User input: Replace with your desired root directory
    #root_directory = "/disk14/sgordon/cirrus-runs-rsync/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSb01/"
    root_directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-start/"
    main(root_directory)
