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
            if age:
                ages.append(float(age.group(1)))
            if avg_dens:
                average_density.append(float(avg_dens.group(1)))
            if avg_cinf:
                average_cinfinity.append(float(avg_cinf.group(1)))
            if avg_vinf:
                average_vinfinity.append(float(avg_vinf.group(1)))
            if gas_mass:
                total_gas_mass.append(float(gas_mass.group(1)))
            if hl:
                hl_radius.append(float(hl.group(1)))
            if bondi:
                bondi_radius.append(float(bondi.group(1)))
            if jeans:
                jeans_length.append(float(jeans.group(1)))
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

def combine_data_from_estd_files(directory):
    estd_files = get_sorted_estd_files(directory)

    # Initialize lists for data columns
    accrates, ages, average_density, total_gas_mass, average_vinfinity = [], [], [], [], []
    average_cinfinity, hl_radius, bondi_radius, jeans_length, mass = [], [], [], [], []

    for file in estd_files:
        with open(os.path.join(directory, file), "r") as infile:
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
                if age:
                    ages.append(float(age.group(1)))
                if avg_dens:
                    average_density.append(float(avg_dens.group(1)))
                if avg_cinf:
                    average_cinfinity.append(float(avg_cinf.group(1)))
                if avg_vinf:
                    average_vinfinity.append(float(avg_vinf.group(1)))
                if gas_mass:
                    total_gas_mass.append(float(gas_mass.group(1)))
                if hl:
                    hl_radius.append(float(hl.group(1)))
                if bondi:
                    bondi_radius.append(float(bondi.group(1)))
                if jeans:
                    jeans_length.append(float(jeans.group(1)))
                if bh_mass:
                    mass.append(float(bh_mass.group(1)))

    return ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass

def main(root_dir):
    sim = extract_simulation_name(root_dir)
    print("Creating combined .csv for simulation", sim)

    # Combine data from the root directory
    ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass, hl_radius, bondi_radius, jeans_length, mass = combine_data_from_estd_files(root_dir)

    # Define header columns
    headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
                  'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass']

    # Write data to CSV file
    write_to = "data_files/data-" + str(sim) + ".csv"
    write_csv_header(write_to, headerList)

    # Append data from the root directory to the CSV file
    append_data_to_csv(write_to, [ages, accrates, average_density, average_vinfinity, average_cinfinity,
                                  total_gas_mass, hl_radius, bondi_radius, jeans_length, mass])
    
    # Combine data from the second directory (if it exists)
    second_directory = ""
    second_directory = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-end/"  # Replace with the correct path to the second directory

    # Check if the second directory exists and append its data to the CSV file
    if os.path.isdir(second_directory):
        print("Creating combined .csv for secondary part of simulation", sim)
        second_ages, second_accrates, second_average_density, second_average_vinfinity, second_average_cinfinity, second_total_gas_mass, second_hl_radius, second_bondi_radius, second_jeans_length, second_mass = combine_data_from_estd_files(second_directory)

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
