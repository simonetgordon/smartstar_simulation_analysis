import os
import re # complex str searches
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

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
    #                                           List of Arrays
    ##########################################################################################################

    all_data = [ages, accrates, average_density, average_vinfinity, average_cinfinity,
                total_gas_mass, hl_radius, bondi_radius, jeans_length, mass]

    return all_data


if __name__ == "__main__":

    # Set root directory(ies) of simulation output files
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.m16-4dx/"
    root_dir = "/ceph/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm04/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSm08/"
    #root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb16/"
    #root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-end/"

    # Extract simulation name
    sim = extract_simulation_name(root_dir)

    # Set output file name
    write_to = "data_files/data-" + str(sim + ".csv")

    # Extract data arrays
    all_data = get_data_arrays(root_dir)

    # Assign header columns corresponding to data arrays
    headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
                'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass']
    
    # Write to csv file (overwrites existing file)
    with open(write_to, "w+") as f:
        dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
        dw.writeheader()
        writer = csv.writer(f, delimiter=',')
        for values in zip_longest(*all_data):
            writer.writerow(values)
    print("Saved data to {}".format(write_to))
    
    # Extract end data from other directory (if it exists) and append to existing csv file
    try:
        all_data_2 = get_data_arrays(root_dir_2)
        with open(write_to, "a+") as f:
            dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
            writer = csv.writer(f, delimiter=',')
            for values in zip_longest(*all_data_2):
                writer.writerow(values)
        print("appended end data to {}".format(write_to))
    except:
        print("No end data found for simulation {}, done.".format(sim))
