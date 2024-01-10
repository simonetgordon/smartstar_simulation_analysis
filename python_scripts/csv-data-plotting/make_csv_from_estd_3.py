import os
import re
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

def _remove_strings(lst):
    for index, i in enumerate(lst):
        try:
            float(i)
        except ValueError:
            i = replace_inf_and_units(i)
            float(i)
        lst[index] = float(i)
    return lst

def replace_inf_and_units(i):
    i = i.replace("inf pc \t Jea", "0")
    i = i.replace("inf pc \t Bon", "0")
    i = i.replace("cm^-3", '')
    i = i.replace("cm^-", '')
    i = i.replace("cm^", '')
    i = i.replace("cm", '')
    i = i.replace("c", '')
    i = i.replace(",", '')
    i = i.replace("y", '')
    return i

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

def collect_and_sort_outfiles(root_dir):
    out_files = []
    misc_files = []
    for file in os.listdir(root_dir):
        if file.endswith(".out"):
            if file.startswith("estd_"):
                out_files.append(file)
            else:
                misc_files.append(file)
    estds = sorted(out_files, key=lambda x: int(x.split("_")[1].split(".")[0]) 
                   if len(x.split("_")[1].split(".")[0]) <= 2 else float('inf'))
    estds.extend(misc_files)
    return estds

def combine_outfiles(root_dir, estds, sim):
    output_combined = 'estd_files/' + sim
    with open(output_combined, "w") as outfile:
        for out_file in estds:
            with open(root_dir + out_file, "r") as infile:
                outfile.write(infile.read())
    print("Created combined .csv for simulation ", sim)
    return output_combined

def extract_from_line(line, pattern):
    match = re.search(pattern, line)
    return match.group(1) if match else None

def extract_from_file(filename, pattern):
    result = []
    with open(filename) as f:
        for line in f:
            extraction = extract_from_line(line, pattern)
            if extraction is not None:
                result.append(extraction)
    return result

def write_to_csv(sim, *data):
    headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
                  'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass']
    all_data = data
    write_to = "data_files/data-" + str(sim + ".csv")
    with open(write_to, "a+") as f:
        dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
        dw.writeheader()
        writer = csv.writer(f, delimiter=',')
        for values in zip_longest(*all_data):
            writer.writerow(values)
    print("saved data to {}".format(write_to))

def main(root_dir):
    estds = collect_and_sort_outfiles(root_dir)
    sim = extract_simulation_name(root_dir)
    print("Creating combined .csv for simulation ", sim)
    output_combined = combine_outfiles(root_dir, estds, sim)

    accrates = extract_from_file(output_combined, r'accrate = (.{8})')
    ages = extract_from_file(output_combined, r'Age = (.{12})')

    accrates = np.array([float(i) for i in accrates])
    ages = np.array([float(i) for i in ages])

    average_density = extract_from_file(output_combined, r'Avg_rho = (.{12})')
    average_vinfinity = extract_from_file(output_combined, r'Avg_vinf = (.{12})')
    average_cinfinity = extract_from_file(output_combined, r'Avg_cinf = (.{12})')
    total_gas_mass = extract_from_file(output_combined, r'TotalGasMass within r_k = (.{12})')

    average_density = np.array(_remove_strings(average_density))
    total_gas_mass = np.array([float(i) for i in total_gas_mass])
    average_vinfinity = np.array([float(i) for i in average_vinfinity])
    average_cinfinity = np.array([float(i) for i in average_cinfinity])

    hl_radius = extract_from_file(output_combined, r'HLRadius = (.{12})')
    bondi_radius = extract_from_file(output_combined, r'BondiRadius = (.{12})')
    jeans_length = extract_from_file(output_combined, r'JeansLengthOfRegion = (.{12})')

    hl_radius = np.array(_remove_strings(hl_radius))
    bondi_radius = np.array(_remove_strings(bondi_radius))
    jeans_length = np.array(_remove_strings(jeans_length))

    mass = extract_from_file(output_combined, r'cmass = (.{12})')
    mass = np.array([float(i) for i in mass])

    write_to_csv(sim, ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass,
                 hl_radius, bondi_radius, jeans_length, mass)

if __name__ == "__main__":
    root_dir = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-start/"
    root_dir_2 = "/ceph/cephfs/sgordon/pleiades/seed2-bh-only/270msun/replicating-beckmann-2/2B.RSb08/estd-data-end/"
    sim = extract_simulation_name(root_dir)
    write_to = "data_files/data-" + str(sim + ".csv")
    if os.path.exists(write_to):
        os.remove(write_to)

    # Extract data and write to csv
    main(root_dir)
    try:
        main(root_dir_2)
    except:
        print("No end data found for simulation ", sim)


