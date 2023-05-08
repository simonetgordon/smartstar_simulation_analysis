import os
import re # complex str searches
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

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


def extract_penultimate_directory(filepath):
    # split the filepath into directory and file components
    directory, filename = os.path.split(filepath)
    # split the directory component into its path elements
    path_elements = directory.split(os.path.sep)
    # return the penultimate element, or None if not found
    return path_elements[-2] if len(path_elements) > 1 else None

##########################################################################################################
#                                    Make csv file of BH/gas properties
#
# to run: python make_csv_from_estd.py
##########################################################################################################

# user input
root_dir = "/cephfs/sgordon/cirrus-runs-rsync/seed2-bh-only/40msun/replicating-beckmann-2/2S.RSmf64/estd-data/"

# Get list of .out files in directory
estds = [f for f in os.listdir(root_dir) if f.endswith(".out")]

# Collect all .out files in the root directory
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
sim = extract_penultimate_directory(root_dir)
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

output = output_combined
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
#                                           Write to csv file
##########################################################################################################

# assign header columns
headerList = ['age', 'accrate', 'average density', 'average vinfinity', 'average cinfinity',
              'total gas mass', 'HL radius', 'Bondi radius', 'Jeans length', 'BH mass']

# open CSV file and assign header
all_data = [ages, accrates, average_density, average_vinfinity, average_cinfinity,
            total_gas_mass, hl_radius, bondi_radius, jeans_length, mass]

write_to = "data_files/data-" + str(sim + ".csv")
with open(write_to, "w+") as f:
    dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
    dw.writeheader()
    writer = csv.writer(f, delimiter=',')
    for values in zip_longest(*all_data):
        writer.writerow(values)

print("saved data to {}".format(write_to))
