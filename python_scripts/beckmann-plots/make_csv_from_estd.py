import os
import re # complex str searches
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

##########################################################################################################
#                                    Make csv file of BH/gas properties
#
# to run: python make_csv_from_estd_2.py
##########################################################################################################

# turn on to combine multiple estd simulation txt output files
MULTIPLE_ESTDS = 1
SMALL = 0

# user input
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
estds = ['estd_0.out', 'estd_1.out', 'estd_2.out',
         'estd_3.out',
         # 'estd_4.out',
         # 'estd_5.out',
         # 'estd_6.out',
         # 'estd_7.out',
         # 'estd_8.out',
         'estd.out']

# naming output file
seed = int(root_dir[43:44])
if seed == 1:
    index = 82
elif seed == 2:
    index = 84
if SMALL:
    index -= 1

# writing data arrays to this file
write_to = "data_files/data-" + str(root_dir[index:] + ".csv")


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


##########################################################################################################
#                                  Read from simulation output txt file(s)
##########################################################################################################

if MULTIPLE_ESTDS:
    output_combined = str(root_dir[index:])
    path = Path(output_combined)

    files = []
    data_list = []
    for estd in estds:
        f = os.path.join(root_dir, estd)
        files.append(f)
        data_list.append("")

    for i, f in enumerate(files):
        # Reading data from file
        with open(f) as fp:
            data_list[i] = fp.read()

    # Merging files
    # To add the data of file2
    # from next line
    data = ""
    for i in range(len(files)):
        data += "\n"
        data += data_list[i]

    Path(output_combined).touch()
    with open(output_combined, 'w') as fp:
        fp.write(data)
    output = output_combined
else:
    output = os.path.join(root_dir, 'estd.out')


##########################################################################################################
#                                  Accretion Arrays: accrates, ages
##########################################################################################################
accrates = []
ages = []
for line in open(output):
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
for line in open(output):
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
for line in open(output):
    hl = re.search(r'HLRadius = (.{12})', line)
    bondi = re.search(r'BondiRadius = (.{12})', line)
    jeans = re.search(r'JeansLengthOfRegion = (.{12})', line)
    if hl:
        hl_radius.append(hl.group(1))
        bondi_radius.append(bondi.group(1))
        jeans_length.append(jeans.group(1))

hl_radius = np.array(_remove_strings(hl_radius))
bondi_radius = np.array(_remove_strings(bondi_radius))
jeans_length = np.array(_remove_strings(jeans_length))

##########################################################################################################
#                                              BH Mass
##########################################################################################################

mass = []
for line in open(output):
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

with open(write_to, "w+") as f:
    dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
    dw.writeheader()
    writer = csv.writer(f, delimiter=',')
    for values in zip_longest(*all_data):
        writer.writerow(values)

print("saved data to {}".format(write_to))
