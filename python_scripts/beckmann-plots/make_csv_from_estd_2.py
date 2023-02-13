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

# turn on if estd file has both AvgValues_MassWeighted and AvgValues printed out
MASS_WEIGHTED = 1

# reading data from this directory
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm04-2"

# writing data arrays to this file
write_to = "data_files/data-" + str(root_dir[82:] + ".csv")


def _remove_strings(lst):
    for index, i in enumerate(lst):
        try:
            float(i)
        except ValueError:
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
    output_combined = str(root_dir[82:])
    path = Path(output_combined)
    if not path.is_file():
        file1 = os.path.join(root_dir, 'estd.out')
        # file2 = os.path.join(root_dir, 'estd_1.out')
        # file3 = os.path.join(root_dir, 'estd_2.out')
        # file4 = os.path.join(root_dir, 'estd_3.out')
        # file5 = os.path.join(root_dir, 'estd.out')
        data = data2 = data3 = data4 = data5 = ""

        # Reading data from file1
        with open(file1) as fp:
            data = fp.read()

        # # Reading data from file2
        # with open(file2) as fp:
        #     data2 = fp.read()
        #
        # # # Reading data from file2
        # with open(file3) as fp:
        #     data3 = fp.read()
        # #
        # # # Reading data from file2
        # with open(file4) as fp:
        #     data4 = fp.read()
        #
        # # # Reading data from file5
        # with open(file5) as fp:
        #     data5 = fp.read()

        # Merging 2 files
        # To add the data of file2
        # from next line
        data += "\n"
        # data += data2
        # data += "\n"
        # data += data3
        # data += "\n"
        # data += data4
        # data += "\n"
        # data += data5

        Path(output_combined).touch()
        with open(output_combined, 'w') as fp:
            fp.write(data)
        output = output_combined
    else:
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

hl_radius = np.array([float(i) for i in hl_radius])
bondi_radius = np.array([float(i) for i in bondi_radius])
jeans_length = np.array([float(i) for i in jeans_length])

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
