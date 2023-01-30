import os
import re # complex str searches
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

# reading data from this directory
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann-fixed-dx-75%"
output_combined = 'output-fixed-dx-75%.out'

# writing data arrays to this file
write_to = "data-fixed-dx-75%.csv"

##########################################################################################################
#                                  Read from simulation output txt files
##########################################################################################################

path = Path(output_combined)
if path.is_file():

    file1 = os.path.join(root_dir, 'estd_0.out')
    file2 = os.path.join(root_dir, 'estd.out')
    data = data2 = ""

    # Reading data from file1
    with open(file1) as fp:
        data = fp.read()

    # Reading data from file2
    with open(file2) as fp:
        data2 = fp.read()

    # Merging 2 files
    # To add the data of file2
    # from next line
    data += "\n"
    data += data2

    Path(output_combined).touch()
    with open(output_combined, 'w') as fp:
        fp.write(data)
    output = output_combined
else:
    output = output_combined


##########################################################################################################
#                                  Accretion Arrays: accrates, accrate_times
##########################################################################################################
accrates = []
accrate_dtimes = []
accrate_times = []
for line in open(output):
    accrate = re.search(r'accrate = (.{8})', line)
    accrate_dtime = re.search(r'deltatime = (.{7})', line)
    if accrate:
        accrates.append(accrate.group(1))
        accrate_dtimes.append(accrate_dtime.group(1))

accrates = np.array([float(i) for i in accrates])
accrate_dtimes = np.array([float(i) for i in accrate_dtimes])

for j in range(len(accrate_dtimes)):
    accrate_times.append(accrate_dtimes[j] + sum(accrate_dtimes[:j]))

accrate_times = np.array(accrate_times)


##########################################################################################################
#                        Average Velocities, Densities + Temperature Arrays
##########################################################################################################

def _avg_densities(avg_densities):
    for index, i in enumerate(avg_densities):
        try:
            float(i)
        except ValueError:
            i = i.replace("cm", '')
            i = i.replace("^-", '')
            try:
                i = float(i)
            except ValueError:
                i = i[:-1]
                i = float(i)
        avg_densities[index] = float(i)
    return avg_densities

avg_densities = []
avg_temperatures = []
avg_vinfinities = []
avg_cinfinities = []
for line in open(output):
    avg_dens = re.search(r'Avg_Density = (.{11})', line)
    avg_temp = re.search(r'AverageTemp = (.{12})', line)
    avg_cinf = re.search(r'Average cInfinity = (.{12})', line)
    avg_vinf = re.search(r'Average vInfinity = (.{12})', line)
    if avg_dens:
        avg_densities.append(avg_dens.group(1))
        avg_temperatures.append(avg_temp.group(1))
        avg_cinfinities.append(avg_cinf.group(1))
        avg_vinfinities.append(avg_vinf.group(1))

avg_densities = np.array(_avg_densities(avg_densities))
avg_temperatures = np.array([float(i) for i in avg_temperatures])
avg_vinfinities = np.array([float(i) for i in avg_vinfinities])
avg_cinfinities = np.array([float(i) for i in avg_cinfinities])

avg_times = np.linspace(accrate_dtimes[0], sum(accrate_dtimes), num=len(avg_cinfinities))


##########################################################################################################
#                                              HL radius
##########################################################################################################

hl_radii = []
for line in open(output):
    hl_radius = re.search(r'to BondiHoyle radius = (.{12})', line)
    if hl_radius:
        hl_radii.append(hl_radius.group(1))

hl_radii = np.array([float(i) for i in hl_radii])


##########################################################################################################
#                                              BH Mass
##########################################################################################################

bh_masses = []
for line in open(output):
    bh_mass = re.search(r'NewMass = (.{12})', line)
    if bh_mass:
        bh_masses.append(bh_mass.group(1))
bh_masses = np.array([float(i) for i in bh_masses])


##########################################################################################################
#                                           Write to csv file
##########################################################################################################

# assign header columns
headerList = ['accrate times', 'accrate', 'parameter times', 'average density', 'average vinfinity',
              'average cinfinity', 'average temperature', 'HL radius', 'BH mass']

# open CSV file and assign header
all_data = [accrate_times, accrates, avg_times, avg_vinfinities, avg_cinfinities, avg_temperatures,
            hl_radii, bh_masses]

with open(write_to, "w+") as f:
    dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
    dw.writeheader()
    writer = csv.writer(f, delimiter=',')
    for values in zip_longest(*all_data):
        writer.writerow(values)

print("saved data to {}".format(write_to))
