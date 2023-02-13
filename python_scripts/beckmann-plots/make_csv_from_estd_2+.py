import os
import re # complex str searches
from pathlib import Path
import numpy as np
import csv
from itertools import zip_longest

# turn on to combine multiple estd simulation txt output files
MULTIPLE_ESTDS = 1

# turn on if estd file has both AvgValues_MassWeighted and AvgValues printed out
MASS_WEIGHTED = 1

# reading data from this directory
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm16"

# writing data arrays to this file
write_to = "data_files/data-1B.RSm16.csv"

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
    output_combined = 'output-combined-temp-3.out'
    path = Path(output_combined)
    if not path.is_file():
        file1 = os.path.join(root_dir, 'estd_0.out')
        file2 = os.path.join(root_dir, 'estd_1.out')
        file3 = os.path.join(root_dir, 'estd.out')
        #file4 = os.path.join(root_dir, 'estd.out')
        data = data2 = data3 = data4 = ""

        # Reading data from file1
        with open(file1) as fp:
            data = fp.read()

        # Reading data from file2
        with open(file2) as fp:
            data2 = fp.read()

        # Reading data from file2
        with open(file3) as fp:
            data3 = fp.read()

        # Reading data from file2
        # with open(file4) as fp:
        #     data4 = fp.read()

        # Merging 2 files
        # To add the data of file2
        # from next line
        data += "\n"
        data += data2
        data += "\n"
        data += data3
        # data += "\n"
        # data += data4

        Path(output_combined).touch()
        with open(output_combined, 'w') as fp:
            fp.write(data)
        output = output_combined
    else:
        output = output_combined
else:
    output = os.path.join(root_dir, 'estd.out')


##########################################################################################################
#                                  Accretion Arrays: accrates, accrate_times
##########################################################################################################
accrates = []
accrate_dtimes = []
accrate_times = []
for line in open(output):
    accrate = re.search(r'accrate = (.{8})', line)
    accrate_dtime = re.search(r'deltatime = (.{7})', line)
    if accrate and accrate_dtime:
        accrates.append(accrate.group(1))
        accrate_dtimes.append(accrate_dtime.group(1))

accrates = np.array([float(i) for i in accrates])
accrate_dtimes = np.array(_remove_strings(accrate_dtimes))

for j in range(len(accrate_dtimes)):
    accrate_times.append(accrate_dtimes[j] + sum(accrate_dtimes[:j]))

accrate_times = np.array(accrate_times)

#accrate_times = np.linspace(10000, 1760000, len(accrates))


##########################################################################################################
#                        Average Velocities, Densities + Temperature Arrays
##########################################################################################################

average_density = []
average_temperature = []
average_vinfinity = []
average_cinfinity = []
for line in open(output):
    avg_dens = re.search(r'Avg_Density = (.{11})', line)
    avg_temp = re.search(r'AverageTemp = (.{12})', line)
    avg_cinf = re.search(r'Average cInfinity = (.{12})', line)
    avg_vinf = re.search(r'Average vInfinity = (.{12})', line)
    if MASS_WEIGHTED:
        mass_weighted = re.search(r'AvgValues_MassWeighted:', line)
    else:
        mass_weighted = 1
    if avg_dens and mass_weighted:
        average_density.append(avg_dens.group(1))
        average_temperature.append(avg_temp.group(1))
        average_cinfinity.append(avg_cinf.group(1))
        average_vinfinity.append(avg_vinf.group(1))

average_density = np.array(_remove_strings(average_density))
average_temperature = np.array([float(i) for i in average_temperature])
average_vinfinity = np.array([float(i) for i in average_vinfinity])
average_cinfinity = np.array([float(i) for i in average_cinfinity])

average_times = np.linspace(accrate_dtimes[0], sum(accrate_dtimes), num=len(average_cinfinity))
print("sum(accrate_dtimes) = ", sum(accrate_dtimes))
#average_times = np.linspace(accrate_times[0], accrate_times[-1], num=len(average_cinfinity))


##########################################################################################################
#                                              HL radius
##########################################################################################################

hl_radius = []
for line in open(output):
    radius = re.search(r'to BondiHoyle radius = (.{12})', line)
    if radius:
        hl_radius.append(radius.group(1))

hl_radius = np.array([float(i) for i in hl_radius])
hl_times = np.linspace(accrate_dtimes[0], sum(accrate_dtimes), num=len(hl_radius))

##########################################################################################################
#                                              BH Mass
##########################################################################################################

mass = []
for line in open(output):
    bh_mass = re.search(r'NewMass = (.{12})', line)
    if bh_mass:
        mass.append(bh_mass.group(1))
mass = np.array([float(i) for i in mass])
mass_times = np.linspace(accrate_dtimes[0], sum(accrate_dtimes), num=len(mass))


##########################################################################################################
#                                           Write to csv file
##########################################################################################################

# assign header columns
headerList = ['accrate times', 'accrate', 'average times', 'average density', 'average vinfinity',
              'average cinfinity', 'average temperature', 'HL radius', 'BH mass', 'hl_times', 'mass_times']

# open CSV file and assign header
all_data = [accrate_times, accrates, average_times, average_density, average_vinfinity, average_cinfinity,
            average_temperature, hl_radius, mass, hl_times, mass_times]

with open(write_to, "w+") as f:
    dw = csv.DictWriter(f, delimiter=',', fieldnames=headerList)
    dw.writeheader()
    writer = csv.writer(f, delimiter=',')
    for values in zip_longest(*all_data):
        writer.writerow(values)

print("saved data to {}".format(write_to))