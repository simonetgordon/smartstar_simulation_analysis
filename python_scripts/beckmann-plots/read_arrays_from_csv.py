import csv
import sys
import pandas as pd

filename = sys.argv[-1]

data = pd.read_csv(filename, sep=',')
columns = data.columns.values

class Black_Hole:
    def __init__(self, accrate_times, accrates, average_times, average_density, average_vinfinity,
                 average_cinfinity, average_temperature, hl_radius, mass):
        self.accrate_times = accrate_times
        self.accrates = accrates
        self.average_times = average_times
        self.average_density = average_density
        self.average_vinfinity = average_vinfinity
        self.average_cinfinity = average_cinfinity
        self.average_temperature = average_temperature
        self.hl_radius = hl_radius
        self.mass = mass

BHL = Black_Hole(data[columns[0]].values, data[columns[1]].values, data[columns[2]].values, data[columns[3]].values,
           data[columns[4]].values, data[columns[5]].values,
           data[columns[6]].values, data[columns[7]].values, data[columns[8]].values)

print(BHL.accrates)