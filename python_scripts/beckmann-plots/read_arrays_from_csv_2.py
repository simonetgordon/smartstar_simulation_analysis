import sys
import pandas as pd

##########################################################################################################
#                               Make BHL.attributes from csv columns
#
# to run: python read_arrays_from_csv.py [csv_filename]
# to call attributes: BHL.accrates
# imported to plot_variables.py
##########################################################################################################


# Black Hole class
class BlackHole:
    def __init__(self, ages, accrates, average_density, average_vinfinity, average_cinfinity, total_gas_mass,
                 hl_radius, bondi_radius, jeans_length, mass):
        self.ages = ages
        self.accrates = accrates
        self.average_density = average_density
        self.average_vinfinity = average_vinfinity
        self.average_cinfinity = average_cinfinity
        self.total_gas_mass = total_gas_mass
        self.hl_radius = hl_radius
        self.mass = mass
        self.bondi_radius = bondi_radius
        self.jeans_length = jeans_length

    def info(self):
        print("------------------------------------------------------------------------------------------------------")
        print("ages (yrs post-formation): {}, \naccrates (Msun/yr): {}, "
              "\naverage_density (cm^-3): {}, "
              "\naverage_vinfinity (km/s): {}, \naverage_cinfinity (km/s): {}, "
              "\ntotal_gas_mass (K): {}, \nhl_radius (pc): {}, \nbondi_radius (pc): {}"
              "\nmass (Msun): {}, \njeans length (pc): {}".format(
            self.ages, self.accrates, self.average_density,
            self.average_vinfinity, self.average_cinfinity, self.total_gas_mass,
            self.hl_radius, self.bondi_radius, self.mass, self.jeans_length))
        print("------------------------------------------------------------------------------------------------------")
        return


# BHL object creator function (to be called in plot_variables)
def make_bhl_object(sys_arg):
    filename = sys_arg
    data = pd.read_csv(filename, sep=',')
    columns = data.columns.values

    # make BHL class object
    output_object = BlackHole(data[columns[0]].values, data[columns[1]].values, data[columns[2]].values,
                     data[columns[3]].values, data[columns[4]].values, data[columns[5]].values,
                     data[columns[6]].values, data[columns[7]].values, data[columns[8]].values,
                     data[columns[9]].values)
    output_object.info()
    return output_object


def make_rolling_bhl_object(sys_arg):
    filename = sys_arg
    data = pd.read_csv(filename, sep=',')
    columns = data.columns.values
    # 1000 * 100 = 100,000 year rolling average?
    #window = int(len(data[columns[7]].values)*0.01)
    window = 1000

    # make BHL class object
    output_object = BlackHole(data[columns[0]].values, data[columns[1]].values,
                              data[columns[2]].values, data[columns[3]].rolling(window=window).mean(),
                              data[columns[4]].rolling(window=window).mean(), data[columns[5]].rolling(window=window).mean(),
                              data[columns[6]].rolling(window=window).mean(), data[columns[7]].rolling(window=window).mean(),
                              data[columns[8]].values, data[columns[9]].rolling(window=window).mean())
    output_object.info()
    return output_object


def bhl_object_list(rolling=0):
    bhl_objects = []
    for i in range(1, (len(sys.argv)-1)):
        if rolling:
            bhl = make_rolling_bhl_object(sys.argv[i])
        else:
            bhl = make_bhl_object(sys.argv[i])
        bhl_objects.append(bhl)
    return bhl_objects


def bhl_object_labels():
    bhl_labels = []
    for i in range(1, (len(sys.argv)-1)):
        bhl = str(sys.argv[i])
        bhl = bhl.replace(".csv", '')
        bhl = bhl.replace("data-", '')
        bhl = bhl.replace("data_files/", '')
        bhl = bhl.replace("s1-", '')
        bhl_labels.append(bhl)
    return bhl_labels


bhl_object_list = bhl_object_list(rolling=0)
bhl_object_labels = bhl_object_labels()
