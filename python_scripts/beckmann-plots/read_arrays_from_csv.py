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
    def __init__(self, accrate_times, accrates, average_times, average_density, average_vinfinity,
                 average_cinfinity, average_temperature, hl_radius, mass, hl_times, mass_times):
        self.accrate_times = accrate_times
        self.accrates = accrates
        self.average_times = average_times
        self.average_density = average_density
        self.average_vinfinity = average_vinfinity
        self.average_cinfinity = average_cinfinity
        self.average_temperature = average_temperature
        self.hl_radius = hl_radius
        self.mass = mass
        self.hl_times = hl_times
        self.mass_times = mass_times

    def info(self):
        print("------------------------------------------------------------------------------------------------------")
        print("accrate_times (yrs post-formation): {}, \naccrates (Msun/yr): {}, \naverage_times (yrs "
              "post-formation): {}, \naverage_density (cm^-3): {}, \naverage_vinfinity (km/s): {}, "
              "\naverage_cinfinity (km/s): {}, \naverage_temperature (K): {}, \nhl_radius (pc): {},"
              "\nmass (Msun): {}".format(
            self.accrate_times,self.accrates, self.average_times, self.average_density,
            self.average_vinfinity, self.average_cinfinity, self.average_temperature,
            self.hl_radius, self.mass))
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
                     data[columns[9]].values, data[columns[10]].values)
    output_object.info()
    return output_object


def make_rolling_bhl_object(sys_arg):
    filename = sys_arg
    data = pd.read_csv(filename, sep=',')
    columns = data.columns.values
    # 1000 * 100 = 100,000 year rolling average?
    window = int(len(data[columns[7]].values)*0.01)

    # make BHL class object
    output_object = BlackHole(data[columns[0]].rolling(window=int(window*0.2)).mean(), data[columns[1]].rolling(window=int(window*0.2)).mean(),
                              data[columns[2]].rolling(window=window).mean(), data[columns[3]].rolling(window=window).mean(),
                              data[columns[4]].rolling(window=window).mean(), data[columns[5]].rolling(window=window).mean(),
                              data[columns[6]].rolling(window=window).mean(), data[columns[7]].rolling(window=window).mean(),
                              data[columns[8]].rolling(window=window).mean(), data[columns[9]].rolling(window=window).mean(),
                              data[columns[10]].rolling(window=window).mean())
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


bhl_object_list = bhl_object_list(rolling=1)
bhl_object_labels = bhl_object_labels()
