import yt
import numpy as np

ds0 = yt.load('DD0121/DD0121')
ds1 = yt.load('DD0122/DD0122')

ds0.print_stats()
ds1.print_stats()
