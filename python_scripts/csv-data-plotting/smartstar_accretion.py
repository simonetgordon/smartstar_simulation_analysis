import yt
import os
import sys

# load data
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann-fixed-dx-75%"
ds = yt.load(os.path.join(root_dir, sys.argv[1]))
dd = ds.all_data()

NumAPs = 1
for i in range(NumAPs):
      AccRateTimes = dd["SmartStar", "AccretionRateTime"][i]
      timeindex = dd["SmartStar", "TimeIndex"][i].d
      accrate = dd["SmartStar", "AccretionRate"][i][1:int(timeindex+1)].in_units("Msun/yr")
      ages = dd["SmartStar", "AccretionRateTime"][i][1:int(timeindex+1)].in_units("kyr") - \
             dd["SmartStar", "creation_time"][i].in_units("kyr")

print(AccRateTimes)
print(accrate)