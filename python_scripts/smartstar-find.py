"""
Find and print SmartStar properties
"""
import sys
import os
import yt
import ytree

# set by use
root_dir = "~/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSm04"
ds = yt.load(os.path.join(root_dir, sys.argv[1]))

# make sphere
width = (0.2, 'unitary')
sp = ds.sphere('max', width)
ad = ds.all_data()

# find ss properties
ss_creation = ad['SmartStar', 'creation_time'].to('yr')
ss_pos = ad['SmartStar', 'particle_position'].to('unitary')[0]
ss_mass = ad['SmartStar', 'particle_mass_'].to('Msun')
ss_class = ad[('SmartStar', 'ParticleClass_')] # 0 = POPIII, 1 = SMS, 2 = BH

# find ss age
time_array = ds.current_time.to('yr')
creation = ss_creation.d # strip units off variables
time = time_array.d
ss_age = time - creation

if __name__ == "__main__":
    print("-----------------------------------------")
    print("current particle position [x,y,z]:")
    print(ss_pos.d)
    print("-----------------------------------------")
    print("time of particle creation =", ss_creation)
    print("-----------------------------------------")
    print("current time =", ds.current_time.to('yr'))
    print("-----------------------------------------")
    print("age of particle in yr =", ss_age)
    print("-----------------------------------------")
    print("particle class_: ", ss_class)
    print("-----------------------------------------")
    print("particle mass_: ", ss_mass)
    print("-----------------------------------------")

    # Load merger tree of dataset (up to DD0118 in gas run)
    a = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/out_0.list')
    # Load halo
    ds_halos = yt.load('/home/sgordon/disk14/pop3/gas+dm-L3/rockstar_halos/halos_36.0.bin')
    ds_data = yt.load('/home/sgordon/disk14/pop3/dm-only-L0/DD0036/DD0036')
    # Load my_tree and find radius
    a1 = ytree.load('/home/sgordon/disk14/pop3/gas+dm-L3/tree_810/tree_810.h5')
    r_halo = a1[0]["virial_radius"].to('pc')
    r = ds.quan(r_halo.d, "pc")  # virial radius
    print("virial radius (pc) ", r)
