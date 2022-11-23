# Print SmartStar info #
import sys
import yt

ds = yt.load(sys.argv[1])
# Star forms in DD0122 so restart from DD0121
center0 = [0.485819995403, 0.490299999714, 0.506219983101]
width = (0.2, 'unitary')

sp = ds.sphere(center0, width)
ad = ds.all_data()

star_creation = ad['SmartStar', 'creation_time'].to('yr')
star_pos1 = ad['SmartStar', 'particle_position_x'].to('unitary')
star_pos2 = ad['SmartStar', 'particle_position_y'].to('unitary')
star_pos3 = ad['SmartStar', 'particle_position_z'].to('unitary')
star_mass = ad['SmartStar', 'particle_mass'].to('Msun')
star_class = ad[('SmartStar', 'ParticleClass')] # 0 = POPIII, 1 = SMS, 2 = BH 
#ad[('SmartStar', 'AccretionRate')] # throws error: 'tuple cannot be interpreted as int'

time_array = ds.current_time.to('yr')
creation = star_creation.d # strip units off variables
time = time_array.d
star_age = time - creation

# Print to terminal
print("-----------------------------------------")
print("current particle position [x,y,z]:")
print(star_pos1.d, star_pos2.d, star_pos3.d)
print("-----------------------------------------")
print("time of particle creation =", star_creation)
print("-----------------------------------------")
print("current time =", ds.current_time.to('yr'))
print("-----------------------------------------")
print("age of particle in yr =", star_age)
print("-----------------------------------------")
print("particle class: ", star_class)
print("-----------------------------------------")
print("particle mass: ", star_mass)
print("-----------------------------------------")
print("accretion rate: ", ad[('SmartStar', 'AccretionRate')].in_units('Msun/yr'))
print("-----------------------------------------")
