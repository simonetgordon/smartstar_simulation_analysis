import yt
import sys
import os
import numpy as np
from smartstar_find import ss_properties
from scipy import stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from derived_fields import add_fields_ds
from find_disc_attributes import _make_disk_L, ds, ss_pos

# set by user
root_dir = "/home/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2"
input = sys.argv[1]
ds = yt.load(os.path.join(root_dir, sys.argv[1]))
add_fields_ds(ds)

# naming plot
seed = int(root_dir[43:44])
print(seed)
if seed == 1:
    index = 82
elif seed == 2:
    index = 84

# make disk data container
disk, L = _make_disk_L(ds, ss_pos, 10, 1)
fn = disk.save_as_dataset(fields=[("gas", "density"), ("gas", "temperature"), ("gas", "H_nuclei_density"),
                                  ("index", "radius"), ("gas", "mass"), ("index", "cylindrical_z"),
                                  ("gas", "angular_frequency"), ("gas", "dynamical_time")])
disk_ds = yt.load(fn)
ad = disk_ds.all_data()


# make density projection and use frb to extract 2D data array
dx = 1.229791e-02  # pc
disc_r_pc = 10
w = disc_r_pc*2*yt.units.pc # pc -  0.00076892 in code units
frb_resolution = int(w/dx)
proj = ds.proj(("gas", "H_nuclei_density"), "x", center=ss_pos, ds=ds, data_source=disk) # z has highest sigma value
disc_frb = proj.to_frb((disc_r_pc*2, "pc"), resolution=(int(frb_resolution)), center=ss_pos)
disc_frb_ds = disc_frb.export_dataset(fields=[("gas", "H_nuclei_density"), ("index", "radius"), ("index", "x")])

# sigma
y = disc_frb_ds.r[("gas", "H_nuclei_density")]

# radius
xleft = np.abs(disc_frb_ds.domain_left_edge[0] - ss_pos[0]).to('pc')
xright = np.abs(disc_frb_ds.domain_right_edge[0] - ss_pos[0]).to('pc')
yleft = np.abs(disc_frb_ds.domain_left_edge[1] - ss_pos[1]).to('pc')
yright = np.abs(disc_frb_ds.domain_right_edge[1] - ss_pos[1]).to('pc')
print(xleft, yleft)
print(xright, yright)
xleft = xright = yleft = yright = 10
xs = np.linspace(-xleft, xright, y.shape[0]) # xwidth/frb_resolution
ys = np.linspace(-yleft, yright, y.shape[0])
r = np.sqrt(xs**2 + ys**2)

# plot sigma v radius
bins = np.logspace(-2, 1, 64)
r_bins = np.histogram_bin_edges(r, bins=bins)
sigma, radius = np.histogram(r, weights=y, bins=bins)
fig = plt.figure()
plt.plot(radius[:-1], sigma)
plt.ylabel('$\Sigma \, (cm^{-2})$')
plt.xlabel('$R \, (pc)$')
plt.yscale('log')
plt.xscale('log')
plot_name = 'disc_sigma_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.png'
fig.savefig('plots/' + plot_name, dpi=100)
print("created plots/" + str(plot_name))


# profile = yt.create_profile(
#     data_source=disk,
#     bin_fields=[("index", "radius")],
#     fields=[("gas", "H_nuclei_density"), ("gas", "dynamical_time"),
#             ("gas", "density")],
#     n_bins=128,
#     units=dict(radius="pc"),
#     logs=dict(radius=False),
#     #weight_field=None,
#     accumulation=False
# )
#
# # plot omega vs radius to find d(omega^2)/dr
# x = profile.x[profile.used]
# y = profile[("gas", "dynamical_time")][profile.used].to('yr')
# plt.plot(x, y)
# print("time: ", x)
# print("radius: ", y)
# slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(y))
# print("slope = ", slope)
# pc_cm = 3.086e+18
# newX = np.logspace(17, 20, 128, base=10)
# linX = np.linspace(0.01*pc_cm, 10*pc_cm, 128)
#
# # Let's fit an exponential function.
# # This looks like a line on a lof-log plot.
# def myExpFunc(x, a, b):
#     return a * np.power(x, b)
#
# popt, pcov = curve_fit(myExpFunc, x, y)
# plt.plot(newX, myExpFunc(newX, 10**intercept, slope), 'green', label='fitted line')
# #plt.plot(newX, myExpFunc(newX, *popt), 'r-', label="({0:.2E}*x**{1:.3f})".format(*popt))
# plt.legend()
# print("Exponential Fit: y = (a*(x**b))")
# print("\ta = popt[0] = {0}\n\tb = slope = {1}".format(*popt))

#plt.ylabel(r'$\rho \, (g \, cm^{-3})$')
# plt.ylabel(r'$t_{dyn} \, (yr)$')
# plt.xlabel('Radius (pc)')
# plt.yscale('log')
# plt.xscale('log')
# plot_name = 'disc_r_tdyn_' + str(root_dir[index:]) + '_' + str(input)[7:] + '.png'
# plt.savefig('plots/' + plot_name, dpi=100)
# print("created plots/" + str(plot_name))