import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yt
from yt.data_objects.level_sets.api import * # for clump finding
from helper_functions import *
from helper_functions import _h2_fraction

# Add information to the master clump
def _mass_weighted_jeans_mass(clump, units="Msun"):
    # Assuming this calculates a scalar value correctly
    jeans_mass = clump.quantities.weighted_average_quantity("jeans_mass", ("gas", "mass"))
    # Return as scalar if that's what's expected or adjust to return a tuple/list if needed
    return "Jeans mass: %s.", jeans_mass.to(units)

def _center_of_mass(clump, units="code_length", **kwargs):
    p = clump.quantities.center_of_mass(**kwargs)
    return "Center of mass: %s.", p.to(units)

def _minimum_gas_mass(clump, min_mass):
    return clump["gas", "mass"].sum() >= min_mass
add_validator("minimum_gas_mass", _minimum_gas_mass)

# Find most massive leaf node
# Function to recursively traverse the clump hierarchy and analyze leaf clumps
def process_leaves(master_clump):
    global most_massive_clump
    global max_leaf_mass
    # Initialize variables to keep track of the most massive clump
    most_massive_clump = None
    max_leaf_mass = -1   
    leaf_clumps = master_clump.leaves
    leaf_masses = []
    for clump in leaf_clumps:
        mass = clump.info["cell_mass"][1]
        #print("Clump leaf mass = {}".format(mass))
        leaf_masses.append(mass)
        if mass > max_leaf_mass:
            max_leaf_mass = mass
            most_massive_clump = clump
    print("Number of leaf clumps = {}".format(len(leaf_clumps)))
    print("Most massive leaf clump mass = {:.2f} (clump id {})".format(max_leaf_mass, most_massive_clump.clump_id))
    print("Least massive leaf clump mass = {:.5f} (clump id {})".format(min(leaf_masses), leaf_clumps[np.argmin(leaf_masses)].clump_id))
    return leaf_masses

def plot_histogram(leaf_masses, sim, dd, nbins=50, xlabel="Mass ($M_\odot$)", ylabel="Number of Clumps", logscale=False):
    # Create the histogram using Seaborn
    plt.figure(figsize=(8, 6))
    leaf_masses_values = np.array([mass.value for mass in leaf_masses])
    ax = sns.histplot(leaf_masses_values, bins=nbins, color='pink', kde=False)
    # # Get the counts and bin edges
    counts, bins = np.histogram(leaf_masses, bins=20)
    # Annotate with the counts in each bin
    for count, bin in zip(counts, bins[:-1]):  # Iterate over counts and left edges of bins
        bin_center = bin + (bins[1] - bins[0]) / 2  # Calculate center of bin
        ax.text(bin_center, count, str(int(count)), ha='center', va='bottom')

    plt.title('Distribution of Leaf Clump Masses in Simulation {} at {}'.format(sim, dd))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    #plt.yscale('log') if logscale else None  # Optional: Use log scale for y-axis
    plt.xscale('log') if logscale else None 
    plt.grid(True)  # Optional: Adds a grid for better readability
    plt.savefig(f'plots/leaf_clump_masses_hist_{sim}_{dd}.png') 
    print(f"Saved plots/leaf_clump_masses_hist_{sim}_{dd}.png")
    return 0

def process_master_clump(master_clump, ds, ss_pos, ss_age, sim, dd):
    # Get a list of just the leaf nodes.
    leaf_clumps = master_clump.leaves

    # Process the leaf clumps  
    leaf_masses = process_leaves(master_clump)

    # Projection plot of clump leaves
    for i in range(3):
        p = yt.ProjectionPlot(ds, i, ("gas", "number_density"), center=ss_pos, width=(6, "pc"))
        p.set_cmap('number_density', 'octarine')
        p.set_zlim('number_density', 4e21, 8e25)
        p.annotate_timestamp(corner='upper_right', redshift=True, draw_inset_box=True)
        p.annotate_clumps(leaf_clumps, cmap="Pastel1")
        p.annotate_clumps([most_massive_clump], text_args={"color": "red"})
        most_massive_clump_mass = most_massive_clump.info["cell_mass"][1]
        most_massive_clump_pos = most_massive_clump.info["position"][1].to('pc')
        text_string = f"Most Massive Clump: {most_massive_clump_mass:.2f}\n Total Clump Mass: {sum(leaf_masses):.0f}\n Number of Clumps: {len(leaf_clumps)}\n"
        p.annotate_text((0.02, 0.85), text_string, coord_system='axis', text_args={'color': 'white'})
        p.annotate_text([0.05, 0.05], sim, coord_system="axis", text_args={"color": "black"}, 
                        inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3, 
                                        "edgecolor": "white", "alpha": 0.5}
                        )
        p.annotate_text((0.82, 0.05), "BH Age = {:.2f} Myr".format(ss_age[0]/1e6), coord_system="axis", text_args={"color": "white"})
        p.annotate_marker(master_clump.info["position"][1], coord_system='data', color='white', s=100)
        p.save(f"plots/clump_projection_{sim}_{dd}_{i}.png")
        print(f"Saved plots/clump_projection_{sim}_{dd}_{i}.png")

    # plot histogram
    plot_histogram(leaf_masses, sim, dd, logscale=False)

    return 0
    

if __name__ == "__main__":

    global LOCATE_MASTER_CLUMP
    LOCATE_MASTER_CLUMP = False

    # load the dataset
    fp = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0446/DD0446"
    #fp = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0522/DD0522"
    #fp = "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0517/DD0517"
    ds = yt.load(fp)
    sim = "1B.b01"
    dd = extract_dd_segment(fp)
    print("Simulation: {} and DD: {} at time: {}".format(sim, dd, ds.current_time.to('Myr')))
    ds.add_field(("gas", "h2_fraction"), function=_h2_fraction, units="dimensionless", display_name="H2 Fraction", sampling_type="cell")

    # find sink particle attribute
    ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
    print("Sink particle properties:")
    print("position = {}, mass = {:.2f}, age = {:.2f} Myr".format(ss_pos, ss_mass, ss_age[0]/1e6))

    # plot the density projection in each direction to locate the master clump by eye
    if LOCATE_MASTER_CLUMP:
        for i in range(3):
            p = yt.ProjectionPlot(ds, i, ('gas','number_density'), center=ss_pos, width=(5, 'pc'))
            p.set_zlim("number_density", 9e-3, 6e2)
            p.set_cmap('number_density', 'cividis')
            p.annotate_marker(ss_pos, coord_system='data', color='red', s=100)
            p.annotate_timestamp(corner='upper_left', redshift=True, draw_inset_box=True)
            p.save("plots/density_projection_{}.png".format(i))
            print("Saved plots/density_projection_{}.png".format(i))

    # Make initial master clump (a disk containing the clump) - takes 30 seconds
    clump_pos = ss_pos.to('pc') + [-1.5, -0.5, -1.5]*yt.units.pc
    data_source = ds.sphere(clump_pos, (2, "pc"))
    master_clump = Clump(data_source, ("gas", "density"))

    # Add clump info items
    add_clump_info("mass_weighted_jeans_mass", _mass_weighted_jeans_mass)
    add_clump_info("position",  _center_of_mass)
    master_clump.add_info_item("position")
    master_clump.add_info_item("mass_weighted_jeans_mass")

    print("------------------------------------")
    print("Master clump created")
    print("Master clump position = {}".format(master_clump.info['position']))
    print("Master clump mass = {}".format(master_clump.info['cell_mass']))
    print("Master clump mass weighted Jeans mass = {}".format(master_clump.info['mass_weighted_jeans_mass']))
    print("------------------------------------")

    # Add clump validators (refine clump properties)
    #master_clump.add_validator("min_cells", 20)
    master_clump.add_validator("gravitationally_bound", use_particles=False)
    master_clump.add_validator("minimum_gas_mass", ds.quan(0.5, "Msun"))

    # Run clump finding
    c_min = data_source["gas", "density"].min()*3
    c_max = data_source["gas", "density"].max()
    step = 2
    find_clumps(master_clump, c_min, c_max, step)

    # Save the clump tree as a reloadable dataset
    fn = master_clump.save_as_dataset(fields=[("gas", "density"), ("gas", "h2_fraction"), ("gas", "cell_mass"), ("gas", "jeans_mass"), ("gas", "number_density"), ("gas", "temperature")])
    print("Saved clump tree to {} for simulation {}".format(fn, sim))

    # Process the master clump
    max_leaf_mass = -1
    process_master_clump(master_clump, ds, ss_pos, ss_age, sim, dd)