import logging
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yt
from yt.data_objects.level_sets.api import *
import yt.units  # for clump finding
from helper_functions import *
from helper_functions import _h2_fraction
import yt.extensions.p2p.clumps
from yt.extensions.p2p.misc import iterate_center_of_mass, reunit
import os
import csv

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

# Conversion factor from g/cm³ to hydrogen nuclei/cm³
mh = 1.6735e-24  # Mass of a hydrogen atom in grams

# Ensure the directory for the CSV file exists
csv_dir = 'csv-data-plotting'
csv_file = os.path.join(csv_dir, 'clump_data.csv')
os.makedirs(csv_dir, exist_ok=True)

# Define the header for the CSV file if it doesn't exist
if not os.path.exists(csv_file):
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['simulation', 'age_myr', 'datadump', 'boundedness', 'future_bound', 'min_density', 'no_clumps', 'max_clump_mass', 'total_clump_mass', 'bh_mass'])

def extract_specific_part(filepath):
    pattern = r'/(.*?)/DD0.+/DD0.+$'
    match = re.search(pattern, filepath)
    if match:
        specific_part = match.group(1).split('/')[-1]
        return specific_part
    else:
        return None

def extract_dd_segment(filepath):
    # Assuming the datadump segment is in the form 'DD####'
    match = re.search(r'DD(\d+)', filepath)
    if match:
        return match.group(1)
    else:
        return None
    
# Add information to the master clump
def _mass_weighted_jeans_mass(clump, units="Msun"):
    jeans_mass = clump.quantities.weighted_average_quantity("jeans_mass", ("gas", "mass"))
    return "Jeans mass: %s.", jeans_mass.to(units)

def _center_of_mass(clump, units="code_length", **kwargs):
    p = clump.quantities.center_of_mass(**kwargs)
    return "Center of mass: %s.", p.to(units)

def _minimum_gas_mass(clump, min_mass):
    return clump["gas", "mass"].sum() >= min_mass
add_validator("minimum_gas_mass", _minimum_gas_mass)

def process_leaves(master_clump):
    global most_massive_clump
    global max_leaf_mass
    most_massive_clump = None
    max_leaf_mass = -1
    leaf_clumps = master_clump.leaves
    leaf_masses = []
    for clump in leaf_clumps:
        mass = clump.info["cell_mass"][1]
        leaf_masses.append(mass)
        if mass > max_leaf_mass:
            max_leaf_mass = mass
            most_massive_clump = clump
    return leaf_masses

def process_master_clump(master_clump, ds, ss_pos, ss_age, ss_mass, sim, dd, c_min_hn):
    leaf_clumps = master_clump.leaves
    leaf_masses = process_leaves(master_clump)
    total_clump_mass = sum(leaf_masses).value
    no_clumps = len(leaf_masses)
    max_clump_mass = max(leaf_masses).value if leaf_masses else 0
    future_bound = " True"

    # Write data to CSV
    with open(csv_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            sim,
            ss_age[0]/1e6,
            dd,
            "True",
            future_bound,
            c_min_hn,
            no_clumps,
            max_clump_mass,
            total_clump_mass,
            ss_mass.value
        ])

    # Generate plots
    plot_histogram(leaf_masses, sim, dd, ss_age, c_min_hn)
    plot_projection(ds, ss_pos, ss_age, ss_mass, leaf_clumps, most_massive_clump, sim, dd, c_min_hn)

def plot_histogram(leaf_masses, sim, dd, ss_age, c_min, nbins=50, xlabel="Mass ($M_\\odot$)", ylabel="Number of Clumps", logscale=False):
    plt.figure(figsize=(8, 6))
    leaf_masses_values = np.array([mass.value for mass in leaf_masses])
    ax = sns.histplot(leaf_masses_values, bins=nbins, color='pink', kde=False)

    min_mass = np.min(leaf_masses_values)
    max_mass = np.max(leaf_masses_values)
    num_leaves = len(leaf_masses_values)

    counts, bins = np.histogram(leaf_masses, bins=20)
    for count, bin in zip(counts, bins[:-1]):
        bin_center = bin + (bins[1] - bins[0]) / 2
        ax.text(bin_center, count, str(int(count)), ha='center', va='bottom')

    text_str = f"Min Mass: {min_mass:.2e} $M_\\odot$\nMax Mass: {max_mass:.2e} $M_\\odot$\nNumber of Leaves: {num_leaves}\n Min Density: {c_min:.2e} cm$^{-3}$"
    plt.text(0.95, 0.95, text_str, transform=ax.transAxes, horizontalalignment='right', verticalalignment='top', fontsize=10, bbox=dict(facecolor='white', alpha=0.5))

    plt.title(f'Distribution of Leaf Clump Masses in Simulation {sim} at {ss_age[0]/1e6:.2f} Myr')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale('log') if logscale else None 
    plt.grid(True)
    plot_name = f"plots/clump_masses_hist_{sim}_{dd}_future_bound.png"
    plt.savefig(plot_name)
    plt.close()
    print(f"Saved histogram to {plot_name}")

def plot_projection(ds, ss_pos, ss_age, ss_mass, leaf_clumps, most_massive_clump, sim, dd, c_min):
    for i in range(3):
        p = yt.ProjectionPlot(ds, i, ("gas", "number_density"), center=ss_pos, width=(6, "pc"))
        p.set_cmap('number_density', 'octarine')
        p.set_zlim('number_density', 4e21, 8e25)
        p.annotate_timestamp(corner='upper_right', redshift=True, draw_inset_box=True)
        p.annotate_clumps(leaf_clumps, cmap="Pastel1")
        p.annotate_clumps([most_massive_clump], text_args={"color": "black"})
        most_massive_clump_mass = most_massive_clump.info["cell_mass"][1]
        text_string = f"Most Massive Clump: {most_massive_clump_mass:.2f}\n Total Clump Mass: {sum([c.info['cell_mass'][1].value for c in leaf_clumps]):.0f}\n Number of Clumps: {len(leaf_clumps)}\n"
        p.annotate_text((0.02, 0.85), text_string, coord_system='axis', text_args={'color': 'white'})
        p.annotate_text([0.05, 0.05], sim, coord_system="axis", text_args={"color": "black"}, 
                        inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3, 
                                        "edgecolor": "white", "alpha": 0.5})
        p.annotate_text((0.68, 0.02), f"BH Age = {ss_age[0]/1e6:.2f} Myr", coord_system="axis", text_args={"color": "white"})
        p.annotate_text((0.68, 0.06), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis", text_args={"color": "white"})
        plot_name = f"plots/clump_projection_{sim}_{dd}_{i}_{c_min:.2e}_future_bound.png"
        p.save(plot_name)
        print(f"Saved projection plot to {plot_name}")

def process_file(fp, c_min, c_max):
    try:
        logging.info(f"Started processing {fp}")
        # Convert c_min and c_max from g/cm³ to hydrogen nuclei/cm³
        c_min_hn = c_min / mh
        c_max_hn = c_max / mh
        
        ds = yt.load(fp)
        sim = extract_specific_part(fp)
        dd = extract_dd_segment(fp)
        logging.info(f"Simulation: {sim} and DD: {dd} at time: {ds.current_time.to('Myr')}")
        
        ds.add_field(("gas", "thermal_energy"), function=specific_thermal_energy, units="erg/g", sampling_type="cell")
        ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
        logging.info(f"BH properties: position = {ss_pos}, mass = {ss_mass:.2f}, age = {ss_age[0]/1e6:.2f} Myr")
        
        field = ("gas", "density")
        data_source = ds.sphere(ss_pos, (4, "pc"))
        master_clump = Clump(data_source, field)
        step = 2
        logging.info(f"Density range in hydrogen nuclei/cm³: {c_min_hn:.2e} to {c_max_hn:.2e}")
        
        output_dir = "clumps/"
        ensure_dir(output_dir)
        master_clump.add_validator("future_bound",
            use_thermal_energy=True,
            truncate=True,
            include_cooling=True
            )
        master_clump.add_validator("minimum_gas_mass", min_mass=2*yt.units.msun)
        add_clump_info("mass_weighted_jeans_mass", _mass_weighted_jeans_mass)
        add_clump_info("position",  _center_of_mass)
        master_clump.add_info_item("position")
        master_clump.add_info_item("mass_weighted_jeans_mass")
        master_clump.add_info_item("center_of_mass")
        master_clump.add_info_item("min_number_density")
        master_clump.add_info_item("max_number_density")
        master_clump.add_info_item("jeans_mass")
        
        find_clumps(master_clump, c_min, c_max, step)
        clump_ds_name = f"_{sim}" + output_dir
        fn = master_clump.save_as_dataset(filename=clump_ds_name, fields=["density", ("gas", "H2_fraction"), ("gas", "cell_mass"), ("gas", "jeans_mass"), ("gas", "number_density"), ("gas", "temperature")])
        logging.info(f"Saved clump tree to {fn} for simulation {sim}")
        
        process_master_clump(master_clump, ds, ss_pos, ss_age, ss_mass, sim, dd, c_min_hn)
        logging.info("Master clump processed")
        logging.info("------------------------------------")
        logging.info(f"Clump finding parameters:")
        logging.info(f"Density range: {c_min_hn:.2e} to {c_max_hn:.2e}")
        logging.info(f"Future Bound Clumps? True")
        logging.info(f"Step: {step}")
        logging.info(f"Master clump position = {master_clump.info['position']}")
        logging.info(f"Master clump mass = {master_clump.info['cell_mass']}")
        logging.info(f"Master clump mass weighted Jeans mass = {master_clump.info['mass_weighted_jeans_mass']}")
        logging.info("------------------------------------")
        return fn
    except Exception as e:
        logging.error(f"Error processing {fp}: {e}")
        raise

def parallel_clump_finding(filepaths, c_min, c_max, num_workers=5):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_file, fp, c_min, c_max) for fp in filepaths]
        results = []
        for future in futures:
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                results.append(e)
                logging.error(f"Future processing error: {e}")
    return results

if __name__ == "__main__":
    # Define the variables directly in the script
    filepaths = [
        # 1B.b01
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0445/DD0445", # 31.70
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0446/DD0446", # 31.80
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0447/DD0447", # 31.90
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0448/DD0448", # 32.00
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0449/DD0449", # 32.10
        "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0450/DD0450", # 32.21.
        # 1B.resim.th.b01-3-eta-0.1
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0445/DD0445", # 31.70
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0573/DD0573", # 31.77
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0603/DD0603", # 31.80
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0623/DD0623", # 31.82
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0445/DD0445", # 31.70
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0574/DD0574", # 31.77
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0604/DD0604", # 31.80

        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0450/DD0450", # 31.70
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0565/DD0565",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0595/DD0595",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0612/DD0612",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0655/DD0655",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695", # 31.90
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD1141/DD1141", # 31.95
        # 1B.resim.th.b01-3-eta-0.01
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0445/DD0445", # 31.70
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0515/DD0515", # 31.77
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0545/DD0545", # 31.80
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0565/DD0565", # 31.82
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0605/DD0605", # 31.86
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0646/DD0646", # 31.90
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0696/DD0696"  # 31.95
         #1B.resim.th.b01-3-eta-0.001
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0445/DD0445", # 31.70
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0514/DD0514", # 31.77
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0545/DD0545", # 31.80
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0561/DD0561", # 31.82
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2745/DD2745"  # 31.86
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2783/DD2783"  # 31.90
        # 1B.resim.th.b01-3-eta-0.0001
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0522/DD0522", # 31.77
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0552/DD0552", # 31.80
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0571/DD0571", # 31.82
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0607/DD0607", # 31.86
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0652/DD0652", # 31.90
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0701/DD0701" 
    ]
    c_min = 1e-21  # g/cm³
    c_max = 1e-13
    num_workers = len(filepaths) + 1

    # Run the parallel clump finding
    results = parallel_clump_finding(filepaths, c_min, c_max, num_workers)
    for result in results:
        print(result)
    logging.info("All clump finding processes completed")