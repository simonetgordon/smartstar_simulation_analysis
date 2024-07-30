import logging
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import yt
import yt.extensions.p2p.clumps
import yt.units
import os
import csv
from yt.data_objects.level_sets.api import *
from concurrent.futures import ProcessPoolExecutor
from helper_functions import *


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s')

# Conversion factor from g/cm³ to hydrogen nuclei/cm³
mh = 1.6735e-24  # Mass of a hydrogen atom in grams

# Ensure the directory for the CSV files exists
csv_dir = 'csv-data-plotting'
csv_file = os.path.join(csv_dir, 'clump_data.csv')
distance_csv_file = os.path.join(csv_dir, 'clump_distances.csv')
os.makedirs(csv_dir, exist_ok=True)

# Define the headers for the CSV files if they don't exist
if not os.path.exists(csv_file):
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['simulation', 'age_myr', 'datadump', 'boundedness', 'future_bound', 'min_density', 'no_clumps', 'max_clump_mass', 'total_clump_mass', 'bh_mass', 'average_clump_distance'])

if not os.path.exists(distance_csv_file):
    with open(distance_csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['simulation', 'age_myr', 'datadump', 'min_density', 'clump_id', 'clump_mass', 'clump_distance'])

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

def calculate_distance_to_bh(clump, ss_pos):
    clump_position = _center_of_mass(clump, units="pc")[1]
    bh_position = ss_pos.to('pc')
    distance = np.sqrt(np.sum((clump_position - bh_position) ** 2))
    return distance

def process_leaves(master_clump, bh_position):
    global most_massive_clump
    global max_leaf_mass
    most_massive_clump = None
    max_leaf_mass = -1
    leaf_clumps = master_clump.leaves
    leaf_masses = []
    distances_to_bh = []
    clump_ids = []
    for clump in leaf_clumps:
        mass = clump.info["cell_mass"][1]
        leaf_masses.append(mass)
        distance_to_bh = calculate_distance_to_bh(clump, bh_position)
        distances_to_bh.append(float(distance_to_bh))
        clump_ids.append(clump.clump_id)
        if mass > max_leaf_mass:
            max_leaf_mass = mass
            most_massive_clump = clump
    return leaf_masses, distances_to_bh, clump_ids

def process_master_clump(master_clump, ds, ss_pos, ss_age, ss_mass, sim, dd, c_min_hn):
    leaf_masses, distances_to_bh, clump_ids = process_leaves(master_clump, ss_pos)
    total_clump_mass = sum(leaf_masses).value
    no_clumps = len(leaf_masses)
    max_clump_mass = max(leaf_masses).value if leaf_masses else 0
    future_bound = "True"
    average_clump_distance = sum(distances_to_bh) / no_clumps if no_clumps > 0 else 0

    # Write data to clump_data.csv
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
            ss_mass.value,
            average_clump_distance
        ])

    # Write individual clump distances to clump_distances.csv
    with open(distance_csv_file, 'a', newline='') as f:
        writer = csv.writer(f)
        for clump_id, mass, distance in zip(clump_ids, leaf_masses, distances_to_bh):
            writer.writerow([
                sim,
                ss_age[0]/1e6,
                dd,
                c_min_hn,
                clump_id,  # clump_id
                mass.value,
                distance
            ])

    # Generate plots
    plot_histogram(leaf_masses, sim, dd, ss_age, c_min_hn)
    plot_projection(ds, ss_pos, ss_age, ss_mass, master_clump.leaves, most_massive_clump, clump_ids, sim, dd, c_min_hn)

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

def plot_projection(ds, ss_pos, ss_age, ss_mass, leaf_clumps, most_massive_clump, clump_ids, sim, dd, c_min):
    for i in range(3):
        p = yt.ProjectionPlot(ds, i, ("gas", "number_density"), center=ss_pos, width=(3, "pc"))
        p.set_cmap('number_density', 'octarine')
        p.set_zlim('number_density', 4e21, 8e25)
        p.annotate_timestamp(corner='upper_right', redshift=True, draw_inset_box=True)
        p.annotate_clumps(leaf_clumps, cmap="Pastel1")
        p.annotate_clumps([most_massive_clump], text_args={"color": "black"})
        most_massive_clump_mass = most_massive_clump.info["cell_mass"][1]
        text_string = f"Most Massive Clump: {most_massive_clump_mass:.2f}\n Total Clump Mass: {sum([c.info['cell_mass'][1].value for c in leaf_clumps]):.0f}\n Number of Clumps: {len(leaf_clumps)}\n"
        p.annotate_text((0.02, 0.85), text_string, coord_system='axis', text_args={'color': 'white'})

        for clump, clump_id in zip(leaf_clumps, clump_ids):
            p.annotate_text(clump.quantities.center_of_mass(), f"{clump_id}", coord_system='data', text_args={"color": "white", "fontsize": 8, "weight": "bold"})

        p.annotate_text([0.05, 0.05], sim, coord_system="axis", text_args={"color": "black"}, 
                        inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3, 
                                        "edgecolor": "white", "alpha": 0.5})
        p.annotate_text((0.68, 0.02), f"BH Age = {ss_age[0]/1e6:.2f} Myr", coord_system="axis", text_args={"color": "white"})
        p.annotate_text((0.68, 0.06), r"BH Mass: {:.2f} $\rm M_\odot$".format(ss_mass.d), coord_system="axis", text_args={"color": "white"})
        plot_name = f"plots/clump_projection_{sim}_{dd}_{i}_mindens={c_min:.2e}_future_bound.png"
        p.save(plot_name)
        print(f"Saved projection plot to {plot_name}")

def process_file(fp, c_min, c_max, sim_name=None):
    try:
        logging.info(f"Started processing {fp}")
        # Convert c_min and c_max from g/cm³ to hydrogen nuclei/cm³
        c_min_hn = c_min / mh
        c_max_hn = c_max / mh
        
        ds = yt.load(fp)
        sim = extract_specific_part(fp) if sim_name is None else sim_name
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
        clump_ds_name = f"_{sim}_" + output_dir
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

def parallel_clump_finding(filepaths, c_min, c_max, num_workers=5, sim_name=None):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_file, fp, c_min, c_max, sim_name=sim_name) for fp in filepaths]
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
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0150/DD0150", # 2.201
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0151/DD0151", # 2.301
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0152/DD0152", # 2.401
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0153/DD0153", # 2.501
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0154/DD0154", # 2.601
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0155/DD0155", # 2.701
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0156/DD0156", # 2.801
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0157/DD0157", # 2.901
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0158/DD0158", # 3.001
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0159/DD0159", # 3.101
    "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0160/DD0160", # 3.201
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0445/DD0445", # 31.70
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0446/DD0446", # 31.80
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0447/DD0447", # 31.90
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0448/DD0448", # 32.00
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0449/DD0449", # 32.10
    # "/Backup01/sgordon/disk14/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0450/DD0450", # 32.21
    # 1B.b04
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0150/DD0150", # 2.201
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0151/DD0151", # 2.301
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0152/DD0152", # 2.401
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0153/DD0153", # 2.501
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0154/DD0154", # 2.601
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0155/DD0155", # 2.701
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0156/DD0156", # 2.801
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0157/DD0157", # 2.901
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0158/DD0158", # 3.001
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0159/DD0159", # 3.101
    "/Backup01/sgordon/disk14/cirrus-runs-rsync/seed1-bh-only/270msun/replicating-beckmann/1B.RSb04/DD0160/DD0160", # 3.201
    # # 1B.resim.th.b01-3-eta-0.1-fb-radius
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0445/DD0445", # 31.70
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0573/DD0573", # 31.77
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0603/DD0603", # 31.80
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx/DD0623/DD0623", # 31.82
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0445/DD0445", # 31.70
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0574/DD0574", # 31.77
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=20dx/DD0604/DD0604", # 31.80
    # # 1B.resim.th.b01-3-eta-0.1
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0450/DD0450", # 31.70
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0565/DD0565",
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermFal-fb/1B.resim.th.b01-3-eta-0.1/DD0595/DD0595",
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0612/DD0612",
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0655/DD0655",
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695", # 31.90
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD1141/DD1141", # 31.95
    # 1B.resim.th.b01-3-eta-0.01
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0445/DD0445", # 31.70
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0515/DD0515", # 31.77
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0545/DD0545", # 31.80
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0565/DD0565", # 31.82
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0605/DD0605", # 31.86
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0646/DD0646", # 31.90
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0696/DD0696",  # 31.95
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0719/DD0719",  # 31.975
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/DD0745/DD0745",  # 32.00

    # #1B.resim.th.b01-3-eta-0.001
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0445/DD0445", # 31.70
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0469/DD0469", # 31.725
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0495/DD0495", # 31.75
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0514/DD0514", # 31.77
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0545/DD0545", # 31.80
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD0561/DD0561", # 31.82
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2745/DD2745",  # 31.86
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2783/DD2783",  # 31.90
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2808/DD2808",  # 31.925
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2833/DD2833",  # 31.95
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2858/DD2858",  # 31.975
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001/DD2883/DD2883",  # 32.00 - not processed
    # # 1B.resim.th.b01-3-eta-0.0001
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0445/DD0445", # 31.70
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0469/DD0469", # 31.725
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0511/DD0511", # 31.75
    #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0522/DD0522", # 31.77
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0552/DD0552", # 31.80
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0571/DD0571", # 31.82
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0607/DD0607", # 31.86
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0652/DD0652", # 31.90
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0676/DD0676", # 31.925
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001/DD0704/DD0704",  # 31.95
    # 1B.resim.th.b04 - end event 2.80
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0153/DD0153", # 2.500
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0177/DD0177", # 2.525
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0203/DD0203", # 2.550
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0227/DD0227", # 2.575
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0253/DD0253", # 2.600
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0277/DD0277", # 2.625
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0303/DD0303", # 2.651
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0327/DD0327", # 2.675
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0353/DD0353", # 2.700
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0377/DD0377", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0403/DD0403", # 2.751
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0427/DD0427", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0453/DD0453", # 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0477/DD0477", # 2.825
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04/DD0503/DD0503", # 2.851
    # 1B.resim.th.b04-eps-0.001
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001/DD0153/DD0153", # 2.500
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001/DD0177/DD0177", # 2.525
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001/DD0202/DD0202", # 2.550
    # # 1B.resim.th.b04-eps-0.001-2
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0155/DD0155", # 2.701
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0179/DD0179", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0204/DD0204", # 2.750
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0229/DD0229", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0254/DD0254", # 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eps-0.001-2/DD0264/DD0264", # 2.810 - temp
    # # 1B.resim.th.b04-eta-0.1
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0153/DD0153", # 2.501
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0177/DD0177", # 2.525
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0202/DD0202", # 2.550
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0227/DD0227", # 2.575
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0252/DD0252", # 2.600
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0277/DD0277", # 2.625
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0302/DD0302", # 2.650
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0327/DD0327", # 2.675
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0352/DD0352", # 2.700
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0377/DD0377", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0402/DD0402", # 2.750
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0408/DD0408", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0411/DD0411", # 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0413/DD0413", # 2.825
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0416/DD0416", # 2.850
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0418/DD0418", # 2.875
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0421/DD0421", # 2.900
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0423/DD0423", # 2.925
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0426/DD0426", # 2.950
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/DD0428/DD0428", # 2.975
    # 1B.resim.th.b04-eta-0.1-fb-r-7dx
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0153/DD0153", # 2.500
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0177/DD0177", # 2.525
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-7dx/DD0187/DD0187", # 2.535 - temp
    # 1B.resim.th.b04-eta-0.1-fb-r-10dx
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0153/DD0153", # 2.500
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0177/DD0177", # 2.525
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0202/DD0202", # 2.550
    "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1/fb-r-10dx/DD0227/DD0227", # 2.575
    # # 1B.resim.th.b04-eta-0.1-2
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0153/DD0153", # 2.700
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0179/DD0179", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0204/DD0204", # 2.750
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0227/DD0227", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0254/DD0254", # 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0279/DD0279", # 2.825
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0304/DD0304", # 2.850
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0329/DD0329", # 2.875
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b04-eta-0.1-2/DD0354/DD0354", # 2.900
    # 1B.resim.th.b01-eps-0.001-2
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0155/DD0155", # 2.701
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0179/DD0179", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0204/DD0204", # 2.750
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0229/DD0229", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0254/DD0254",# 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0279/DD0279", # 2.825
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0304/DD0304", # 2.850
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0329/DD0329", # 2.875
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0354/DD0354", # 2.900
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0379/DD0379", # 2.925
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0404/DD0404", # 2.950
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0429/DD0429", # 2.975
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eps-0.001-2/DD0440/DD0440", # temp - 2.990
    # # 1B.resim.th.b01-eta-0.1-2
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0155/DD0155", # 2.701
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0179/DD0179", # 2.725
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0204/DD0204", # 2.750
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0229/DD0229", # 2.775
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0254/DD0254", # 2.800
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0279/DD0279", # 2.825
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0304/DD0304", # 2.850
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0329/DD0329", # 2.875
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0354/DD0354", # 2.900
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0379/DD0379", # 2.925
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0404/DD0404", # 2.950
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0429/DD0429", # 2.975
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0454/DD0454", # 3.000
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0479/DD0479", # 3.025
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0504/DD0504", # 3.050
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0529/DD0529", # 3.075
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0554/DD0554", # 3.100
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0579/DD0579", # 3.125
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0604/DD0604", # 3.150
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0629/DD0629", # 3.175
    # "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-eta-0.1-2/DD0654/DD0654", # 3.200
    # 2B.resim.th.b01
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01/DD0561/DD0561", # 6.530
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01/DD0721/DD0721", # 6.546
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01/DD0821/DD0821", # 6.621
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01/DD0866/DD0866", # 6.666
    # # 2B.resim.th.b01-eps-0.001
    #"/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD0263/DD0263", # 6.500
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD0363/DD0363", # 6.510
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD0463/DD0463", # 6.520
    #"/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD0563/DD0563", # 6.530
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD0962/DD0962", # 6.570
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD1002/DD1002", # 6.570
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD1062/DD1062", # 6.580
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD1014/DD1014", # 6.675
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001/DD2330/DD2330", # 6.707
    # # 2B.resim.th.b01-eps-0.01
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0263/DD0263", # 6.500
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0334/DD0334", # 6.525
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0360/DD0360", # 6.550
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0384/DD0384", # 6.575
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0410/DD0410", # 6.600
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.01/DD0439/DD0439", # 6.630
    # 2B.resim.th.b01-eps-0.001-2
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0267/DD0267", # 6.900
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0292/DD0292", # 6.925
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0305/DD0305", # 6.950
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0317/DD0317", # 6.975
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0330/DD0330", # 7.000
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0342/DD0342", # 7.025
    #"/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eps-0.001-2/DD0355/DD0355", # 7.050
    # # 2B.resim.th.b01-eta-0.1-2 
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0267/DD0267", # 6.900
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0292/DD0292", # 6.925
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0317/DD0317", # 6.950
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0342/DD0342", # 6.975
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0367/DD0367", # 7.000
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0392/DD0392", # 7.025
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0417/DD0417", # 7.050
    # "/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0442/DD0442", # 7.075
    #"/disk01/sgordon/pleiades-18-03-24/seed2-bh-only/270msun/thermal-fb/2B.resim.th.b01-eta-0.1-2/DD0467/DD0467", # 7.100


    ]
    c_min = 1e-21  # g/cm³
    c_max = 1e-13
    num_workers = len(filepaths) + 1
    sim_name = None

    # Run the parallel clump finding
    results = parallel_clump_finding(filepaths, c_min, c_max, num_workers, sim_name=sim_name)
    for result in results:
        print(result)
    logging.info("All clump finding processes completed")
    logging.info("End of script")