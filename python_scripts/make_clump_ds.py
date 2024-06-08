import logging
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
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
    future_bound = "True"

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
        fn = master_clump.save_as_dataset(filename=output_dir, fields=["density", ("gas", "H2_fraction"), ("gas", "cell_mass"), ("gas", "jeans_mass"), ("gas", "number_density"), ("gas", "temperature")])
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
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0450/DD0450",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0565/DD0565",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0595/DD0595",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0612/DD0612",
        "/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0655/DD0655",
        #"/disk01/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/DD0695/DD0695"
    ]
    c_min = 1e-18  # g/cm³
    c_max = 1e-13
    num_workers = len(filepaths) + 1

    # Run the parallel clump finding
    results = parallel_clump_finding(filepaths, c_min, c_max, num_workers)
    for result in results:
        print(result)
