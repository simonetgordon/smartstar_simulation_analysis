import yt
import numpy as np
import csv
from datetime import datetime
from yt.fields.field_detector import FieldDetector
from yt.utilities.physical_constants import c
from helper_functions import ss_properties, extract_simulation_name, setup_plot_env, eddington_rate


# Calculate Specific Heating Rate

def _bh_specific_heating_rate(field, data):
    field_data = data.ds.arr(np.zeros(data["gas", "density"].shape), "erg/s/g")
    if isinstance(data, FieldDetector):
        return field_data

    # get all black holes in this data object
    cell_pos = data.ds.arr([data["index", ax] for ax in "xyz"]).T.to("cm")
    bh_pos = data.get_field_parameter("bh_centers").to("cm")
    bh_r = data.get_field_parameter("bh_feedback_radii").to("cm")
    mdot = data.get_field_parameter("accretion_rate").to('g/s')
    epsilon = data.get_field_parameter("epsilon")
    print(f"bh_pos: {bh_pos}, bh_r: {bh_r.to('pc'):.2e}, mdot: {mdot[0].to('Msun/yr'):.2e}, epsilon: {epsilon}")

    # loop over each bh and add heating where d < r
    my_r2 = bh_r**2
    d2 = ((cell_pos - bh_pos)**2).sum(axis=1)
    my_filter = d2 <= my_r2
    field_data[my_filter] = epsilon * c**2 * mdot / data["gas", "cell_mass"] / my_filter.sum()
    # total_mass_in_region = data["gas", "cell_mass"][my_filter].sum()
    # field_data[my_filter] = epsilon * c**2 * mdot / data["gas", "cell_mass"] * (data["gas", "cell_mass"] / total_mass_in_region)

    return field_data

yt.add_field(("gas", "BH_specific_heating_rate"),
             _bh_specific_heating_rate,
             sampling_type="local",
             units="erg/s/g")

def _bh_accretion_cooling_time(field, data):
    # assume here that cooling is positive
    edot = data["gas", "specific_thermal_energy"] / data["gas", "cooling_time"]
    edot_bh = data["gas", "BH_specific_heating_rate"]
    edot_net = edot - edot_bh
    return data["gas", "specific_thermal_energy"] / edot_net

yt.add_field(("gas", "BH_accretion_cooling_time"),
             _bh_accretion_cooling_time,
             sampling_type="local",
             units="s")

def _bh_accretion_heating_to_cooling_ratio(field, data):
    edot = data["gas", "specific_thermal_energy"] / data["gas", "cooling_time"]
    edot_bh = data["gas", "BH_specific_heating_rate"]
    return edot_bh / edot

yt.add_field(("gas", "BH_accretion_heating_to_cooling_ratio"),
             _bh_accretion_heating_to_cooling_ratio,
             sampling_type="local", units="")

def bh_prepare_object(data, epsilon=None):
    spos = data["SmartStar", "particle_position"]
    sdx = data.ds.find_field_values_at_points([("index", "dx")], spos)
    try:
        r_cells = ds.parameters['SmartStarBHThermalFeedbackRadiusinCells']
        if r_cells == 0:
            r_cells = 5
    except KeyError:
        r_cells = 5
    data.set_field_parameter("bh_centers", spos)
    data.set_field_parameter("bh_feedback_radii", r_cells * sdx)
    data.set_field_parameter("epsilon", epsilon)
    mdot_now = data["SmartStar", "AccretionRate"][:, -1]
    data.set_field_parameter("accretion_rate", mdot_now)

# Simulation groups dictionary with paths and parameters
simulation_groups = {
    # 1: {},
    #     # Add more simulations for group 1
    2: {
        "radius_pc": 0.6,
        "simulations": [
            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001", "name": "Initial State", "epsilon": 1e-11,
            # "datadumps": ["DD0445"]},
            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.0001", "name": "eps-0.0001-fbr=5dx", "epsilon": 1e-6, 
            #  "datadumps": ["DD0653", "DD0655", "DD0657", "DD0659", "DD0661", "DD0663", "DD0665", "DD0667", "DD0669", "DD0671", "DD0673", "DD0675", "DD0677", "DD0679", "DD0681"]},
            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eps-0.001", "name": "eps-0.001-fbr=5dx", "epsilon": 1e-5,
            # "datadumps": ["DD2783", "DD2785", "DD2787", "DD2789", "DD2791", "DD2793", "DD2795", "DD2797", "DD2799", "DD2801", "DD2803", "DD2805", "DD2807", "DD2809", "DD2811"]},
            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/fb-radius=10dx", "name": "eta-0.01-fbr=10dx", "epsilon": 5e-4,
            # "datadumps": ["DD0704", "DD0706", "DD0708", "DD0710", "DD0712", "DD0714", "DD0716", "DD0718", "DD0720", "DD0722", "DD0724", "DD0726", "DD0728", "DD0730", "DD0732"]},

            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01/fb-radius=7dx", "name": "eta-0.01-fbr=7dx", "epsilon": 5e-4,
            #  "datadumps": ["DD0696", "DD0698", "DD0700", "DD0702", "DD0704", "DD0706", "DD0708", "DD0710", "DD0712", "DD0714", "DD0716", "DD0718", "DD0720", "DD0722", "DD0724"]},

            # {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.01", "name": "eta-0.01-fbr=5dx", "epsilon": 5e-4,
            # "datadumps": ["DD0641", "DD0643", "DD0645", "DD0647", "DD0649", "DD0651", "DD0653", "DD0655", "DD0657", "DD0659", "DD0661", "DD0663", "DD0665", "DD0667", "DD0669"]},
            {"path": "/ceph/cephfs/sgordon/pleiades-18-03-24/seed1-bh-only/270msun/thermal-fb/1B.resim.th.b01-3-eta-0.1/fb-radius=10dx", "name": "eta-0.1-fbr=10dx", "epsilon": 5e-3,
            "datadumps": ["DD0680", "DD0682", "DD0684", "DD0686", "DD0688", "DD0690", "DD0692", "DD0694", "DD0695"]},
            ],
        },
    # 3: {},
    # 4: {},
}

# CSV setup
csv_file = "specific_heating_rates.csv"
with open(csv_file, mode="a", newline="") as file:
    writer = csv.writer(file)
    #writer.writerow(["event", "simulation_name", "datadump_id", "age_myr", "bh_specific_heating_rate", "mass_enclosed_msun"])

   # Iterate through each simulation group
    for event, group_data in simulation_groups.items():
        radius_pc = group_data["radius_pc"]  # Get the radius for mass calculation
        for sim in group_data["simulations"]:
            for dd_id in sim["datadumps"]:
                try:
                    ds = yt.load(f"{sim['path']}/{dd_id}/{dd_id}")

                    # Add and calculate fields
                    ds.add_field(("gas", "BH_specific_heating_rate"), function=_bh_specific_heating_rate, sampling_type="local", units="erg/s/g", force_override=True)
                    ds.add_field(("gas", "BH_accretion_cooling_time"), function=_bh_accretion_cooling_time, sampling_type="local", units="s", force_override=True)
                    
                    # Black hole feedback sphere
                    pos, mass, age = ss_properties(ds)
                    r_cells = ds.parameters.get('SmartStarBHThermalFeedbackRadiusinCells', 5)
                    dx = ds.find_field_values_at_point([("index", "dx")], pos)[0].to('pc')
                    sphere = ds.sphere(pos, dx * r_cells)
                    
                    # bh_prepare_object
                    sphere.set_field_parameter("bh_centers", pos.to('pc'))
                    sphere.set_field_parameter("bh_feedback_radii", r_cells * dx.to('pc'))
                    sphere.set_field_parameter("epsilon", sim["epsilon"])
                    mdot_now = sphere["SmartStar", "AccretionRate"][:, -1].to('g/s')
                    sphere.set_field_parameter("accretion_rate", mdot_now)
                    
                    # Calculate the heating rate
                    mean_rate = np.mean(sphere[("gas", "BH_specific_heating_rate")])

                    # Enclosed mass calculation within radius_pc
                    enclosed_sphere = ds.sphere(pos, (radius_pc, "pc"))
                    mass_enclosed = enclosed_sphere[("gas", "cell_mass")].sum().to("Msun")
                    
                    # Write to CSV
                    writer.writerow([event, sim["name"], dd_id, age[0]/1e6, mean_rate.d, mass_enclosed.d])

                    # Print progress
                    num_cells = len(sphere[("gas", "BH_specific_heating_rate")])
                    print(f"Processed {sim['name']} {dd_id} at age {age[0]/1e6:.2f} Myr with mean heating rate {mean_rate:.2e} over {num_cells} cells and enclosed mass {mass_enclosed:.2e}")
                
                except Exception as e:
                    print(f"Error processing {sim['name']} {dd_id}: {e}")
                    raise e
