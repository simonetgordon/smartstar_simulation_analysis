from make_clump_ds import *

# Load the datasets
fp = "/disk14/sgordon/pleiades-11-12-23/seed1-bh-only/270msun/replicating-beckmann/1B.RSb01-2/DD0446/DD0446"
ds = yt.load(fp)
master_clump = yt.load("DD0446_clump_0.h5")
sim = extract_simulation_name(fp)
dd = extract_dd_segment(fp)
leaf_clumps = master_clump.leaves

# Find sink particle attributes
ss_pos, ss_mass, ss_age = ss_properties(ds, velocity=False)
print("Sink particle properties:")
print("position = {}, mass = {:.2f}, age = {:.2f} Myr".format(ss_pos, ss_mass, ss_age[0]/1e6))

# Make a projection plot of the clump leaves
for i in range(3):
    p = yt.ProjectionPlot(ds, i, ("gas", "number_density"), center=ss_pos, width=(6, "pc"))
    p.set_cmap('number_density', 'octarine')
    p.set_zlim('number_density', 4e21, 8e25)
    p.annotate_timestamp(corner='upper_right', redshift=True, draw_inset_box=True)
    p.annotate_clumps(leaf_clumps, cmap="Pastel1")
    #p.annotate_clumps([most_massive_clump], text_args={"color": "red"})
    #most_massive_clump_mass = most_massive_clump.info["cell_mass"][1]
    #most_massive_clump_pos = most_massive_clump.info["position"][1].to('pc')
    #text_string = f"Most Massive Clump: {most_massive_clump_mass:.2f}\n Total Clump Mass: {sum(leaf_masses):.0f}\n Number of Clumps: {len(leaf_clumps)}\n"
    #p.annotate_text((0.02, 0.85), text_string, coord_system='axis', text_args={'color': 'white'})
    #p.annotate_text([0.05, 0.05], sim, coord_system="axis", text_args={"color": "black"}, 
     #               inset_box_args={"boxstyle": "square,pad=0.3", "facecolor": "white", "linewidth": 3, 
     #                               "edgecolor": "white", "alpha": 0.5}
     #               )
    p.annotate_text((0.78, 0.05), "BH Age = {:.2f} Myr".format(ss_age[0]/1e6), coord_system="axis", text_args={"color": "white"})
    #p.annotate_marker(master_clump.info["position"][1], coord_system='data', color='white', s=100)
    p.save(f"plots/clump_projection_{sim}_{dd}_{i}_s.png")
    print(f"Saved plots/clump_projection_{sim}_{dd}_{i}_s.png")

for i in range(3):
    print("created density_projection_{}_{}_s.png".format(i, sim))