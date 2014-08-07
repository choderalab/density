import density_sim as ds

molecule_name = 'dmso'
num = 200
steps = 500000
ff_mod = True

infile, topfile, outfile, outdata, lastframe_file, dcd_file, mol2_filename = ds.name_outputs(molecule_name, num, ff_mod)
ds.add_bonds_to_pdb(infile, topfile)
forcefield = ds.generate_ff(molecule_name, infile, mol2_filename, ff_mod, topfile)
traj, xyz, top = ds.set_topology(infile, num)
system, integrator, temperature = ds.openmm_system(forcefield, top)
system = ds.modify_mass(system, traj, topfile)
ds.openmm_simulation(system, top, integrator, xyz, temperature, outfile, lastframe_file, outdata, dcd_file, steps)
ds.avg_density(dcd_file,lastframe_file)

