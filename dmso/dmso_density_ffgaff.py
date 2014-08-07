import density_sim as ds

molecule_name = 'dmso'
num = 200
steps = 500000
ff_mod = False

ds.caseFromPDB(molecule_name, num, steps, ff_mod)
