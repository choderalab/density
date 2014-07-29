import density_sim as ds

molecule_name = 'ethanol'
num = 200
steps = 500000
ff_mod = True

ds.case(molecule_name, num, steps, ff_mod)
