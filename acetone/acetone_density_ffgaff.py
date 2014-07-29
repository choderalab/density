import density_sim as ds

molecule_name = 'acetone'
num = 200
steps = 500000
ff_mod = False

ds.case(molecule_name, num, steps, ff_mod)
