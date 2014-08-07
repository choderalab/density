import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as u
from sys import stdout
import sys
import gaff2xml
from gaff2xml.utils import run_antechamber
from gaff2xml.utils import create_ffxml_file
from gaff2xml.utils import convert_molecule
import mdtraj as md
import numpy as np
import simtk.openmm.app.element as elem
import pymbar.timeseries as ts
import rdkit.Chem as chem
from rdkit.Chem import AllChem as allchem


def caseFromPDB(molecule_name, num, steps, ff_mod):
##################################
## How to write a specific case ##
##################################

#    molecule_name = name of molecule
#    num = number of molecules in box
#    ff_mod = True if using charge values from gromacs
#    steps = number of simulation steps (500000)
##################################
    infile, topfile, outfile, outdata, lastframe_file, dcd_file, mol2_filename = name_outputs(molecule_name, num, ff_mod)
    add_bonds_to_pdb(infile, topfile)
    forcefield = generate_ff(molecule_name, infile, mol2_filename, ff_mod, topfile)
    traj, xyz, top = set_topology(infile, num)
    top = assign_element(top)
    system, integrator, temperature = openmm_system(forcefield, top)
    openmm_simulation(system, top, integrator, xyz, temperature, outfile, lastframe_file, outdata, dcd_file, steps)
    avg_density(dcd_file, lastframe_file, outdata)

def caseFromPDBMassMod(molecule_name, num, steps, ff_mod):
    infile, topfile, outfile, outdata, lastframe_file, dcd_file, mol2_filename = name_outputs(molecule_name, num, ff_mod)
    add_bonds_to_pdb(infile, topfile)
    forcefield = generate_ff(molecule_name, infile, mol2_filename, ff_mod, topfile)
    traj, xyz, top = set_topology(infile, num)
    top = assign_element(top)
    system, integrator, temperature = openmm_system(forcefield, top)
    system = modify_mass(system, traj, topfile)
    openmm_simulation(system, top, integrator, xyz, temperature, outfile, lastframe_file, outdata, dcd_file, steps)
    avg_density(dcd_file, lastframe_file, outdata)

def caseFromSmiles(molecule_name, smiles, num, steps, ff_mod):
    infile, topfile, outfile, outdata, lastframe_file, dcd_file, mol2_filename = name_outputs(molecule_name, num, ff_mod)
    files_from_smiles(smiles, infile)
    forcefield = generate_ff(molecule_name, infile, mol2_filename, ff_mod, topfile)
    traj, xyz, top = set_topology(infile, num)
    system, integrator, temperature = openmm_system(forcefield, top)
    system = modify_mass(system, traj, topfile)
    openmm_simulation(system, top, integrator, xyz, temperature, outfile, lastframe_file, outdata, dcd_file, steps)
    avg_density(dcd_file, lastframe_file, outdata)

def name_outputs(molecule_name, num, ff_mod):
    infile = molecule_name+'.pdb'
    topfile = molecule_name+'.top'
    if ff_mod == True:
        outfile = molecule_name+'_'+str(num)+'_out_ff.pdb'
    else:
        outfile = molecule_name+'_'+str(num)+'_out.pdb'
    outdata = outfile[:-4]+'.dat'
    lastframe_file = outfile[:-4]+'_last.pdb'
    dcd_file = outfile[:-4]+'.dcd'
    mol2_filename = molecule_name + '.mol2'
    return infile, topfile, outfile, outdata, lastframe_file, dcd_file, mol2_filename 

def files_from_smiles(smiles, infile):
    molecule = chem.MolFromSmiles(smiles)
    molecule = chem.AddHs(molecule)
    allchem.EmbedMolecule(molecule)
    allchem.UFFOptimizeMolecule(molecule)
    chem.MolToPDBFile(molecule, infile)

def add_bonds_to_pdb(infile, topfile):
    t = md.load_pdb(infile)
    if t.n_residues == 1:
        return
    first = t.top.residue(0)
    t.restrict_atoms(range(first.n_atoms))
    gtop = app.GromacsTopFile(topfile)._currentMoleculeType
    top, bonds = t.top.to_dataframe()
    bonds = np.array([(row[0],row[1]) for row in gtop.bonds],'int')
    bonds = bonds -1
    t.top = md.Topology.from_dataframe(top, bonds)
    t.save(infile)

def generate_ff(molecule_name, infile, mol2_filename, ff_mod, topfile):
    convert_molecule(infile,mol2_filename)
    gaff_mol2_filename, frcmod_filename = run_antechamber(molecule_name, mol2_filename, charge_method="bcc")
    ffxml = create_ffxml_file(gaff_mol2_filename, frcmod_filename)
    forcefield = app.ForceField(ffxml)
    if ff_mod == True:
        ff_mod_funct(forcefield, topfile)
    return forcefield

def ff_mod_funct(forcefield, topfile):
    g = forcefield.getGenerators()[3]
    gtop = app.GromacsTopFile(topfile)._currentMoleculeType
    lock = dict((row[4],float(row[6])) for row in gtop.atoms)
    for key, val in g.typeMap.iteritems():
        f = lambda x: x.split("-")[1]
        lkey = f(key)
        new_charge = lock[lkey]
        g.typeMap[key] = (new_charge, val[1], val[2])

def set_topology(infile, num):
    traj = gaff2xml.packmol.pack_box([infile],[num])
    xyz = traj.openmm_positions(0)
    top = traj.top.to_openmm()
    top.setUnitCellDimensions(mm.Vec3(*traj.unitcell_lengths[0])*u.nanometer)
    return traj, xyz, top

def openmm_system(forcefield, top):
    system = forcefield.createSystem(top, nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer, constraints=app.HBonds)
    temperature = 298.15*u.kelvin
    integrator = mm.LangevinIntegrator(temperature, 1/u.picosecond, 0.002*u.picoseconds)
    barostat = mm.MonteCarloBarostat(1.0*u.atmospheres, temperature, 25)
    system.addForce(barostat)
    return system, integrator, temperature

def openmm_simulation(system, top, integrator, xyz, temperature, outfile, lastframe_file, outdata, dcd_file, steps):
    simulation = app.Simulation(top, system, integrator)
    simulation.context.setPositions(xyz)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(temperature)
    simulation.step(60000)
    simulation.reporters.append(app.PDBReporter(outfile, 1000))
    simulation.reporters.append(app.PDBReporter(lastframe_file,steps-1))
    simulation.reporters.append(app.StateDataReporter(outdata, 1000, step=True, temperature=True, density=True, potentialEnergy=True))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, temperature=True, density=True, potentialEnergy=True))
    simulation.reporters.append(app.DCDReporter(dcd_file,1000))
    simulation.step(steps)

def avg_density(dcd_file,lastframe_file, outdata):
    trj = md.load(dcd_file, top=lastframe_file)
    volume = trj.unitcell_lengths.prod(1)
    mass = sum([a.element.mass for a in trj.top.atoms]) / 6.0221413e23
    density_nounit = mass / volume
    density = density_nounit * u.gram / u.nanometer**3
    A_t = np.array(density_nounit )
    indices = ts.subsampleCorrelatedData(A_t)
    ind_density = density[indices]
    avg_ind_density = ind_density.mean().in_units_of(u.gram / u.liter)
    std_ind_density = ind_density.std().in_units_of(u.gram / u.liter)
    N = len(indices)
    stderr_ind_density = std_ind_density / (N**0.5)
    temps = []
    fid = open(outdata, 'r')
    fid.next()
    for line in fid:
        dtemp = float(line.split(',')[1])
        temps.append(dtemp)
    fid.close()
    temps = np.array(temps)
    avg_temp = temps.mean()
    density_file = 'density_'+dcd_file[:-4]+'_indstd.dat'
    f = open(density_file, 'w')
    f.write("Average density of the system:\n")
    f.write(str(avg_ind_density))
    f.write("\nStandard Deviation of density:\n")
    f.write(str(std_ind_density))
    f.write("\nStandard Error of the density:\n")
    f.write(str(stderr_ind_density))
    f.write("\nAverage Temperature of the system:\n")
    f.write(str(avg_temp))
    f.close()

def modify_mass(system, traj, topfile):
    n_atoms = traj.n_atoms
    n_residues = traj.n_residues
    mlc_atoms = n_atoms / n_residues
    gtop = app.GromacsTopFile(topfile)._currentMoleculeType
    lock = dict((row[4],float(row[7])) for row in gtop.atoms)
    print("\nModified atomic masses:")
    for i in range(mlc_atoms):
        mass = lock[traj.top.atom(i).name]
        for x in np.linspace(i, n_atoms+i, n_residues, endpoint=False):
            system.setParticleMass(int(x), mass *u.dalton)
        print(str(traj.top.atom(i).name)+' '+str(system.getParticleMass(i)))
    print("\n##############################################################\n")
    return system

def assign_element(top):
    for chain in top.chains():
        for res in chain.residues():
            for atom in res.atoms():
                ## this is only going to work as long as there are fewer than 10 of an atom type
                atom.element = elem.get_by_symbol(atom.name[:-1])
    return top

if __name__ == '__main__':
    molecule_name = raw_input('Name of molecule\n')
    num = int(raw_input('Number of molecules in box\n'))
    ff_mod = ('True' == raw_input('Use gromacs charges? (True or False)\n'))
    steps = int(raw_input('Number of simulation steps (Use 500000)\n'))
    case(molecule_name, num, steps, ff_mod)
