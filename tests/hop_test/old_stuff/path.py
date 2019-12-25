# This is an example file to compute the zindo coupling. The location of the ZINDOBIN must be set in dimers_pair.py

from dimers_pair import *
from ocelot.routines.geometry import Fitter
from ocelot.schema.config import Config


config = Config.from_file("POSCAR") # load the structre from file
dimer_array = config.get_dimers_array(maxfold=2, fast=True) # get the dimer array
mols = len(config.omols)  # find the numeber of molecules in unit cell
dimer_unique = get_unique_dimers(dimer_array, mols) # get unique dimer pairs for zindo calculation

couplings = run_zindo(dimer_unique, '/home/vinayak/OCELOT/ocelot_api/') # set the work dir here
couplings = [abs(val) for val in couplings]
hop_vector = dimer_unique[couplings.index(max(couplings))][1] # specify the hopping vector
print("Hopping along", hop_vector)
com, mol_id = hopping(dimer_array,hop_vector)
vector, mean, error = Fitter.linear_fit(com)
print('The error is', error)
print('The vector is', vector)
print('The mean is', mean)
