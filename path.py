from dimers_pair import *
from ocelot.routines.geometry import Fitter
from ocelot.schema.config import Config


config = Config.from_file("POSCAR")
dimer_array = config.get_dimers_array(maxfold=2, fast=True)
mols = len(config.omols)
dimer_unique = get_unique_dimers(dimer_array, mols)

couplings = run_zindo(dimer_unique, '/home/vinayak/OCELOT/ocelot_api/')
couplings = [abs(val) for val in couplings]
hop_vector = dimer_unique[couplings.index(max(couplings))][1]
print("Hopping along", hop_vector)
com, mol_id = hopping(dimer_array,hop_vector)
vector, mean, error = Fitter.linear_fit(com)
print('The error is', error)
print('The vector is', vector)
print('The mean is', mean)