from dimers_pair import *

config = Config.from_file("POSCAR")
dimer_array = config.get_dimers_array(maxfold=2, fast=True)
#dimer_close = get_close_molecules(dimer_array,1)
#write_dimer_close(dimer_close)
print(dimer_array[1][55])