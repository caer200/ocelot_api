from ocelot.schema.config import Config, np
from ocelot.routines.pbc import CIFparser
from ocelot.schema.dimercollection import DimerCollection

def inspect_cif(ciffile):
    with open(ciffile, "r") as myfile:
        data = myfile.read()
    cp = CIFparser.from_cifstring(data)
    clean_cif_strings = cp.get_clean_cifs_stringlist()

    c = Config.from_cifstring(clean_cif_strings[0])
    return c


# fname = './tipge-bw.cif'
# c = inspect_cif(fname)
# dimer_arrays, trans = c.get_dimers_array(maxfold=1, fast=False, symm=False)
# dc = DimerCollection(dimer_arrays[0][0])
# dc.to_xyz('dc.xyz')

fname = 'POSCAR'
c = Config.from_file(fname)
dimer_arrays, trans = c.get_dimers_array(maxfold=1, fast=False, symm=False)
dc = DimerCollection(dimer_arrays[0][0])
dc.to_xyz('dc.xyz')
