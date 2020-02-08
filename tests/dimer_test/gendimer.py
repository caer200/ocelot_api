from ocelot.schema.conformer import DimerCollection
from ocelot.schema.configuration import Config
from ocelot.routines.disparser import DisParser

def inspect_cif(ciffile):
    dp = DisParser.from_ciffile(ciffile)
    confs = dp.to_configs()
    conf = confs[0][0]
    return Config(conf)


# fname = './tipge-bw.cif'
# c = inspect_cif(fname)
# dimer_arrays, trans = c.get_dimers_array(maxfold=1, fast=False, symm=False)
# dc = DimerCollection(dimer_arrays[0][0])
# dc.to_xyz('dc.xyz')

# fname = 'POSCAR'
# c = Config.from_file(fname)
c = inspect_cif('tipge-bw.cif')

dimer_arrays, trans = c.get_dimers_array(maxfold=1, fast=True, symm=False)
dc = DimerCollection(dimer_arrays[0][0])
dc.to_xyz('dc.xyz')
