from ocelot.schema.conformer import DimerCollection
from ocelot.schema.configuration import Config
from ocelot.routines.disparser import DisParser
from ocelot.task.hop import Hop

def inspect_cif(ciffile):
    dp = DisParser.from_ciffile(ciffile)
    confs = dp.to_configs()
    conf = confs[0]
    return Config.from_labeled_clean_pstructure(conf)


# fname = './tipge-bw.cif'
# c = inspect_cif(fname)
# dimer_arrays, trans = c.get_dimers_array(maxfold=1, fast=False, symm=False)
# dc = DimerCollection(dimer_arrays[0][0])
# dc.to_xyz('dc.xyz')

# fname = 'POSCAR'
# c = Config.from_file(fname)
#c = inspect_cif('tipge-bw.cif')
c = inspect_cif('86.cif')

hp = Hop(config=c)
dimers = hp.screen_dimers(hp.dimer_array)
for d in dimers:
    print(d.pslip, d.qslip, d.oslip)

c = inspect_cif('100.cif')

hp = Hop(config=c)
dimers = hp.screen_dimers(hp.dimer_array)
for d in dimers:
    print(d.pslip, d.qslip, d.oslip)

# dc = DimerCollection(dimers)
# for d in dc.dimers:
#     print(d)
