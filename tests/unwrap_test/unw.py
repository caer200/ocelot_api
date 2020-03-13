from ocelot.routines.pbc import PBCparser, Structure
import glob
import datetime

# patt = 'c17013.cif'
# patt = '*.cif'
# patt = 'k06071.cif'
# patt = 'x15029.cif'
patt = 'dum.cif'
# patt = 'k03017.cif'
# patt = 'k14136.cif'

for fn in glob.glob(patt):
    structure = Structure.from_file(fn)
    identifier = fn[:-4]
    print('{} try to read {}'.format(datetime.datetime.now(), fn))
    _, unwrap, _ = PBCparser.unwrap(structure)
    unwrap.to('cif', 'unw.cif')
