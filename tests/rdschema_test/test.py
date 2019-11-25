from ocelot.schema.rdschema import OmolRd
from ocelot.schema.msitelist import MSitelist
msl = MSitelist.from_file('2b.xyz')
m = msl.rdkit_mol
rdomol = OmolRd(m)
import pprint

d = rdomol.as_dict()
pprint.pprint(d['rings'])
