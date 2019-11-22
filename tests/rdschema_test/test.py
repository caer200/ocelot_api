from ocelot.schema.rdschema import OmolRd
from ocelot.schema.msitelist import MSitelist
msl = MSitelist.from_file('2b.xyz')
m = msl.rdkit_mol
rdomol = OmolRd(m)

print(rdomol.scs)
print(rdomol.scs[0].rankmap)
print(rdomol.scs[0].atoms)
print(rdomol.backbone.unsat)
print(rdomol.backbone)
print(rdomol.rings)
