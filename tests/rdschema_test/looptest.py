from ocelot.schema.msitelist import MSitelist

msl = MSitelist.from_file('2b.xyz')
sites = msl.msites
# sites.reverse()
msl = MSitelist(sites)
m = msl.rdkit_mol
atoms = m.GetAtoms()
for a in atoms:
    print(a.GetSymbol(), a.GetIdx())