from ocelot.schema.configuration import *
from ocelot.schema.rdfunc import *
# TODO this doesnt work as one would expect, see todo in conformerparser

config = Config.from_file('test.cif')

mc = config.molconformers[0]
# mc.to('xyz', 'test.xyz')
rdmol: Chem.Mol
rdmol, smiles, siteid2atomidx, atomidx2siteid = mc.to_rdmol()
print(len(rdmol.GetAtoms()))
print(Descriptors.NumRadicalElectrons(rdmol))
for a in rdmol.GetAtoms():
    radicals = AllChem.Atom.GetNumRadicalElectrons(a)
    if radicals > 0 :
        AllChem.Atom.SetNumRadicalElectrons(a, 0)
        AllChem.Atom.SetNumExplicitHs(a, radicals)
# print(smiles)
rdmolh: Chem.Mol = Chem.AddHs(rdmol, addCoords=True)
print(len(rdmolh.GetAtoms()))
# print(Chem.MolToSmiles(rdmolh))
RdFunc.mol2xyz_by_confid(rdmolh, 'testrd', 0)

