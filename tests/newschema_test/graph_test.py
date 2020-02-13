# molgraph smiles

from ocelot.schema.graph import MolGraph
from ocelot.schema.rdfunc import RdFunc
from rdkit.Chem import Draw

smiles = 'CC[C+](C)Cc3c1ccccc1c(c2ccccc2)c4c[c+]ccc34'
rdmol = RdFunc.from_string('smiles', smiles)
mg = MolGraph.from_rdmol(rdmol)
# rdmol, s1, _, _ = mg.to_rdmol(charge=2, charged_fragments=True)
# tmp=AllChem.Compute2DCoords(rdmol)
Draw.MolToFile(rdmol,'m.ps')
mg.draw('m.eps')



