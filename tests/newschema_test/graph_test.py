
# test 1 nodecollection
from ocelot.schema.graph.elementnode import ElementNode
from ocelot.schema.graph.nodecollection import NodeCollection

print('test 1')
n1 = ElementNode('bq', 2)
n2 = ElementNode('bq', 3)
n3 = ElementNode('bq', 4)
n4 = ElementNode('ck', 2)
nc = NodeCollection([n2, n1, n4])
print(nc)
print(n2 in nc)  # True
print(n3 in nc)  # False
print(n4 in nc)  # False

# test 2 molgraph smiles

from ocelot.schema.graph.molgraph import MolGraph
from ocelot.schema.rdfunc import RdFunc
from rdkit.Chem import Draw

smiles = 'CC[C+](C)Cc3c1ccccc1c(c2ccccc2)c4c[c+]ccc34'
rdmol = RdFunc.from_string('smiles', smiles)
mg = MolGraph.from_rdmol(rdmol)
# rdmol, s1, _, _ = mg.to_rdmol(charge=2, charged_fragments=True)
# tmp=AllChem.Compute2DCoords(rdmol)
Draw.MolToFile(rdmol,'m.ps')
mg.draw('m.eps')



