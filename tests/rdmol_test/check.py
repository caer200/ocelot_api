from ocelot.schema.conformer import Molecule
from ocelot.routines.conformerparser import pmgmol_to_rdmol
from ocelot.schema.rdfunc import RdFunc

# rdkit is quite tricky in handling nitro groups
# see https://github.com/rdkit/rdkit/issues/1979
# and http://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization
m = Molecule.from_file('m.xyz')
rm, smiles = pmgmol_to_rdmol(m)
RdFunc.draw_mol(rm, 'm.ps')
print(rm)