from rdkit.Chem import AllChem
import rdkit.Chem as Chem
from ocelot.schema.graph import MolGraph
from pymatgen.core.structure import Molecule
from ocelot.routines.conformerparser import pmgmol_to_rdmol
from rdkit.Chem import Draw
from rdkit.Chem.rdmolops import AddHs

def draw(rdmol, name):
    rdmol = AddHs(rdmol)
    AllChem.Compute2DCoords(rdmol)
    Draw.MolToFile(rdmol, '{}.png'.format(name))
    AllChem.EmbedMolecule(rdmol)
    print(Chem.MolToMolBlock(rdmol),file=open('{}.mol'.format(name),'w+'))


def parse(fn):
    print('--- start {}'.format(fn))
    name = fn[:-4]
    m = Molecule.from_file(fn)
    rdmol, smiles = pmgmol_to_rdmol(m)
    print('molecule smiles: {}'.format(smiles))
    om = MolGraph.from_rdmol(rdmol)

    omrdmol, _, _, _ = om.to_rdmol()
    draw(omrdmol, name)
    bone, scs = om.partition_to_bone_frags()
    bone_rdmol, bone_smiles, _, _ = bone.to_rdmol()
    print('bone smiles: {}'.format(bone_smiles))
    draw(bone_rdmol, name + '-bone')
    i = 0
    for sc in scs:
        if len(sc.graph) == 1:
            continue
        rdmol, smiles, _, _ = sc.to_rdmol()
        print('sc #{} smiles: {}'.format(i, smiles))
        draw(rdmol, name + '-sc{}'.format(i))
        print('sc #{} in gm? {}'.format(i, sc.is_subgraph(om)))
        i += 1
    print('bone in gm? {}'.format(bone.is_subgraph(om)))
    print('bone graph iso gm? {}'.format(bone == om))
    return om

omdict = {}
for i in [
    # 'bp.xyz',
    # 'penfrag.xyz',
    # 'rubbish.xyz',
    # 'tipg.xyz',
    'tipgpn.xyz',
    # 'tipgpn1.xyz',
    # 'rub.xyz'
]:
    omdict[i] = parse(i)

# print(omdict['tipgpn.xyz'] == omdict['tipgpn1.xyz'])
