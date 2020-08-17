
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Scaffolds import MurckoScaffold
from ocelot.schema.graph_without_partition import MolGraph
#from ocelot.schema.conformer_no_original_partition import MolConformer as fragment_helper
from ocelot.tasks.newfragment import fragment as f



TEST_SMILES = ["C1(=CC5=C(C2=C1C(=C3C(=C2)C=CC=C3CCCCC4=CC=CC=C4)CCCC)C=CC=C5)CCCC"]
def draw_rdmol(rdmol, name):
    print(name)
    # rdmol = AddHs(rdmol)
    AllChem.Compute2DCoords(rdmol)
    Draw.MolToFile(rdmol, '{}.png'.format(name))
    AllChem.EmbedMolecule(rdmol)
    print(Chem.MolToMolBlock(rdmol), file=open('{}.mol'.format(name), 'w+'))


def test_smiles(smiles: str, name):
    name = "{}".format(name)

    #draw_rdmol(rm, name + '_rdmol')

    mg = MolGraph.from_smiles(smiles)
    mg.draw("original molecule")
    bg = f(mg, "partition", False) #params for frag are molecule, partition method, and a bool to show all backbones or just BM scaffolds
    #the abolve params mean it will return the partitioned graph with only the BM scaffold backbone, and fragments
    
    identifier = 0
    if bg is not None: #None is returned with a fragment error
        for item in bg:#the backbone and fragments are returned as a list, which can be iterated over to view contents (probably not ideal)
            identifier+=1 #numerical identification
            if ((not isinstance(item, list)) and( not isinstance(item, tuple ))):
                item.draw("{}".format(identifier))
            else:
                for obj in item:
                    identifier+=1
                    if ((not isinstance(obj, list)) and( not isinstance(obj, tuple ))):
                        obj.draw("{}".format(identifier))
                    else:
                        for other_obj in obj:
                            identifier+=1
                            if ((not isinstance(other_obj, list)) and( not isinstance(other_obj, tuple ))):
                                other_obj.draw("{}".format(identifier))



import os
def testa():
    name = 0
    
    for smiles in TEST_SMILES:

        test_smiles(smiles, name)
        name += 1
   
testa()
