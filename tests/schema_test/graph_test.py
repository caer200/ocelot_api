import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

from ocelot.schema.graph import MolGraph

"""
testing graph schema:
1. io
2. partition a conjugate molecule into bone and sidechain
3. partition a conjugate molecule into chromophore and sidechain
4. convert a graph to rdmol, joints (if any) are represented by radicals
"""
TEST_SMILES = [
    "c1ccc(c2ccccc2)cc1",
    "Oc1cccc(c2ccccc2)c1",
    "O=C1C=CC(=O)C(=C1)c1ccccc1",
    "O=C1C=C(c2ccccc2)C(=O)c2sccc12",
    "C1=CC(c2ccccc2)C=C1",
    "Clc1ccc2ccccc2c1",
    "C=C/C(=C\C)/CCc1cc2ccccc2cc1Cl",
    "C=Nc1cc(OCC)cc2cc(C#N)c(F)cc12"
]


def draw_rdmol(rdmol, name):
    # rdmol = AddHs(rdmol)
    AllChem.Compute2DCoords(rdmol)
    Draw.MolToFile(rdmol, '{}.png'.format(name))
    AllChem.EmbedMolecule(rdmol)
    print(Chem.MolToMolBlock(rdmol), file=open('{}.mol'.format(name), 'w+'))


def test_smiles(smiles: str, name):
    name = "{}".format(name)
    rm = Chem.MolFromSmiles(smiles)
    print('draw rdmol...')
    draw_rdmol(rm, name + '_rdmol')

    mg = MolGraph.from_smiles(smiles)
    schemes = ['lgfr', 'lgcr', 'chrom']
    for scheme in schemes:
        bg, scgs = mg.partition_to_bone_frags(scheme)
        print('draw {} bone graph...'.format(scheme))
        bg.draw('{}_bonegraph_{}.eps'.format(name, scheme))

        bgrdmol, bgsmiles, _, _ = bg.to_rdmol()
        print('draw {} bone rdmol...'.format(scheme))
        draw_rdmol(bgrdmol, '{}_bonemol_{}'.format(name, scheme))

        scrdmol, scsmiles, _, _ = scgs[0].to_rdmol()
        print('draw {} sc0...'.format(scheme))
        draw_rdmol(scrdmol, '{}_sc0_{}'.format(name, scheme))


import os
def testa():
    name = 0
    os.chdir('./graph_results')
    os.system('rm *.mol *.eps *.png')
    for smiles in TEST_SMILES:
        test_smiles(smiles, name)
        name += 1
    os.chdir('../')

testa()

