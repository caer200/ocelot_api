from ocelot.schema.conformer_no_redundance import MolConformer
from tqdm import tqdm
import sys

"""
testing conformer schema:
1. io
2. partition a conjugate molecule into bone and sidechain
3. partition a conjugate molecule into chromophore and sidechain
4. convert a conformer to rdmol, joints are represented by radicals
5. compare rms
"""


def test_from_xyz(xyzfn):
    # m = Molecule.from_file(xyz)
    # rdmol, smiles = pmgmol_to_rdmol(m)
    name = xyzfn.split('/')[-1][:-4]
    mc = MolConformer.from_file(xyzfn)
    #mc.sccs
    mc.scgs
    mc.backbone
    mc.chromobone
    mc.chromobone_graph
    mc.chromsccs
    mc.chromscgs
    mc.molgraph
    mc.sccs
    


def test_compare(xyzfn):
    # m = Molecule.from_file(xyz)
    # rdmol, smiles = pmgmol_to_rdmol(m)
    name = xyzfn.split('/')[-1][:-4]
    mc = MolConformer.from_file(xyzfn)
    return mc, name


import os
import glob


def testa():
    os.chdir('put your directory here')
    os.system('rm *chrom* *geo*')
    mcs: [MolConformer] = []
    for xyzfn in tqdm(glob.glob('*.xyz')):
        mc = test_from_xyz(xyzfn)
        mcs.append(mc)
    os.chdir('../')




testa()
#testb()
