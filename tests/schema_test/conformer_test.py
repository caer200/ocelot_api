from ocelot.schema.conformer import MolConformer
from tqdm import tqdm

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
    mc.backbone.to('xyz', '{}_geobone.xyz'.format(name))
    mc.sccs[0].to('xyz', '{}_geosc0.xyz'.format(name))
    bh, scchs = mc.to_addhmol('geo')
    bh.to('xyz', '{}_geoboneh.xyz'.format(name))
    scchs[0].to('xyz', '{}_geosc0h.xyz'.format(name))
    mc.chrombone.to('xyz', '{}_chrombone.xyz'.format(name))
    mc.chromsccs[0].to('xyz', '{}_chromsc0.xyz'.format(name))
    bh, scchs = mc.to_addhmol('chrom')
    bh.to('xyz', '{}_chromboneh.xyz'.format(name))
    scchs[0].to('xyz', '{}_chromsc0h.xyz'.format(name))
    return mc


def test_compare(xyzfn):
    # m = Molecule.from_file(xyz)
    # rdmol, smiles = pmgmol_to_rdmol(m)
    name = xyzfn.split('/')[-1][:-4]
    mc = MolConformer.from_file(xyzfn)
    return mc, name


import os
import glob


def testa():
    os.chdir('./ixyzs')
    os.system('rm *chrom* *geo*')
    mcs: [MolConformer] = []
    for xyzfn in tqdm(glob.glob('*.xyz')):
        mc = test_from_xyz(xyzfn)
        mcs.append(mc)
    os.chdir('../')


def testb():
    os.chdir('./ixyzs')
    os.system('rm *chrom* *geo*')
    mcs: [MolConformer] = []
    names = []
    for xyzfn in glob.glob('*.xyz'):
        mc, name = test_compare(xyzfn)
        mcs.append(mc)
        names.append(name)
    os.chdir('../')
    import time
    for i in range(len(mcs)):
        for j in range(i, len(mcs)):
            ts1 = time.time()
            rms = mcs[i].compare(mcs[j])
            ts2 = time.time()
            print(rms, "{} vs {}".format(names[i], names[j]), ts2 - ts1)

testa()
# testb()
