from ocelot.schema.conformer_no_original_partition import MolConformer
from tqdm import tqdm
from ocelot.tasks.newfragment import fragment as f
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

    
 
    retfrags = f(mc, "chrom partition", False) #params for frag are molecule, partition method, and a bool to show all backbones or just BM scaffolds
    #the abolve params mean it will return the fragments with the chrom partition method
    
    identifier = 0
    if retfrags is not None: #None is returned with a fragment error
        for item in retfrags:#the backbone and fragments are returned as a list, which can be iterated over to view contents (probably not ideal)
            identifier+=1 #numerical identification
            if ((not isinstance(item, list)) and( not isinstance(item, tuple ))):
                print(item)
                
            else:
                for obj in item:
                    identifier+=1
                    if ((not isinstance(obj, list)) and( not isinstance(obj, tuple ))):
                        print(obj)
                    else:
                        for other_obj in obj:
                            identifier+=1
                            if ((not isinstance(other_obj, list)) and( not isinstance(other_obj, tuple ))):
                                print(other_obj)


def test_compare(xyzfn):
    # m = Molecule.from_file(xyz)
    # rdmol, smiles = pmgmol_to_rdmol(m)
    name = xyzfn.split('/')[-1][:-4]
    mc = MolConformer.from_file(xyzfn)
    return mc, name


import os
import glob


def testa():
    os.chdir('directory for xyz files')
    os.system('rm *chrom* *geo*')
    mcs: [MolConformer] = []
    for xyzfn in tqdm(glob.glob('*.xyz')):
        mc = test_from_xyz(xyzfn)
        mcs.append(mc)
    os.chdir('../')




testa()
