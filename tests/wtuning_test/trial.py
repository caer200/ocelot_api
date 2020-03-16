from pymatgen.core.structure import Molecule
from sys import argv
from ocelot.task.wtuning import WtuningJob
import os

fn = argv[1]
name = fn.split('.')[0]
iwd = os.getcwd()
pymol = Molecule.from_file(fn)
os.mkdir(name)
os.chdir(str(iwd)+'/'+name)
Mol = WtuningJob(func='uLC-wHPBE', basis='6-31G', name=name, nproc=16, mem=50,
                 n_charge=0, n_spin=1, wdir='./', scheme='Jh', wbmin=0.05, wbmax=0.5)
Mol.mol=pymol
Mol.wtuning_cycle()
