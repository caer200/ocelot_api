from Routines import *
from Schema import *
from Gauss import *
import glob
import os

for inp in glob.glob('./optgeo/*.gjf'):
    mol = Mol.from_gjf(inp)
    name = inp.split('/')[-1].split('.')[0]
    os.chdir('./nicsxyinp/')
    gp = GaussProj(mol, 'nicsxy-backbone-total', name, addkwds=mol.comment[2:], level='')
    gp.write_inp({'step_size':0.2, 'height':1.7, 'maxnbq':50, 'nrings':3, 'normaldirection':0})
    gp = GaussProj(mol, 'nicsxy-backbone-sigma', name, addkwds=mol.comment[2:], level='')
    gp.write_inp({'step_size':0.2, 'height':1.7, 'maxnbq':50, 'nrings':3, 'normaldirection':0})
    os.chdir('../')

