"""
calculate indo/1 transfer integral
"""
from ocelot.routines.zindo import ZindoJob
from pymatgen.core.structure import Molecule

ZINDOLIB = '/home/ai/ZINDO/dummy/'
ZINDOBIN = '/home/ai/ZINDO/dummy/zindo'
ZINDOCTBIN = '/home/ai/ZINDO/dummy/zindo-ct'
wdir = './'
mol_A = Molecule.from_file('a.xyz')
mol_D = Molecule.from_file('d.xyz')

data = ZindoJob.dimer_run('test', './', ZINDOBIN, ZINDOCTBIN, ZINDOLIB, mol_A, mol_D)
print(data)

