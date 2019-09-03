import sys

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/ai/wd/oscar'])
<<<<<<< HEAD
from api.schema.omol import Mol
=======
from api.schema.Mol import Mol
>>>>>>> 88b3e96e2d37f403133c7972e5e02f735a27c716
from api.tasks.Fitbox import Boxfitter
from api.routines.geometry import abcabg2pbc
import time

ts1 = time.time()

mol = Mol.from_xyz('../sample/ge.xyz')

abcabg_sc_ta = [8.91, 7.51, 17.16, 86.87, 79.75, 65.04]
abcabg_sc_sva = [9.43, 8.24, 12.31, 87.39, 88.43, 93.39]
pbc_sc_ta = abcabg2pbc(abcabg_sc_sva)
bf = Boxfitter(mol, pbc_sc_ta)
bone_configs, config_pqos = bf.gen_bone_configs(steps=9)

for i in range(len(bone_configs)):
    good_configs = bf.grow(bone_configs[i], config_pqos[i], 9, 12, i)
    print('bone {} / {} good config {}'.format(i, len(bone_configs), len(good_configs)))
#    j = 0
#    for c in good_configs:
#         c.to_poscar('bone-{}-{}.poscar'.format(i, j), bf.pbc)
#        j += 1
#    if j > 20:
#        break


# total = len(configs)
# for i in range(len(configs)):
#     good_configs += bf.grow(configs[i], config_pqos[i], 9)
#     print(len(good_configs), i, '/', total)
