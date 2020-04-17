from ocelot.task.emtensor import *
from pymatgen.io.vasp.outputs import Eigenval, Spin, Outcar, Poscar
import glob
from pprint import pprint
kfc = (0, 0, 0)
iband = 15
ss = 0.01
ss = rb2ra(0.01)
# outcar = Outcar('OUTCAR')
poscar = Poscar.from_file('POSCAR', check_for_POTCAR=False)
real_latt = poscar.structure.lattice.matrix
real_latt = np.array(real_latt)
reci_matrix = 2 * np.pi * np.linalg.inv(real_latt).T  # 1/A
eigens = Eigenval('./EIGENVAL').eigenvalues[Spin.up]
flat_eigens = [x[iband][1] for x in eigens]
flat_eigens = np.array(flat_eigens)

emt = EmTensor(kfc, iband, ss, reci_matrix, flat_eigens, 3)

emt.write_kmesh('KPOINTS_ocelot')
ems, es, eigenvs_frac, eigenvs_cart = emt.cal_emtensor()
pprint(ems)
"""
[[-2.90687151  0.          0.        ]
 [ 0.         -2.90687151  0.        ]
 [ 0.          0.         -2.90687151]]
array([-0.34401245, -0.34401245, -0.34401245])
"""