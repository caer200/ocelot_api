from ocelot.task.emtensor import *

emt = EmTensor.from_vasprun((0, 0, 0), 118, rb2ra2pi(0.035), 'vasprun.xml', st=3)
emt.write_kmesh('KPOINTS_ocelot')
print(len(emt.kmesh))
ems, es, eigenvs_frac, eigenvs_cart = emt.cal_emtensor()
from pprint import pprint


pprint(ems)
pprint(eigenvs_cart)
pprint(eigenvs_frac)
