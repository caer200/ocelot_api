from ocelot.task.emtensor import EmTensor

emt = EmTensor.from_vasprun((0, 0, 0), 118, 0.035, 'vasprun.xml', st=3)
print(len(emt.kmesh))
ems, es, eigenvs_frac, eigenvs_cart = emt.cal_emtensor()
from pprint import pprint

pprint(ems)
pprint(eigenvs_cart)
pprint(eigenvs_frac)
