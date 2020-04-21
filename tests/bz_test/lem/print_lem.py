from ocelot.task.bzcal import *



# DispersionRelationVASP.read_vasprun('vasprun.xml')
dr = DispersionRelationLine.from_klines_files('vasprun.xml', 'KPOINTS')
# dr.plotlinebs()
ems = dr.get_line_ems('vb')
from pprint import pprint
pprint(ems)
