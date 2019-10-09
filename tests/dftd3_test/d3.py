from pymatgen.core.structure import Structure
from ocelot.task.dftd3 import DFTD3

d3cmd = 'dftd3'
structure = Structure.from_file('POSCAR')
d3 = DFTD3('dummy', structure)
result = d3.run(d3cmd, './')
r = d3.parse_results(result)
print(r)
