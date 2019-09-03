from api.schema.omol import OMol
from api.schema.config import Config
from api.routines.pbc import CIFparser, PBCparser
from pymatgen.io.cif import CifParser


with open('../cifs/0.cif', 'r') as f:
    s = f.read()

cp = CIFparser.from_cifstring(s)
d = cp.as_dict()
print(d.keys())
# cp.write_clean_cifs()
# cpp = CifParser('../cifs/0.cif')
# print(cpp.as_dict().keys())
# print(cp.pymatgen_dict.keys())


