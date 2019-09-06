import sys
from schema.config import Config
from routines.pbc import CIFparser
import glob
from task.pkid import PackingIdentifier

sys.path.append('../../')
"""
this reads cif files from tipge-*.cif, identify packing patterns and write backbone-only configs
"""


def inspect_cif(ciffile):
    name = ciffile.split('/')[-1][:-4]
    with open(ciffile, "r") as myfile:
        data = myfile.read()
    cp = CIFparser.from_cifstring(data)
    clean_cif_strings = cp.get_clean_cifs_stringlist()

    c = Config.from_cifstring(clean_cif_strings[0])
    omol = c.omols[0]
    lbone = omol.backbone.lp

    bc, bcpstructure, terminated_backbones = c.get_bone_config()
    bc.unwrap_structure.to('cif', 'boneonlyconfig-{}.cif'.format(name))
    pid = PackingIdentifier(bc)
    packingd = pid.identify_heuristic()
    return name, packingd, lbone


for cif in sorted(glob.glob('./tipge-*.cif')):
    fname = cif.split('/')[-1][:-4]
    n, packingdata, lbone = inspect_cif(cif)
    print("{} {} {:.3f} \n".format(fname, packingdata[0]['packing'], lbone))
