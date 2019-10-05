import sys
sys.path.append('../../')
from ocelot.schema import Config
from ocelot.routines import CIFparser
import glob
from ocelot.task.pkid import PackingIdentifier
"""
write pk params to 'packing.txt' for all cifs in ../sample/cifs/
'./japacking.txt' shows the packing pattern identified by JA manually
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
    pid = PackingIdentifier(bc)
    packingd = pid.identify_heuristic()
    return name, packingd, lbone

f = open('packing.txt', 'w')
f.write('# cifname packing  n_close_azm_and_parallel  n_close_azm_and_notparallel  n_close_vertical_and_parallel  n_close_vertical_and_notparallel  n_parallel_and_overlap  n_notparallel_and_overlap  \n')
for cif in sorted(glob.glob('../sample/cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
    print('working on: ' + cif)
    n, packingdata, lbone = inspect_cif(cif)
    packingdata_list = [p for p in packingdata.values()]
    outstring = ''
    for refp in packingdata_list:
        packingstring = '-'.join([str(va) for va in refp.values()])
        outstring += packingstring + '@'
    f.write("{}\t\t{}\t\t{:.3f} \n".format(n, outstring, lbone))
f.close()
