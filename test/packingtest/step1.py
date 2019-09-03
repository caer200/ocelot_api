from api.schema.omol import OMol
from api.schema.config import Config
from api.routines.pbc import CIFparser, PBCparser
from pymatgen.core.structure import Structure
import os
import glob
from api.task.packing import PackingIdentifier

def dostuff(cif):
    name = cif.split('/')[-1][:-4]
    with open (cif, "r") as myfile:
        data=myfile.read()
    cp = CIFparser.from_cifstring(data)
    cp.write_clean_cifs()

    c = Config.from_file('cleanconfig-0.cif')
    omol = c.omols[0]
    lbone = omol.backbone.lp
    dsidegroup = max([sc.get_sphere_diameter() for sc in omol.sidechains])

    bc, bcpstructure, terminated_backbones= c.get_bone_config()
    bc.unwrap_structure.to('cif', 'boneonlyconfig-{}.cif'.format(name))
    pid = PackingIdentifier(bc)
    packingdata = pid.identify_heuristaic()



    return name, packingdata, lbone, dsidegroup

# for cif in ['ss.cif', 'sh.cif', 'rb.cif']:
# f = open('packing.txt', 'w')
# f.write('# cifname packing  n_close_azm_and_parallel  n_close_azm_and_notparallel  n_close_vertical_and_parallel  n_close_vertical_and_notparallel  n_parallel_and_overlap  n_notparallel_and_overlap  \n')
# f.close()
for cif in sorted(glob.glob('/home/ai/mnt/dlx4/polytip/raw/tipge-*.cif')):
# for cif in sorted(glob.glob('/home/ai/wd/oscar/cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
    print('working on: ' + cif)
    # name = int(cif.split('/')[-1][:-4])
    name = cif.split('/')[-1][:-4]
    # if name == 47:
    if name != None:
        n, packingdata, lbone, dsidegroup = dostuff(cif)
        f = open('packing-tipge.txt', 'a')
        f.write("{} {} {:.3f} {:.3f}\n".format(name, '-'.join(packingdata), lbone, dsidegroup))
        # for pk in range(len(packingdata)):
            # packing = packingdata[pk]['packing']
            # keys = packingdata[pk].keys()
            # f.write("{} {} {}\n".format(n, packing, '-'.join([str(packingdata[pk][kk]) for kk in keys])))
        f.close()




# c.unwrap_structure.to('cif', 'unwrap.cif')
# dimers = bc.dimers
# dimer = dimers[2]
# dimer.to_xyz('dimer2.xyz')
# print(dimer.bone_overlap())

# terminated_backbones[0].to('xyz', 'bone0.xyz')
# terminated_backbones[1].to('xyz', 'bone1.xyz')

# bc.unwrap_structure.to('cif', 'bconfig.cif')
# dimers = c.dimers
# for i in range(len(dimers)):
#     print(dimers[i].jmol, dimers[i].label)
#     dimers[i].to_xyz('d_{}.xyz'.format(i))

