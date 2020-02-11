from ocelot.routines.disparser import DisParser
from ocelot.task.pkid import PackingIdentifier
from ocelot.schema.configuration import Config
import glob
"""
this reads cif files from tipge-*.cif, identify packing patterns and write backbone-only configs
"""


def pkid_ciffile(ciffile):
    name = ciffile.split('/')[-1][:-4]
    dp = DisParser.from_ciffile(ciffile)
    dp.to_configs(write_files=True)  # writes conf_x.cif

    # backbone config
    config = Config.from_file('conf_0.cif')
    bc, boneonly_pstructure, terminated_backbone_hmols = config.get_bone_config()
    bc.unwrap_structure.to('cif', 'boneonlyconfig-{}.cif'.format(name))
    pid = PackingIdentifier(bc)
    packingd = pid.identify_heuristic()
    return packingd

for cif in sorted(glob.glob('./tipge-*.cif')):
    fname = cif.split('/')[-1][:-4]
    packingdata = pkid_ciffile(cif)
    print("{} {}\n".format(fname, packingdata[0]['packing']))
