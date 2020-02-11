from ocelot.routines.disparser import DisParser
from ocelot.task.pkid import PackingIdentifier
from ocelot.schema.configuration import Config
import glob
"""
this reads cif files from tipge-*.cif, identify packing patterns and write backbone-only configs
"""


def pkid_ciffile(ciffile):
    dp = DisParser.from_ciffile(ciffile)
    dp.to_configs(write_files=True)  # writes conf_x.cif

    # backbone config
    config = Config.from_file('conf_0.cif')
    bc, boneonly_pstructure, terminated_backbone_hmols = config.get_bone_config()
    boneonly_pstructure.to('cif', 'boneonly.cif')

    pid = PackingIdentifier(bc)
    packingd = pid.identify_heuristic()
    return packingd




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
