from ocelot.routines.pbc import CIFparser

with open('tipgebw.cif', 'r') as f:  # read cif file into string
    fs = f.read()

cp = CIFparser.from_cifstring(fs)  # create a parser object from cif string
clean_cif_strings = cp.get_clean_cifs_stringlist()  # a list of cif strings without disorder
with open('tipgebw_clean_0.cif', 'w') as f:  # write the first cleaned cif string into file
    f.write(clean_cif_strings[0])


from ocelot.schema.config import Config

# init Config object from cleaned cif string
tipge_config = Config.from_cifstring(clean_cif_strings[0])

# we can remove pbc and focus on the organic molecule in this configuration
# OMol -- organic molecule -- is a list of MSites
omol_0 = tipge_config.omols[0]
print(omol_0)
# OMol:
# MSite: Ge (1.8797, 1.8744, 14.3509) siteid: 0
# MSite: C (2.7783, 2.5445, 12.7950) siteid: 1
# MSite: C (0.0097, 2.4446, 14.2340) siteid: 2
# ...

# most of the time we just care conjugate backbone, is also a list of MSites
bone = omol_0.backbone
print(bone)
# Backbone:
# MSite: C (3.8686, 3.3564, 10.5508) siteid: 15
# MSite: C (4.5668, 2.4190, 9.7409) siteid: 34
# MSite: C (3.7152, 4.7050, 10.1166) siteid: 35
# ...

# we can also check sidechains, notice even a single H is considered as one sidechain
# again, sidechain is a list of MSites
sidechains = omol_0.sidechains
print(len([sc for sc in sidechains if not sc.ishydrogen]))  # how many non-H side chains?
# 2










#
#
# print()
#
# lbone = omol_0.backbone.lp




# import sys
# from schema.config import Config
# from routines.pbc import CIFparser
# import glob
# from task.pkid import PackingIdentifier
#
# sys.path.append('../../')
# """
# this reads cif files from tipge-*.cif, identify packing patterns and write backbone-only configs
# """
#
#
# def inspect_cif(ciffile):
#     name = ciffile.split('/')[-1][:-4]
#     with open(ciffile, "r") as myfile:
#         data = myfile.read()
#     cp = CIFparser.from_cifstring(data)
#     clean_cif_strings = cp.get_clean_cifs_stringlist()
#
#     c = Config.from_cifstring(clean_cif_strings[0])
#     omol = c.omols[0]
#     lbone = omol.backbone.lp
#
#     bc, bcpstructure, terminated_backbones = c.get_bone_config()
#     bc.unwrap_structure.to('cif', 'boneonlyconfig-{}.cif'.format(name))
#     pid = PackingIdentifier(bc)
#     packingd = pid.identify_heuristic()
#     return name, packingd, lbone
