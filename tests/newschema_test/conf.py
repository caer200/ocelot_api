from ocelot.schema.conformer import *

# # m = Molecule.from_file('tipgpn.xyz')
# # bc = BasicConformer.from_file('tipgpn.xyz')
# # print(bc.siteids)
# # bc = BoneConformer.from_file('tipgpn.xyz')
# mc = MolConformer.from_file('rub.xyz')
# # backbone_hmol, schmols = mc.to_addhmol()
# for sc in mc.sccs:
#     print(sc.conformer_properties)
# # pri
# # backbone_hmol.to('xyz', 'bone.xyz')
# # i = 0
# # for sch in schmols:
# #     sch.to('xyz', 'sc-{}.xyz'.format(i))
# #     i += 1
#
# # print(mc.siteids)
conf_a = MolConformer.from_file('a.xyz')
conf_b = MolConformer.from_file('d.xyz')
dimer = ConformerDimer(conf_a, conf_b)
dimer.plt_bone_overlap('convex')
