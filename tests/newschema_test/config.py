from ocelot.schema.configuration import Config

config = Config.from_file('tipge-ss.cif')
bc, boneonly_pstructure, terminated_backbone_hmols = config.get_bone_config()
# boneonly_pstructure.to('cif', 'boneonly.cif')
# for s in boneonly_pstructure:
#     print(s, s.properties)
for s in config.pstructure:
    print(s, s.properties)

