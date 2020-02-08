# from ocelot.routines.disparser import DisParser
#
# # disorder
# ciffile = 'x17059.cif'
# dp = DisParser.from_ciffile(ciffile)
# dp.to_configs(write_files=True)  # writes conf_x.cif
#
# # backbone config
# from ocelot.schema.configuration import Config
# config = Config.from_file('conf_1.cif')
# bc, boneonly_pstructure, terminated_backbone_hmols = config.get_bone_config()
# boneonly_pstructure.to('cif', 'boneonly.cif')
#
# from ocelot.task.pkid import PackingIdentifier
# pid = PackingIdentifier(bc)
# packingd = pid.identify_heuristic()
# print(packingd[0]['packing'])

# rubrene
from rdkit.Chem import MolFromSmiles
from ocelot.schema.graph import MolGraph

smiles = 'c1ccc(cc1)c7c2ccccc2c(c3ccccc3)c8c(c4ccccc4)c5ccccc5c(c6ccccc6)c78'
rdmol = MolFromSmiles(smiles)
mg = MolGraph.from_rdmol(rdmol)
backbone, sidegroups = mg.partition_to_bone_frags('lgfr')
print(backbone)
for sg in sidegroups:
    print(sg)

from ocelot.schema.conformer import MolConformer

mc = MolConformer.from_file('rub.xyz')
bone_conformer, sccs, bg, scgs = mc.partition(coplane_cutoff=20)
print(bone_conformer)
for sc in sccs:
    print(sc)
