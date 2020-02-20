from ocelot.routines.disparser import DisParser
from ocelot.schema.conformer import BasicConformer, MolConformer
from ocelot.schema.configuration import Config

dp = DisParser.from_ciffile("AHEYUI.cif")
# disorder_class = dp.classify()
cleanconfig, occu, mols = dp.to_configs()[0]
cleanconfig = Config(cleanconfig)
molconformers = cleanconfig.molconformers
# print(molconformers[0].to_rdmol())


print(mols[0])
rdmol = MolConformer.from_pmgmol(mols[0]).to_rdmol()
