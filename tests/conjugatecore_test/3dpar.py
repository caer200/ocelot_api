from ocelot.schema.conformer import MolConformer
from ocelot.schema.configuration import Config
# mc = MolConformer.from_file('cnbr.xyz')
# mc.chrombone.to_file('xyz', 'chrombone.xyz')
# backbone_hmol, schmols = mc.to_addhmol('geo')
# backbone_hmol.to('xyz', 'geo-h.xyz')
# backbone_hmol, schmols = mc.to_addhmol('chrom')
# backbone_hmol.to('xyz', 'chrom-h.xyz')
#
conf = Config.from_file('AHEYUI.cif')
molconf: MolConformer = conf.molconformers[0]
print(molconf.chrombone_graph)
