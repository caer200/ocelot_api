from ocelot.schema.conformer import MolConformer

mc = MolConformer.from_file('cnbr.xyz')
mc.chrombone.to_file('xyz', 'chrombone.xyz')
backbone_hmol, schmols = mc.to_addhmol('geo')
backbone_hmol.to('xyz', 'geo-h.xyz')
backbone_hmol, schmols = mc.to_addhmol('chrom')
backbone_hmol.to('xyz', 'chrom-h.xyz')

