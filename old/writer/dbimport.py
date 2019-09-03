import glob
import os
from api.writer.analyser import Analyser
"""
write db entry

crystal
clean cif file, cif name, contributor, crystal family, a, b, c, cell volume, occupied ratio, z per cell,
packing pattern, density, space group, bands, effective mass

dimer (each)
slip, electronic coupling, xyz

molecule
weight, backbone, sidechain, reorganization energy, xyz, svg

backbone
length, width, linearity, # rings, density

sidechain
position, volume, density,

"""

findir = '/scratch/qai222/database/fin'

alldirs = glob.glob(findir+'*/')

for dir in alldirs:
    whereami = os.getcwd()
    os.chdir(dir)
    ana = Analyser()
    os.chdir(whereami)

