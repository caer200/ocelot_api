"""
config init from clean (disorder-free) structure
"""
import glob
import os
from ocelot.schema.configuration import Config
from ocelot.schema.conformer import DimerCollection, ConformerDimer

def testa():
    cifs = glob.glob('*.cif')
    for ciffile in cifs:
        name = ciffile.split('/')[-1][:-4]
        config = Config.from_file(ciffile)
        # mcs = config.molconformers
        bc, boneonly_pstructure, terminated_backbone_hmols = config.get_bone_config()
        boneonly_pstructure.to('cif', '{}_boneonly.cif'.format(name))
        # for s in boneonly_pstructure:
        #     print(s, s.properties)

def testb():
    cifs = glob.glob('*.cif')
    for ciffile in cifs:
        name = ciffile.split('/')[-1][:-4]
        config = Config.from_file(ciffile)
        # mcs = config.molconformers
        dimers, transv_fcs = config.get_dimers_array(fast=True, symm=False)
        for i in range(len(dimers)):
            dc = DimerCollection(dimers[i].flatten())
            dc.to_xyz('{}_dimercoll_{}.xyz'.format(name, i), lalabel=True)
        dimers[0][0][1].to_xyz('{}_dimer001.xyz'.format(name))
        dimers[0][0][1].plt_bone_overlap('convex', output='{}_001convex.eps'.format(name))
        dimers[0][0][1].plt_bone_overlap('concave', output='{}_001concave.eps'.format(name))

os.chdir('./iconfs')
testa()
# testb()
os.chdir('../')


