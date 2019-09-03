import subprocess

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

class Analyser:

    def __init__(self):
        self.data = dict()
        self.data['crystal'] = self.crystal_data()

    def crystal_data(self):
        with open('./s0/config0.cif', 'r') as f:
            clean_cif_file = f.read()
        cif_name = subprocess.run(['head', '-n1','./s0/r.cif'], stdout=subprocess.PIPE).stdout.decode('utf-8')
        cif_name = cif_name.strip().split('_')[1]
        contributor = 'John Anthony'
        crystal_family = Analyser.sg_crystal_family()

    @staticmethod
    def sg_crystal_family(sg):
        if 2 >= sg >= 1:
            return 'triclinic'
        elif 15 >= sg >= 3:
            return 'monoclinic'
        elif 16 >= sg >= 74:
            return 'Orthorhombic'
        elif 142 >= sg >= 75:
            return 'Tetragonal'
        elif 194 >= sg >= 143:
            return 'Hexagonal'
        else:
            return ''
