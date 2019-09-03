from api.schema.omol import OMol
from api.schema.config import Config
from api.routines.pbc import CIFparser, PBCparser
from pymatgen.core.structure import Structure
import os
import glob
import subprocess

from api.task.packing import PackingIdentifier


def dostuff(cif):
    name = cif.split('/')[-1][:-4]
    with open (cif, "r") as myfile:
        data=myfile.read()
    cp = CIFparser(data)
    cp.write_clean_cifs()

    c = Config.from_file('cleanconfig-0.cif')
    bc, bcpstructure, terminated_backbones= c.get_bone_config()
    bc.mols[0].to('xyz','x.xyz')
    # with open('x.xyz', 'r') as f:
    #     lines = f.readlines()
    # m = pybel.from_lines(lines)
    # print(m)
    # p = subprocess.Popen("babel -ixyz x.xyz -ocan", shell=True, stdout=subprocess.PIPE)
    # out = p.stdout.read()
    # return out
    # omols = [m for m in bc.omols if not m.is_solvent]
    # return [m.canonical_smiles for m in omols]

# # for cif in ['ss.cif', 'sh.cif', 'rb.cif']:
# f = open('output.txt', 'a')
# for cif in sorted(glob.glob('/home/ai/wd/oscar/cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
#     # print('working on: ' + cif)
#     name = int(cif.split('/')[-1][:-4])
#     if name >= 0:
#         o = dostuff(cif)
#         f.write(str(o) + '\n')
# f.close()




