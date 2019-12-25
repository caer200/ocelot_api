from ocelot.schema.config import Config
from ocelot.task.hop import Hop, os
import numpy as np
from ocelot.routines.pbc import PBCparser, CIFparser
import json
from pprint import pprint


ZINDOLIB = '/home/ai/zindo/'
ZINDOBIN = '/home/ai/zindo/zindo'
ZINDOCTBIN = '/home/ai/zindo/zindo-ct'

# config = Config.from_file("ben.cif")
# config.unwrap_structure.to('poscar', 'POSCAR_ben')
# wdir = './ben'
#
# hop = Hop(config, zindobin=ZINDOBIN, zindoctbin=ZINDOCTBIN, zindolib=ZINDOLIB, wdir=wdir)
# mesh, hopdata, symdata, augdata, network = hop.get_hopping_network(10.0)
# pprint(network)

ciffile = 'gebw.cif'
with open(ciffile, "r") as myfile:
    data = myfile.read()
cp = CIFparser.from_cifstring(data)
clean_cif_strings = cp.get_clean_cifs_stringlist()
c = Config.from_cifstring(clean_cif_strings[0])
c.unwrap_structure.to('poscar', 'POSCAR_gebw')
wdir = './gebw'
hop = Hop(c, zindobin=ZINDOBIN, zindoctbin=ZINDOCTBIN, zindolib=ZINDOLIB, wdir=wdir)
hopdata, symdata = hop.get_hopping_network_s1()
mesh_hh, augdata_hh, network_hh = hop.get_hopping_network_s2(hopdata, 10, (1, 1, 1), 'hh')
mesh_ll, augdata_ll, network_ll = hop.get_hopping_network_s2(hopdata, 10, (1, 1, 1), 'll')
