from ocelot.schema.configuration import Config
from ocelot.task.hop import Hop, os
import numpy as np
from ocelot.routines.pbc import PBCparser
from ocelot.routines.disparser import DisParser
import json
from pprint import pprint


ZINDOLIB = '/home/ai/zindo/'
ZINDOBIN = '/home/ai/zindo/zindo'
ZINDOCTBIN = '/home/ai/zindo/zindo-ct'

def test(ciffile, wdir):
    dp = DisParser.from_ciffile(ciffile)
    c = dp.to_configs()[0][0]
    c = Config(c)
    c.unwrap_structure.to('poscar', 'POSCAR_gebw')
    os.system('mkdir -p {}'.format(wdir))
    hop = Hop(c, zindobin=ZINDOBIN, zindoctbin=ZINDOCTBIN, zindolib=ZINDOLIB, wdir=wdir)
    hopdata, symdata = hop.get_hopping_network_s1()
    mesh_hh, augdata_hh, network_hh = hop.get_hopping_network_s2(hopdata, 10, (1, 1, 1), 'hh')
    pprint(network_hh)
    mesh_ll, augdata_ll, network_ll = hop.get_hopping_network_s2(hopdata, 10, (1, 1, 1), 'll')
    pprint(network_ll)

test('ben.cif', './ben')
test('gebw.cif', './gebw')
