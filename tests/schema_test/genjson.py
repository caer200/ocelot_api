from ocelot.routines.pbc import CIFparser
from ocelot.schema.config import Config
from ocelot.task.pkid import PackingIdentifier
import glob
import json
import numpy as np


# def get_tips():
#     with open ("../cifs/128.cif", "r") as myfile:
#         data=myfile.read()
#     cp = CIFparser.from_cifstring(data)
#     cp.write_clean_cifs(prefix='tips')
#     c = Config.from_file('{}-0.cif'.format('tips'))
#     m = c.omols[0]
#     tips_sc = sorted(m.sidechains, key=lambda x: len(x), reverse=True)[0]
#     return tips_sc


def parser(cif, verbose=0):
    name = cif.split('/')[-1][:-4]
    with open(cif, "r") as myfile:
        data = myfile.read()
    cp = CIFparser.from_cifstring(data)
    cp.write_clean_cifs(prefix=name)
    c = Config.from_file('{}-0.cif'.format(name))
    c.unwrap_structure.to('cif', 'uconfig-{}.cif'.format(name))
    bc, bcpstructure, terminated_backbones = c.get_bone_config()
    bc.unwrap_structure.to('cif', 'bconfig-{}.cif'.format(name))
    pid = PackingIdentifier(bc)
    packingdata = pid.identify_heuristic()
    m = [m for m in c.omols if not m.is_solvent][0]
    d = dict()
    if verbose:
        d['config'] = c.as_dict()
        d['cifname'] = name
        d['packingdata'] = packingdata
        return d
    d['bulky_scs'] = [sc for sc in m.as_dict()['sidechains'] if not len(sc['msites']) > 3]
    d['cifname'] = name
    d['packingdata'] = packingdata
    return d


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


tipcif = '../cifs/128.cif'

for cif in sorted(glob.glob('../cifs/*.cif'), key=lambda x: int(x.split('/')[-1][:-4])):
    # if int(cif.split('/')[-1][:-4]) >= 196:
    print('working on {}'.format(cif))
    # d = parser(tipcif, verbose=1)
    # with open('./full.json'.format(d['cifname']), 'w') as fp:
    #     json.dump(d, fp, cls=NumpyEncoder)
    d = parser(cif)
    with open('./jsons/{}.json'.format(d['cifname']), 'w') as fp:
        json.dump(d, fp, cls=NumpyEncoder)
