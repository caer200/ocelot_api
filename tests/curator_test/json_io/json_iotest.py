from ocelot.curator.Contribution import *

ACCESS_PROVIDER = JohnAnthony_UKY
ACCESS_DATE = datetime.date(2018, 10, 7)
ACCESS_LIC = None

dir_path = os.path.dirname(os.path.realpath(__file__))
data_access = DataAccess(ACCESS_PROVIDER, ACCESS_DATE, ACCESS_LIC)

ciffile = 'k03053.cif'
rawjsonfile = 'raw_k03053.json'
curatejsonfile = 'curate_k03053.json'

def reads(fn):
    with open(fn, 'r') as f:
        s= f.read()
    return s

import numpy as np

cifstring = reads(ciffile)
raw_data = RawData(cifstring, data_access, 'lalal', {'basename': 'k03053', 'test':np.zeros((3, 5))})
raw_data.to_jsonfile(rawjsonfile)
reproduce = RawData.from_jsonfile(rawjsonfile)
print(type(reproduce.data_properties['test']))  # <class 'numpy.ndarray'>

contri = Contribution(data_access, '{}/test'.format(dir_path))
c_data = contri.curate_one(raw_data, '{}/test'.format(dir_path))
c_data.to_jsonfile(curatejsonfile)
reproduce = CuratedData.from_jsonfile(curatejsonfile)

c: Config = reproduce.data_content['configuration']
mc: MolConformer = c.molconformers[0]
print(type(mc.backbone.pfit_vp))   # <class 'numpy.ndarray'>

