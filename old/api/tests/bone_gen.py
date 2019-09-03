import sys

print('Python %s on %s' % (sys.version, sys.platform))
sys.path.extend(['/home/ai/wd/oscar'])
<<<<<<< HEAD
from api.schema.omol import Mol
=======
from api.schema.Mol import Mol
>>>>>>> 88b3e96e2d37f403133c7972e5e02f735a27c716
from api.tasks.Fitbox import Boxfitter
import time
import numpy as np

ts1 = time.time()

mol = Mol.from_xyz('../sample/ge.xyz')

pbc = [
    [10, 0, 0],
    [0, 5, 10],
    [5, 5, 5]
]
bf = Boxfitter(mol, np.array(pbc))
configs = bf.gen_bone_configs(steps=9)
i=0
for c in configs:
    c.to_xyz(str(i) + '.xyz')
    i+=1