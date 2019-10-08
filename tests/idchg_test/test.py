from ocelot.task.idchg import IdChg
from pymatgen.core.structure import Molecule


cmd = 'singularity exec /home/ai/bin/mopac.sif /opt/mopac/MOPAC2016.exe'

xyzfn = 'bf4.xyz'


def idchg(xyz):
    name = xyz.split('.')[0]
    mol = Molecule.from_file(xyz)
    idchg = IdChg(mol, name, cmd + ' ' + name)
    idchg.write_inp('./')
    return idchg.run('./')

d = idchg(xyzfn)
print(d)

