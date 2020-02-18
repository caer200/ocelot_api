import numpy as np

from ocelot.routines.geometry import cart2frac
from ocelot.routines.geometry import frac2cart
"""
modified based on sasha's code @
https://github.com/afonari/emc
it looks like in there the reciprocal lattice basis was calculated following the physicist convention, so the 2pi was not
dropped in the deltas in fd, we're going to stick to 1/A (crystallographer convention)
"""

st3 = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0],
       [0.0, 0.0, 1.0], [-1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, -1.0],
       [1.0, 0.0, 1.0], [1.0, 0.0, -1.0], [-1.0, 0.0, 1.0], [0.0, -1.0, -1.0], [0.0, 1.0, 1.0], [0.0, 1.0, -1.0],
       [0.0, -1.0, 1.0], ]


def fd_effmass_st3(e, h):
    m = np.zeros((3, 3))
    m[0][0] = (e[1] - 2.0 * e[0] + e[2]) / h ** 2
    m[1][1] = (e[3] - 2.0 * e[0] + e[4]) / h ** 2
    m[2][2] = (e[5] - 2.0 * e[0] + e[6]) / h ** 2

    m[0][1] = (e[7] + e[8] - e[9] - e[10]) / (4.0 * h ** 2)
    m[0][2] = (e[11] + e[12] - e[13] - e[14]) / (4.0 * h ** 2)
    m[1][2] = (e[15] + e[16] - e[17] - e[18]) / (4.0 * h ** 2)
    # symmetrize
    m[1][0] = m[0][1]
    m[2][0] = m[0][2]
    m[2][1] = m[1][2]
    return m


st = np.array(st3)


class EmTensor:
    eigens: np.ndarray

    def __init__(self, kcoord, iband: int, stepsize: float, reci_mat: np.ndarray, eigens=None):
        self.kcoord = kcoord
        self.iband = iband
        self.stepsize = stepsize
        self.reci_mat = reci_mat
        self.eigens = eigens

    @property
    def kcart(self):
        return frac2cart(self.kcoord, self.reci_mat)

    @property
    def kmesh(self):
        h = self.stepsize
        kpoints = []
        for i in range(len(st)):
            k_c = self.kcoord + st[i] * h
            k_f = cart2frac(k_c, self.reci_mat)
            kpoints.append(k_f)
        return np.array(kpoints)

    def write_kmesh(self, fn):
        kpts = self.kmesh
        with open(fn, 'w') as f:
            f.write("EMC at {}\n".format(self.kcoord))
            f.write("%d\n" % len(st))
            f.write("Reciprocal\n")
            for kpt in kpts:
                f.write("{:.6f} {:.6f} {:.6f}\n".format(*kpt))

    @classmethod
    def from_vasprun(cls, kcoord, iband, stepsize, file, channel="up"):
        from ocelot.task.bzcal import DispersionRelationLine
        vdata = DispersionRelationLine.read_vasprun(file)
        if channel == "down":
            eigens = vdata['eigens_down']
            if eigens is None:
                raise KeyError('Spin down eigens not found!')
        else:
            eigens = vdata['eigens_up']
        flat_eigens = [x[iband][0] for x in eigens]
        flat_eigens = np.array(flat_eigens)
        reci_mat = vdata['reci_mat']
        return cls(kcoord, iband, stepsize, reci_mat, flat_eigens)

    def cal_emtensor(self):
        fdm = fd_effmass_st3(self.eigens, self.stepsize)
        e, v = np.linalg.eig(fdm)

        eigenvs_cart = np.zeros((3, 3))
        eigenvs_frac = np.zeros((3, 3))
        es = np.zeros(3)
        ems = np.zeros(3)
        for i in range(3):
            eigenvs_cart[i] = v[:, i]
            eigenvs_frac[i] = cart2frac(v[:, i], self.reci_mat)
            es[i] = e[i]
            ems[i] = 1 / e[i]

        return ems, es, eigenvs_frac, eigenvs_cart

def get_reci_mat(file, filetype='poscar'):
    if filetype == 'poscar':
        from pymatgen.io.vasp.inputs import Poscar
        poscar = Poscar.from_file(file, check_for_POTCAR=False)
        reci_mat = poscar.structure.lattice.reciprocal_lattice_crystallographic.matrix
    elif filetype == 'vasprun':
        from ocelot.task.bzcal import DispersionRelationLine
        vdata = DispersionRelationLine.read_vasprun(file)
        reci_mat = vdata['reci_mat']
    else:
        raise NotImplementedError('filetype {} not implemented'.format(filetype))
    return reci_mat


def rb2ra2pi(rb:float):
    return rb/0.529177/2/np.pi

def rb2ra(rb:float):
    return rb/0.529177

if __name__ == '__main__':
    # emt = EmTensor.from_vasprun((0, 0, 0), 118, rb2ra2pi(0.035), 'vasprun.xml')
    emt = EmTensor.from_vasprun((0, 0, 0), 118, 0.035, 'vasprun.xml')
    # print(emt.reci_mat*2*np.pi)
    # print(emt.kmesh)
    ems, es, eigenvs_frac, eigenvs_cart = emt.cal_emtensor()
    print(ems)
    # reci_mat = get_reci_mat('POSCAR')
    # emt = EmTensor((0, 0, 0), 119, rb2ra2pi(0.035), reci_mat)

