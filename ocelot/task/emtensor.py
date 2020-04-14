import numpy as np

from ocelot.routines.geometry import cart2frac
from ocelot.routines.geometry import frac2cart
from ocelot.routines.mathop import ev_to_ha

"""
modified based on sasha's code @
https://github.com/afonari/emc
it looks like in there the reciprocal lattice basis was calculated following the physicist convention, so the 2pi was not
dropped in the deltas in fd, we're going to stick to 1/A (crystallographer convention)

choose stepsize:
if we know the lowest line effective mass was obtained at a point along GX, this point can be G or X
calculate the length of G-X say it is l_gx
the initial stepsize should be l_gx/4 (for st3) or l_gx/8 (for st5)
then scan downwards till converge, stepsize for scan can be 10% of init stepsize, e.g. 0.04, 0.036, 0.032...
"""

st3 = [[0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0],
       [0.0, 0.0, 1.0], [-1.0, -1.0, 0.0], [1.0, 1.0, 0.0], [1.0, -1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, 0.0, -1.0],
       [1.0, 0.0, 1.0], [1.0, 0.0, -1.0], [-1.0, 0.0, 1.0], [0.0, -1.0, -1.0], [0.0, 1.0, 1.0], [0.0, 1.0, -1.0],
       [0.0, -1.0, 1.0], ]

st5 = [[0.0, 0.0, 0.0], [-2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, -2.0, 0.0],
       [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, -2.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0],
       [0.0, 0.0, 2.0], [-2.0, -2.0, 0.0], [-1.0, -2.0, 0.0], [1.0, -2.0, 0.0], [2.0, -2.0, 0.0], [-2.0, -1.0, 0.0],
       [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [2.0, -1.0, 0.0], [-2.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [1.0, 1.0, 0.0],
       [2.0, 1.0, 0.0], [-2.0, 2.0, 0.0], [-1.0, 2.0, 0.0], [1.0, 2.0, 0.0], [2.0, 2.0, 0.0], [-2.0, 0.0, -2.0],
       [-1.0, 0.0, -2.0], [1.0, 0.0, -2.0], [2.0, 0.0, -2.0], [-2.0, 0.0, -1.0], [-1.0, 0.0, -1.0], [1.0, 0.0, -1.0],
       [2.0, 0.0, -1.0], [-2.0, 0.0, 1.0], [-1.0, 0.0, 1.0], [1.0, 0.0, 1.0], [2.0, 0.0, 1.0], [-2.0, 0.0, 2.0],
       [-1.0, 0.0, 2.0], [1.0, 0.0, 2.0], [2.0, 0.0, 2.0], [0.0, -2.0, -2.0], [0.0, -1.0, -2.0], [0.0, 1.0, -2.0],
       [0.0, 2.0, -2.0], [0.0, -2.0, -1.0], [0.0, -1.0, -1.0], [0.0, 1.0, -1.0], [0.0, 2.0, -1.0], [0.0, -2.0, 1.0],
       [0.0, -1.0, 1.0], [0.0, 1.0, 1.0], [0.0, 2.0, 1.0], [0.0, -2.0, 2.0], [0.0, -1.0, 2.0], [0.0, 1.0, 2.0],
       [0.0, 2.0, 2.0]]


def fd_effmass_st5(e, h):
    m = np.zeros((3, 3))
    m[0][0] = (-(e[1] + e[4]) + 16.0 * (e[2] + e[3]) - 30.0 * e[0]) / (12.0 * h ** 2)
    m[1][1] = (-(e[5] + e[8]) + 16.0 * (e[6] + e[7]) - 30.0 * e[0]) / (12.0 * h ** 2)
    m[2][2] = (-(e[9] + e[12]) + 16.0 * (e[10] + e[11]) - 30.0 * e[0]) / (12.0 * h ** 2)
    #
    m[0][1] = (-63.0 * (e[15] + e[20] + e[21] + e[26]) + 63.0 * (e[14] + e[17] + e[27] + e[24]) + 44.0 * (
            e[16] + e[25] - e[13] - e[28]) + 74.0 * (e[18] + e[23] - e[19] - e[22])) / (600.0 * h ** 2)
    m[0][2] = (-63.0 * (e[31] + e[36] + e[37] + e[42]) + 63.0 * (e[30] + e[33] + e[43] + e[40]) + 44.0 * (
            e[32] + e[41] - e[29] - e[44]) + 74.0 * (e[34] + e[39] - e[35] - e[38])) / (600.0 * h ** 2)
    m[1][2] = (-63.0 * (e[47] + e[52] + e[53] + e[58]) + 63.0 * (e[46] + e[49] + e[59] + e[56]) + 44.0 * (
            e[48] + e[57] - e[45] - e[60]) + 74.0 * (e[50] + e[55] - e[51] - e[54])) / (600.0 * h ** 2)
    # symmetrize
    m[1][0] = m[0][1]
    m[2][0] = m[0][2]
    m[2][1] = m[1][2]
    return m


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


# st = np.array(st3)
# st = np.array(st5)


class EmTensor:
    eigens: np.ndarray

    def __init__(self, kcoord, iband: int, stepsize: float, reci_mat: np.ndarray, eigens=None, st=3):
        """
        init a effective tensor calculation

        :param kcoord: in frac, center of the stencil
        :param iband: starts from 0
        :param stepsize: in 1/A
        :param reci_mat: in 1/A
        :param eigens: in eV
        :param int st: size of the stencil, 3 or 5
        """
        self.kcoord = kcoord
        self.iband = iband
        self.stepsize = stepsize
        self.reci_mat = reci_mat
        self.eigens = eigens
        self.real_mat = np.linalg.inv(self.reci_mat).T
        if st == 3:
            self.st = np.array(st3)
        else:
            self.st = np.array(st5)

    @property
    def kcart(self):
        return frac2cart(self.kcoord, self.reci_mat)

    @property
    def kmesh(self):
        h = self.stepsize  # in 1/A
        kpoints = []
        for i in range(len(self.st)):
            k_c = self.kcart + self.st[i] * h
            k_f = cart2frac(k_c, self.reci_mat)
            kpoints.append(k_f)
        return np.array(kpoints)

    def write_kmesh(self, fn):
        kpts = self.kmesh
        with open(fn, 'w') as f:
            f.write("EMC at frac {}\n".format(self.kcoord))
            f.write("%d\n" % len(self.st))
            f.write("Reciprocal\n")
            for kpt in kpts:
                f.write("{:.6f} {:.6f} {:.6f}  1.01\n".format(*kpt))

    @classmethod
    def from_vasprun(cls, kcoord, iband, stepsize, file, channel="up", st=3):
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
        return cls(kcoord, iband, stepsize, reci_mat, flat_eigens, st=st)

    def cal_emtensor(self):
        eigens = ev_to_ha(self.eigens)
        if len(self.st) == 19:
            fdm = fd_effmass_st3(eigens, ra2rb(self.stepsize))
        else:
            fdm = fd_effmass_st5(eigens, ra2rb(self.stepsize))
        e, v = np.linalg.eig(fdm)

        eigenvs_cart = np.zeros((3, 3))
        eigenvs_frac = np.zeros((3, 3))
        es = np.zeros(3)
        ems = np.zeros(3)

        for i in range(3):
            eigenvs_cart[i] = v[:, i]
            eigenvs_frac[i] = cart2frac(v[:, i], self.real_mat)
            eigenvs_frac[i] = eigenvs_frac[i] / max(np.abs(eigenvs_frac[i]))
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

def ra2rb(ra:float):
    return ra * 0.529177

def rb2ra2pi(rb:float):
    return rb/0.529177/2/np.pi

def rb2ra(rb:float):
    return rb/0.529177


if __name__ == '__main__':
    emt = EmTensor.from_vasprun((0, 0, 0), 118, 0.035, 'vasprun.xml', st=3)
    print(len(emt.kmesh))
    ems, es, eigenvs_frac, eigenvs_cart = emt.cal_emtensor()
    from pprint import pprint

    pprint(ems)
    pprint(eigenvs_cart)
    pprint(eigenvs_frac)
