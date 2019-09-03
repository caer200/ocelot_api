import itertools
import math

import numpy as np
from pymatgen.core.structure import Molecule
from api.routines.geometry import norm, rotate_along_axis
from api.schema.msite import MSite


class MSitelist:

    def __init__(self, sites, description=''):
        """
        this is the parent obj for Mol, Ring, Bond, Sidechain, Backbone
        :param sites: a *list* of msites
        :param description:
        """
        self.sites = sites
        self.description = description

    def __len__(self):
        return len(self.sites)

    def __contains__(self, site):
        return site in self.sites

    def __iter__(self):
        return self.sites.__iter__()

    def __getitem__(self, ind):
        return self.sites[ind]

    def __repr__(self):
        outs = ['MSitelist:']
        for s in self.sites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    def get_coordmat(self):
        """
        coordinates matrix
        """
        coordmat = np.empty((len(self.sites), 3))
        for i in range(len(self.sites)):
            for j in range(3):
                coordmat[i][j] = self.sites[i].coords[j]
        return coordmat

    @property
    def geoc(self):
        """
        geometric center
        :return:
        """
        v = np.zeros(3)
        for s in self.sites:
            v += s.coords
        return v / len(self.sites)

    def intersection(self, other):
        """
        :param other:
        :return: a list of msites in self that belong to both Sitelists, there is no copy!
        """
        r = []
        for s in self.sites:
            if s in other:
                r.append(s)
        return r

    def get_site(self, siteid):
        for s in self.sites:
            if s.siteid == siteid:
                return s
        return None

    def issubset(self, other):
        return len(self.sites) == len(self.intersection(other))

    @classmethod
    def from_coordmat(cls, mat, names, ids=None):
        """

        :param mat:
        :param names: a list of string as element_name
        :param ids:
        :return:
        """
        if ids is None:
            ids = np.ones(len(names))
            ids[:] = -1
        ss = []
        for i in range(len(names)):
            ss.append(MSite(names[i], mat[i], ids[i]))
        return cls(ss)

    @classmethod
    def from_file(cls, fname):
        m = Molecule.from_file(fname)
        return cls([MSite.from_pymatgen_site(s) for s in m.sites])

    def to(self, ftype, fname):
        pymatgen_sites = [s.to_pymatgen_site() for s in self.sites]
        m = Molecule.from_sites(pymatgen_sites)
        m.to(ftype, fname)

    def rotate(self, theta, axis, origin, unit='degree'):
        """
        rotate the vectors defined the origin along given axis
        notice the coords are changed in-situ
        :param theta:
        :param axis:
        :param origin:
        :param unit:
        :return:
        """
        for s in self.sites:
            v = s.coords - origin
            s.coords = origin + rotate_along_axis(v, axis, theta, thetaunit=unit)

    def rotate_with_matrix(self, matrix, origin):
        """
        a quicker version rotate if we konw the rotation matrix
        :param matrix:
        :param origin:
        :return:
        """
        for s in self.sites:
            v = s.coords - origin
            s.coords = origin + np.dot(matrix, v)

    def volume(self, boxdensity=0.5):
        """
        # TODO profile
        http://wiki.bkslab.org/index.php/Calculate_volume_of_the_binding_site_and_molecules
        First, Lay a grid over the spheres.
        Count the number or points contained in the spheres (Ns).
        Count the number of points in the grid box (Ng).
        Calculate the volume of the grid box (Vb).
        :return:
        """
        mat = np.empty((len(self.sites), 4))
        for i in range(len(self.sites)):
            for j in range(3):
                mat[i][j] = self.sites[i].coords[j]
            mat[i][3] = self.sites[i].element.covrad

        box_min = math.floor(min(itertools.chain.from_iterable(mat))) - 2
        box_max = math.ceil(max(itertools.chain.from_iterable(mat))) + 2
        axis = np.arange(box_min, box_max + boxdensity, boxdensity)
        grid = np.array(np.meshgrid(axis, axis, axis)).T.reshape(-1, 3)
        ngps_total = len(grid)
        ngps_occu = 0
        for igp in range(len(grid)):
            for iap in range(len(mat)):
                dist = norm(grid[igp] - mat[iap][:3])
                if dist < mat[iap][3]:
                    ngps_occu += 1
                    break
        v = (ngps_occu / ngps_total) * ((box_max - box_min) ** 3)
        return v
