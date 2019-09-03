import sys

import numpy as np

from routines.Geometry import angle_btw, Fitter
from schema.Bond import Bond
from schema.Sitelist import Sitelist


class Ring(Sitelist):

    def __init__(self, sites, mol):
        super().__init__(sites)
        self.mol = mol
        self.number = len(self.sites)
        self.idx = [s.id for s in self.sites]
        if self.number < 3:
            sys.exit('you are initializing a ring with less than 3 msites!')

    @property
    def normals(self):
        coordmat = self.mol.coordmat[[s.id for s in self.sites]]
        normal, ptsmean, pl_error = Fitter.plane_fit(coordmat)
        return np.array([normal, -normal])

    @property
    def bonds(self):
        bonds = []
        for si in self.idx:
            for j in self.mol.nbrmap[si]:
                b = Bond([self.mol.sites[si], self.mol.sites[j]])
                if b not in bonds:
                    bonds.append(b)
        return bonds

    def normal_along(self, refnormal, tol=60.0):
        for n in self.normals:
            if abs(angle_btw(n, refnormal)) < np.radians(tol):
                return n
        sys.exit('no normal along this direction!')

    def iscoplane_with(self, other, tol=15.0):
        angles = []
        for v1 in self.normals:
            for v2 in other.normals:
                angles.append(abs(angle_btw(v1, v2)))
        if min(angles) < np.pi * tol / 180.0:
            return 1
        return 0

    def isfused_with(self, other):
        if len(self.interscet(other)) > 1:
            return 1
        return 0

    def interscet(self, other):
        # only used after init id
        if isinstance(other, Ring):
            return list(set(self.idx) & set(other.idx))
        else:
            sys.exit('intersceting ring with another non-ring obj!')

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __eq__(self, other):
        match = 0
        if len(self.sites) != len(other.sites):
            return 0
        for ss in self.sites:
            for so in other.sites:
                if np.allclose(ss.coords, so.coords):
                    match += 1
                    break
        if match == len(self.sites):
            return 1
        else:
            return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        outs = ['Ring:']
        for s in self.sites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @classmethod
    def from_idxlst(cls, idxlst, mol):
        rs = []
        for s in mol.sites:
            if s.id in idxlst:
                rs.append(s)
        return cls(rs, mol)

    def insaturation(self):
        # TODO get unsaturation
        pass
