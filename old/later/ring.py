import numpy as np
import warnings
import itertools
from ocelot.routines.geometry import angle_btw, Fitter
from ocelot.schema.conformation.atomlist import AtomList
from ocelot.schema.conformation.bond import Bond

class RingInitError(Exception):
    pass

class Ring(AtomList):

    def __init__(self, sites):
        """
        try not to init it outside an omol

        :var n1: first normal vector
        :var n2: second normal vector
        :var pl_error: plane fitting error
        :var ptsmean: cart geo center of fitted plane
        :var int n_member: # of sites
        :var idx: a list of siteid
        :param list sites: a list of atomsites
        """
        super().__init__(sites)
        self.checkstatus('all assigned', 'unique ids')
        self.n_member = len(self)
        if self.n_member < 3:
            raise RingInitError('you are initializing a ring with less than 3 sites!')
        normal, self.ptsmean, self.pl_error = Fitter.plane_fit(self.coordmat)
        self.n1 = normal
        self.n2 = -normal

    def as_dict(self):
        """
        keys are

        can, sites, n_member, volume, bonds, ring_insaturation
        :return:
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "sites": [s.as_dict() for s in self.atomsites],
            "n_member": self.n_member,
            "bonds": [b.as_dict() for b in self.bonds_in_ring],
        }
        return d

    @property
    def bonds_in_ring(self):
        """
        bond objects can be extracted from this ring
        e.g. for benzene the C-H bonds are NOT here

        :return: a list of bonds
        """
        bmat = self.get_bondmat(self.atomsites, self.get_distmat(self.coordmat), co=1.3)
        bonds = []
        for ij in itertools.combinations(range(len(self)), 2):
            i, j = ij
            if bmat[i][j]:
                b = Bond(self[i], self[j])
                bonds.append(b)
        return bonds


    def normal_along(self, refnormal, tol=60.0):
        """
        get the n1 or n2 that is along the direction of refnormal within a certain tol
        this is useful to identify 2 plane normals for a partially-bent structure

        :param refnormal: 3x1 array
        :param float tol: default 60 in degree
        :return: None if no normal along refnormal found
        """
        for n in [self.n1, self.n2]:
            if abs(angle_btw(n, refnormal)) < np.radians(tol):
                return n
        warnings.warn('W: no normal along this direction!')
        return None

    def iscoplane_with(self, other, tol=20.0, tolunit='degree'):
        """
        whether two rings are on the same plane with tol

        :param Ring other:
        :param tol: degree default 20
        :param tolunit: degree/radian
        :return: bool
        """
        angles = []
        if tolunit == 'degree':
            tol = np.radians(tol)
        for v1 in [self.n1, self.n2]:
            for v2 in [other.n1, other.n2]:
                angles.append(abs(angle_btw(v1, v2)))
        if min(angles) < tol:
            return True
        return False

    def __eq__(self, other):
        match = 0
        if isinstance(other, Ring):
            return 0
        if len(self) != len(other):
            return 0
        for ss in self:
            for so in other:
                if np.allclose(ss.coords, so.coords):
                    match += 1
                    break
        if match == len(self):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        outs = ['Ring:']
        for s in self:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @classmethod
    def from_idxlst(cls, idxlst, sites):
        """
        build up a ring based on a list of siteid and a list of atomsites or AtomList

        :param idxlst: a list of siteid
        :param sites: all sites (in an omol)
        """
        rs = []
        sites_idx = [s.siteid for s in sites]
        if not set(idxlst).issubset(set(sites_idx)):
            warnings.warn('idxlst is not a subset of sitelist when init ring')
            return None
        for s in sites:
            if s.siteid in idxlst:
                rs.append(s)
        return cls(rs)

