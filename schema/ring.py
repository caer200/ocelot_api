import numpy as np
import sys
import warnings
from api.routines.geometry import angle_btw, Fitter
from api.schema.msitelist import MSitelist
from api.schema.bond import Bond


class Ring(MSitelist):

    def __init__(self, msites):
        for ms in msites:
            if ms.siteid == -1:
                warnings.warn('you are init a ring with sites not in an omol obj')
        super().__init__(msites)
        self.n_member = len(self)
        self.idx = [s.siteid for s in self.msites]
        if self.n_member < 3:
            sys.exit('you are initializing a ring with less than 3 msites!')
        normal, self.ptsmean, self.pl_error = Fitter.plane_fit(self.coordmat)
        self.n1 = normal
        self.n2 = -normal

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "can": self.canonical_smiles,
            "msites": [s.as_dict() for s in self.msites],
            "n_member": self.n_member,
            "volume": self.volume,
            "bonds": [b.as_dict() for b in self.bonds],
            "ring_insaturation": self.ring_insaturation,
        }
        return d

    @property
    def bonds(self):
        bonds = []
        for si in self.idx:
            sa = self.get_site_byid(si)
            for sb in sa.nbs:
                b = Bond(sa, sb)
                if b not in bonds:
                    bonds.append(b)
        return bonds

    def normal_along(self, refnormal, tol=60.0):
        """
        get the n1 or n2 that is along the direction of refnormal within a certain tol
        :param refnormal:
        :param tol:
        :return:
        """
        for n in [self.n1, self.n2]:
            if abs(angle_btw(n, refnormal)) < np.radians(tol):
                return n
        warnings.warn('no normal along this direction!')
        return None

    def iscoplane_with(self, other, tol=20.0):
        angles = []
        for v1 in [self.n1, self.n2]:
            for v2 in [other.n1, other.n2]:
                angles.append(abs(angle_btw(v1, v2)))
        if min(angles) < np.radians(tol):
            return True
        return False

    def isfused_with(self, other):
        if len(self.interscet(other)) > 1:
            return 1
        return 0

    def isconnected_with(self, other):
        nbs_idx_of_ring = []
        for s in self.msites:
            nbs_idx_of_ring += s.nbs_idx
        nbs_idx_of_other = []
        for s in other.msites:
            nbs_idx_of_other += s.nbs_idx
        intersection = set(nbs_idx_of_ring).intersection(set(nbs_idx_of_other))
        if len(intersection) > 0:
            return True
        return False

    def interscet(self, other):
        # only used after init id
        if isinstance(other, Ring):
            return list(set(self.idx) & set(other.idx))
        else:
            warnings.warn('intersceting ring with another non-ring obj!')
            return None

    def __eq__(self, other):
        match = 0
        if len(self.msites) != len(other.sites):
            return 0
        for ss in self.msites:
            for so in other.msites:
                if np.allclose(ss.coords, so.coords):
                    match += 1
                    break
        if match == len(self.msites):
            return True
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        outs = ['Ring:']
        for s in self.msites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @classmethod
    def from_idxlst(cls, idxlst, msites):
        rs = []
        msites_idx = [s.siteid for s in msites]
        if not set(idxlst).issubset(set(msites_idx)):
            warnings.warn('idxlst is not a subset of sitelist when init ring')
            return None
        for s in msites:
            if s.siteid in idxlst:
                rs.append(s)
        return cls(rs)

    @property
    def ring_insaturation(self):
        insaturated_sites = 0
        for s in self.msites:
            if s.insaturation > 0:
                insaturated_sites += 1
        return insaturated_sites / self.n_member
