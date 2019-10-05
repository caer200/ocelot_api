from ocelot.schema import Sitelist
from ocelot.routines import Fitter, angle_btw, unify, coord_transform
from copy import deepcopy
import numpy as np


class Backbone(Sitelist):

    def __init__(self, sites, mol):
        super().__init__(sites)
        self.mol = mol
        self.rings = mol.largest_fr()
        self.nrings = len(self.rings)
        self.vp, self.center_line_b, self.linearity = Fitter.linear_fit([r.geoc for r in self.rings])
        self.vp = unify(self.vp)
        self.vo = np.zeros(3)
        for i in range(len(self.sites)-1):
            vdummy = self.sites[i].coords - self.sites[i+1].coords
            if angle_btw(self.vp, vdummy) > 0.01:
                self.vo = unify(np.cross(self.vp, vdummy))
        self.vq = unify(np.cross(self.vo, self.vp))
        self.pqomat = np.array([self.vp, self.vq, self.vo])

    @property
    def data(self):
        d = {}
        d['lq'] = self.lq
        d['lp'] = self.lp
        d['linearity'] = self.linearity
        d['nrings'] = self.nrings
        return d

    @property
    def lq(self):
        """
        mdiff( (s.coord - ref) proj at vq )
        :return:
        """
        ref = self.rings[0].geoc
        projs = [np.dot(s.coords - ref, self.vq) for s in self.sites]
        return max(projs) - min(projs)

    @property
    def lp(self):
        ref = self.rings[0].geoc
        projs = [np.dot(s.coords - ref, self.vp) for s in self.sites]
        return max(projs) - min(projs)

    @classmethod
    def from_mol(cls, mol):
        bsites = []
        for r in mol.largest_fr():
            for s in r.sites:
                if s not in bsites:
                    bsites.append(s)
        return cls(bsites, mol)

    @classmethod
    def orient(cls, obb):
        nss = []
        for s in obb.sites:
            ns = deepcopy(s)
            ns.coords = coord_transform(obb.vp, obb.vq, obb.vo, s.coords-obb.geoc)
            nss.append(ns)
        return cls(nss, obb.mol)





