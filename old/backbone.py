import sitesop
import copy
import sys
import warnings
import numpy as np
import mathop
from scipy.spatial.distance import pdist, squareform


class Backbone:

    def __init__(self, omol, sites):
        self.sites = copy.deepcopy(sites)  # necessary?
        self.omol = omol
        self.comat = sitesop.coordmat(self.sites)
        self.distmat = squareform(pdist(self.comat))
        self.joints = [s for s in sites if s.properties['relative_position'] == 'bone-joint']
        self.centers = self.ring_centers()
        self.nrings = len(self.centers)
        self.center_line_v, self.center_line_b, self.linearity = mathop.linear_fit_3d(self.centers)
        self.vp = self.center_line_v
        self.vo = np.zeros(3)
        for i in range(len(self.sites)):
            vdummy = self.comat[i] - self.comat[i+1]
            if mathop.angle_btw(self.vp, vdummy) > 0.01:
                self.vo = np.cross(self.vp, vdummy)
                break
        self.vq = mathop.unify(np.cross(self.vo, self.vp))
        self.data = {
            'lp': self.lp,
            'lq': self.lq,
            'nrings': self.nrings,
            'linearity': self.linearity,
        }

    @property
    def lq(self):
        """
        mdiff( (s.coord - ref) proj at vq )
        :return:
        """
        ref = self.centers[0]
        projs = [np.dot(s.coords - ref, self.vq) for s in self.sites]
        return max(projs) - min(projs)

    @property
    def gom(self):
        """
        geometric center
        :return:
        """
        center = np.zeros(3)
        for s in self.sites:
            center += s.coords
        return center/len(self.sites)

    @property
    def lp(self):
        mdistance = max(self.distmat.flatten())
        mi, mj = np.where(self.distmat == mdistance)[0][0], np.where(self.distmat == mdistance)[1][0]
        v = self.comat[mi] - self.comat[mj]
        return abs(np.dot(v, self.center_line_v))

    def ring_centers(self):
        centers = []
        ss = self.sites
        nrnbrtable = sitesop.nrnbrmap(ss)
        loops4 = mathop.loop_finding(nrnbrtable, loopsize=4)
        loops5 = mathop.loop_finding(nrnbrtable, loopsize=5)
        loops6 = mathop.loop_finding(nrnbrtable, loopsize=6)
        if len(loops4 + loops5 + loops6) == 0:
            warnings.warn('no ring found!')
            sys.exit(1)
        for loop_idx_list in loops4 + loops5 + loops6:
            loop_coords = [ss[i].coords for i in loop_idx_list]
            center = np.average(loop_coords, axis=0)
            centers.append(center)
        return np.array(centers)

    def to_xyz(self, xyzname):
        sitesop.sites_toxyz(self.sites, xyzname)

    @staticmethod
    def project(self, com, p, q):
        # TODO project the backbone such that it can be defined by an arbitrary coord sys
        pass
