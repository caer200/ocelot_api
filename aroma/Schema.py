import Constants
import Routines
import warnings
import copy
import numpy as np
import sys
from scipy.spatial.distance import pdist, squareform


class Element:

    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        if self.name == other.name:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.name)

    @property
    def covrad(self):
        try:
            return Constants.COVRAD[self.name]
        except KeyError:
            return None

    @property
    def vanrad(self):
        try:
            return Constants.VANRAD[self.name]
        except KeyError:
            return None

    @property
    def atomic_number(self):
        try:
            return Constants.ATOMN[self.name]
        except KeyError:
            return None

    @property
    def valence(self):
        try:
            return Constants.VALENCE[self.name]
        except KeyError:
            return None


    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

class Site:
    def __init__(self, element_name, coords, siteid=-1):
        self.element = Element(element_name)
        self._id = int(siteid)
        self._coords = np.empty(3)
        for i in range(3):
            self._coords[i] = float(coords[i])

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    @property
    def id(self):
        return self._id

    @id.setter
    def id(self, v):
        try:
            self._id = int(v)
        except ValueError:
            warnings.warn('id must be int, but you tried to set it as:', v)

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, vs):
        for i in range(3):
            self._coords[i] = float(vs[i])

    def distance(self, other):
        # TODO cythonize
        return np.linalg.norm(self.coords - other.coords)

    def __eq__(self, other):
        if self.element == other.element and np.allclose(self.coords, other.coords) and self.id == other.id:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.element)

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f}) id: {}".format(
            self.element, *self.coords, self.id)

    def __str__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f}) id: {}".format(
            self.element, *self.coords, self.id)


class SiteList:

    def __init__(self, sites):
        self.sites = copy.deepcopy(sites)

    @property
    def coordmat(self):
        """
        coordinates matrix
        """
        mat = np.empty((len(self.sites), 3))
        for i in range(len(self.sites)):
            for j in range(3):
                mat[i][j] = self.sites[i].coords[j]
        return mat

    @property
    def distmat(self):
        return squareform(pdist(self.coordmat))

    @property
    def nbrmap(self):
        # TODO profile
        bmat = self.bondmat
        ma = {}
        numsites = len(bmat)
        for i in range(numsites):
            nbs = []
            for j in range(numsites):
                if bmat[i][j] and j != i:
                    nbs.append(j)
            ma[i] = nbs
        return ma

    @property
    def bondmat(self, co=1.3):
        """
        Bij = whether there is a bond between si and sj
        i is NOT bonded with itself
        :param co: coefficient for cutoff
        :return:
        """
        numsites = len(self.sites)
        mat = np.ones((numsites, numsites), dtype=bool)
        for i in range(numsites):
            mat[i][i] = False
            for j in range(i + 1, numsites):
                if self.distmat[i][j] > (self.sites[i].element.covrad + self.sites[j].element.covrad) * co:
                    mat[i][j] = False
                    mat[j][i] = False
        return mat

    def to_xyz(self, fname):
        with open(fname, 'w') as f:
            f.write(str(len(self.sites)) + '\r\n\r\n')
            for s in self.sites:
                f.write(s.element.name + ' ' + ' '.join([str(c) for c in s.coords]) + '\r\n')

    def __repr__(self):
        outs = ['SiteList:']
        for s in self.sites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

class Bond(SiteList):
    def __init__(self, sites):
        super().__init__(sites)
        if len(self.sites) != 2:
            sys.exit('you are initializing a bond with <2 or >2 sites')
        self.center = (self.sites[0].coords + self.sites[1].coords)*0.5
        self.length = np.linalg.norm(self.sites[0].coords - self.sites[1].coords)
        self.type = {self.sites[0].element, self.sites[1].element}
        # TODO add insaturation

    def __eq__(self, other):
        if self.type == other.type and np.allclose(self.center, other.center):
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "Bond: {} ({:.4f}, {:.4f}, {:.4f}) length: {}".format(
            self.type, *self.center, self.length)

    def __str__(self):
        return "Bond: {} ({:.4f}, {:.4f}, {:.4f}) length: {}".format(
            self.type, *self.center, self.length)

class Ring(SiteList):

    def __init__(self, sites):
        super().__init__(sites)
        self.number = len(self.sites)
        if self.number < 3:
            sys.exit('you are initializing a ring with less than 3 sites!')
        normal, self.ptsmean, self.pl_error = Routines.GeoFitter3d.plane_fit(self.coordmat)
        self.normals = np.array([normal, -normal])
        self.idx = [s.id for s in self.sites]

    def normal_along(self, refnormal, tol=60.0):
        for n in self.normals:
            if abs(Routines.Geo.angle_btw(n, refnormal)) < np.radians(60.0):
                return n
        sys.exit('no normal along this direction!')

    @property
    def bonds(self):
        bonds = []
        for si in range(self.number):
            for j in self.nbrmap[si]:
                b = Bond([self.sites[si], self.sites[j]])
                if b not in bonds:
                    bonds.append(b)
        return bonds

    @property
    def geoc(self):
        """
        geometric center
        :return:
        """
        v = np.zeros(3)
        for s in self.sites:
            v += s.coords
        return v/len(self.sites)

    def iscoplane_with(self, other, tol=15.0):

        angles = []
        for v1 in self.normals:
            for v2 in other.normals:
                angles.append(abs(Routines.Geo.angle_btw(v1, v2)))
        if min(angles) < np.pi*tol/180.0:
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
    def from_idxlst(cls, idxlst, sites):
        rs = []
        for s in sites:
            if s.id in idxlst:
              rs.append(s)
        return cls(rs)

    def unsaturation(self):
        # TODO get unsaturation
        pass

class Mol(SiteList):

    def __init__(self, sites, name=None, charge=0, multiplicity=1, comment=''):
        super().__init__(sites)
        for i in range(len(self.sites)):
            self.sites[i].id = i
            self.sites[i].nbs = [self.sites[j] for j in self.nbrmap[i]]  # will the elements in nbs also change?
            self.sites[i].insaturation = abs(len(self.sites[i].nbs) - self.sites[i].element.valence)
        self.name = name
        self.comment = comment
        self.charge = charge
        self.multiplicity = multiplicity
        rs = Routines.LoopSearcher(self.nbrmap)
        self.rings = []
        for ring_size in range(3, 9):
            idxlsts = rs.alex_method(ring_size)
            self.rings += [Ring.from_idxlst(idxlst, self.sites) for idxlst in idxlsts]

        # TODO add parser for tips-like mol

    def nics_sigma_structure(self, normaldirection=0):

        lgfr = self.largest_fr()
        geocs = [r.geoc for r in lgfr]  # nx3
        dmat = squareform(pdist(np.array(geocs)))
        maxi, maxj = np.unravel_index(dmat.argmax(), dmat.shape)
        fr = sorted(lgfr, key=lambda r: np.linalg.norm(r.geoc - lgfr[maxi].geoc))
        available_normals = fr[0].normals
        normals = [r.normal_along(available_normals[normaldirection]) for r in fr]

        frsites = []
        for r in lgfr:
            for s in r.sites:
                if s not in frsites:
                    frsites.append(s)

        sigma_sites = copy.deepcopy(self.sites)

        terminated_sites_id = []
        added_hydrogens = []

        for i in range(len(lgfr)):
            for s in lgfr[i].sites:
                if s.insaturation == 1 and s.id not in terminated_sites_id:
                    terminated_sites_id.append(s.id)
                    hs = Site('H', 1.0*normals[i]+s.coords)
                    added_hydrogens.append(hs)
                for nb in s.nbs:
                    if nb.id not in terminated_sites_id and nb.insaturation == 1 and nb not in frsites and all([nb not in ring.sites for ring in self.rings]) :
                        terminated_sites_id.append(nb.id)
                        hs = Site('H', 1.0*normals[i]+nb.coords)
                        added_hydrogens.append(hs)
        return SiteList(sigma_sites + added_hydrogens)

    def nics_line_scan_path(self, step_size, nrings, height=1.7, normaldirection=0):
        lgfr = self.largest_fr()[:nrings]
        geocs = [r.geoc for r in lgfr]  # nx3
        dmat = squareform(pdist(np.array(geocs)))
        maxi, maxj = np.unravel_index(dmat.argmax(), dmat.shape)
        fr = sorted(lgfr, key=lambda r: np.linalg.norm(r.geoc - lgfr[maxi].geoc))
        available_normals = fr[0].normals
        normals = [r.normal_along(available_normals[normaldirection]) for r in fr]

        pts = []
        ring_idx = []
        xnumbers = []
        xticks = [0]
        cross_point = np.zeros(3)
        for i in range(len(fr)-1):
            ring = fr[i]
            nb_ring = fr[i+1]
            segend1 = ring.geoc
            segend2 = nb_ring.geoc
            # nb_ring center -- cross_pt -- ring center -- prev_cross_pt -- prev ring center
            if i == 0:
                bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: Routines.Geo.angle_btw(segend2-segend1, bc-segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
                subpath = Routines.Geo.genlinepts(start_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(ring.geoc - start_point)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)


                subpath = Routines.Geo.genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(ring.geoc - cross_point)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)
            elif i != 0 and i != len(fr)-2:
                prev_cross_point = cross_point
                bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: Routines.Geo.angle_btw(segend2-segend1, bc-segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
                subpath = Routines.Geo.genlinepts(prev_cross_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)

                subpath = Routines.Geo.genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)
            elif i == len(fr)-2:
                prev_cross_point = cross_point
                bond_centers = sorted([b.center for b in ring.bonds], key=lambda bc: Routines.Geo.angle_btw(segend2-segend1, bc-segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
                subpath = Routines.Geo.genlinepts(prev_cross_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)

                subpath = Routines.Geo.genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
                subpath = [p + normals[i]*height for p in subpath]
                pts += subpath
                ring_idx += [i]*len(subpath)

                subpath = Routines.Geo.genlinepts(cross_point, nb_ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - nb_ring.geoc)]
                subpath = [p + normals[i+1]*height for p in subpath]
                pts += subpath
                ring_idx += [i+1]*len(subpath)

                bond_centers = sorted([b.center for b in nb_ring.bonds], key=lambda bc: Routines.Geo.angle_btw(segend1-segend2, bc-segend2))
                cross_point, end_point = sorted([bond_centers[0], bond_centers[-1]], key=lambda c: np.linalg.norm(c-0.5*(segend2+segend1)))
                subpath = Routines.Geo.genlinepts(nb_ring.geoc, end_point, step_size)
                xnumbers += [np.linalg.norm(p-subpath[0])+xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(end_point - nb_ring.geoc)]
                subpath = [p + normals[i+1]*height for p in subpath]
                pts += subpath
                ring_idx += [i+1]*len(subpath)
        return pts, ring_idx, xnumbers, xticks

    def largest_fr(self):
        """
        reture the largest fused ring system
        :return: a list of ring objs
        """
        frlst = self.frlst()
        if len(frlst) < 0:
            sys.exit('no ring identified while you want to grep the largest fused ring system!')
        frlst.sort(key=lambda x: len(x), reverse=True)
        return frlst[0]

    @property
    def backbone(self):
        sites = []
        for r in self.largest_fr():
            for s in r.sites:
                if s not in sites:
                    sites.append(s)
        return SiteList(sites)


    def frlst(self):
        # fused rings, [[r1, r2, r3], [r4], [r5, r6], ...]
        indices = range(len(self.rings))
        block_list = []
        visited = []
        while len(visited) != len(self.rings):
            unvisited = [idx for idx in indices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in indices if idx not in block and idx not in visited]
                for i in outside:
                    if self.rings[block[pointer]].isfused_with(self.rings[i]):
                        block.append(i)
                visited.append(block[pointer])
                pointer += 1
            block_list.append([self.rings[j] for j in block])
        return block_list

    @classmethod
    def from_xyz(cls, xyz, name=None, charge=0, multiplicity=1):
        with open(xyz, 'r') as f:
            lines = f.readlines()[2:]
        sites = []
        for line in lines:
            items = line.strip().split()
            if len(items) == 4:
                s = Site(items[0], items[1:])
                sites.append(s)
        return cls(sites, name, charge, multiplicity)

    @classmethod
    def from_gjf(cls, gjf, name=None):
        with open(gjf, 'r') as f:
            lines = f.readlines()
        idx = 0
        hashline = ''
        for i in range(len(lines)):
            if lines[i].strip().startswith('#'):
                hashline = lines[i].strip()
                idx = i
                break
        if hashline == '':
            sys.exit('no hashline found in this gjf file!')

        charge, multiplicity = [int(s) for s in lines[idx+4].strip().split()]
        sites = []
        for i in range(idx + 5, len(lines)):
            items = lines[i].strip().split()
            if len(items) == 4:
                sites.append(Site(items[0], items[1:]))
        return cls(sites, name, charge, multiplicity, comment=hashline)



    #     self.nsites = len(self.sites)
    #
    #     for i in range(self.nsites):
    #         self.sites[i].mid = i
    #         self.sites[i].properties['num_neighbors'] = len(self.nrnbrmap[i])
    #         self.sites[i].properties['hybrid'] = OMol.check_hybid(self.sites[i].symbol, len(self.nrnbrmap[i]))
    #
    #     self.side_idx = []
    #     self.bone_idx = []
    #     # partition mol into backbone and side chains based on rings
    #     # TODO the group ring method cannot handle rubrene
    #     rings = mathop.loop_finding(self.nrnbrmap, 6) + mathop.loop_finding(self.nrnbrmap, 5) + \
    #             mathop.loop_finding(self.nrnbrmap, 4)
    #     ring_site_idx = [i for ring in rings for i in ring]
    #     largest_group = sorted(self.group_bonded_sites(self.bondmat, ring_site_idx), key=lambda x: len(x))[-1]
    #
    #     for i in range(len(self.sites)):
    #         # TODO not sound here, need a better way to identify backbone
    #         if len([idx for idx in self.nrnbrmap[i] if
    #                 self.sites[idx].properties['hybrid'] in ['sp2'] or self.sites[idx].symbol not in ['C', 'H']]) > 1 \
    #                 and i in largest_group:
    #             self.sites[i].properties['relative_position'] = 'bone'
    #             self.sites[i].properties['sidechain_label'] = [None, None]
    #             self.bone_idx.append(i)
    #         else:
    #             self.sites[i].properties['relative_position'] = 'side'
    #             self.sites[i].properties['sidechain_label'] = [-1, -1]
    #             self.side_idx.append(i)
    #
    #     # here we assume there's only one atom extended out from each backbone joint
    #     # TODO 68.cif  this is not the case
    #     # the property sidechain_label denotes [the label of the side chain, the level of that site on side chain]
    #     sidechain_counter = 0
    #     for i in range(len(self.sites)):
    #         if not set(self.nrnbrmap[i]).issubset(set(self.bone_idx)) and i in self.bone_idx:
    #             self.sites[i].properties['relative_position'] = 'bone-joint'
    #             for j in self.nrnbrmap[i]:
    #                 if j not in self.bone_idx:
    #                     self.sites[j].properties['relative_position'] = 'side-joint'
    #                     self.sites[i].properties['sidechain_label'] = [sidechain_counter, 0]
    #                     self.sites[j].properties['sidechain_label'] = [sidechain_counter, 1]
    #                     # break
    #             sidechain_counter += 1
    #     self.num_sidechains = sidechain_counter
    #
    #     # sitesop.sites_toxyz([self.sites[i] for i in self.bone_idx], 'bone.xyz')
    #
    #     # now assign sidechain_label to non-joint sites on side chains
    #     bonded_side_sites = OMol.group_bonded_sites(self.bondmat, self.side_idx)
    #     idxtrees = [None] * self.num_sidechains
    #     for grp in bonded_side_sites:
    #         grp_sidechain_counter = -1
    #
    #         # print(grp)
    #         # for i in grp:
    #         #     print(self.sites[i])
    #
    #         for idx in grp:
    #             if self.sites[idx].properties['sidechain_label'][0] != -1 and \
    #                     self.sites[idx].properties['sidechain_label'][1] == 1:
    #                 grp_sidechain_counter = self.sites[idx].properties['sidechain_label'][0]
    #                 break
    #
    #         bone_joint = [i for i in self.bone_idx if
    #                       self.sites[i].properties['relative_position'] == 'bone-joint'
    #                       and self.sites[i].properties['sidechain_label'][0] == grp_sidechain_counter][0]
    #
    #         for idx in grp:
    #             self.sites[idx].properties['sidechain_label'][0] = grp_sidechain_counter
    #
    #         # now assign the 2nd inx of sidechain label
    #         bone_joint_node = OMol.chain2sidetree(bone_joint, self.bondmat, grp)
    #         # this is redundant as levels have been calculated in chain2sidetree
    #         levels = [[node.name for node in children] for children in anytree.LevelOrderGroupIter(bone_joint_node)]
    #         for ilevel in range(len(levels)):
    #             for idx in levels[ilevel]:
    #                 self.sites[idx].properties['sidechain_label'][1] = ilevel
    #         idxtrees[grp_sidechain_counter] = bone_joint_node
    #
    #
    #     self.sitetrees = [OMol.idxtree2sitetree(tree, self) for tree in idxtrees]
    #
    #     self.backbone = Backbone(self, [self.sites[i] for i in self.bone_idx])
    #
    #     self.sidechains = [SideChain(self, tree) for tree in self.sitetrees]  # this muse be after backbone!
    #
    #     self.data = {
    #         'backbone' : self.backbone.data,
    #         'sidechains' : [sc.data for sc in self.sidechains],
    #     }
    #
    #
    #
    # @classmethod
    # def from_xyz(cls, xyz):
    #     with open(xyz, 'r') as f:
    #         lines = f.readlines()[2:]
    #     sites = []
    #     for line in lines:
    #         items = line.strip().split()
    #         if len(items) == 4:
    #             s = MSite(items[0], items[1:])
    #             sites.append(s)
    #     return cls(sites)
    #
    # @staticmethod
    # def idxtree2sitetree(idxtree, mol):
    #     symboltree = copy.deepcopy(idxtree)
    #     for node in anytree.LevelOrderIter(symboltree):
    #         node.site = mol.sites[node.name]
    #     return symboltree
    #
    # @staticmethod
    # def chain2sidetree(root_idx, bondmat, idxrange):
    #     """
    #     root_idx should only have one bonded idx in idxrange, this is the grow direction
    #     notice nodes are mutable
    #     # TODO add compatibility of closed loop
    #     :param root_idx: bone-joint
    #     :param bondmat:
    #     :param idxrange: the idx for one group, this does not include root_idx which is on the bone
    #     :return:
    #     """
    #     idxrange = [root_idx] + idxrange  # root_idx must come first
    #     # initialization
    #     classified = [[root_idx], ]
    #     pointer = 0
    #     while pointer != len(classified):
    #         outside = [idx for idx in idxrange if idx not in list(itertools.chain.from_iterable(classified))]
    #         tobeincluded = []
    #         for j in classified[pointer]:
    #             for i in outside:
    #                 if bondmat[j][i]:
    #                     tobeincluded.append(i)
    #         if len(tobeincluded) != 0:
    #             classified += [tobeincluded]
    #         pointer += 1
    #
    #     classified_nodes = []
    #     for i in classified:
    #         level = []
    #         for j in i:
    #             level.append(anytree.Node(j))
    #         classified_nodes.append(level)
    #
    #     # now convert to tree
    #     for i in reversed(range(1, len(classified_nodes))):
    #         for j in range(len(classified_nodes[i])):
    #             upper = i - 1
    #             for k in range(len(classified_nodes[upper])):
    #                 if bondmat[classified_nodes[i][j].name][classified_nodes[upper][k].name]:
    #                     classified_nodes[i][j].parent = classified_nodes[upper][k]
    #     return classified_nodes[0][0]
    #
    # @staticmethod
    # def group_bonded_sites(bondmat, idxrange):
    #     """
    #     :param bondmat: bondmat[i][j] returns whether ij are bonded
    #     :param idxrange:
    #     :return:
    #     """
    #     visited = []
    #     block_list = []
    #     while len(visited) != len(idxrange):
    #         # initialization
    #         unvisited = [idx for idx in idxrange if idx not in visited]
    #         ini_idx = unvisited[0]
    #         block = [ini_idx]
    #         pointer = 0
    #         while pointer != len(block):
    #             outside = [idx for idx in idxrange if idx not in block and idx not in visited]
    #             for i in outside:
    #                 if bondmat[block[pointer]][i]:
    #                     block.append(i)
    #             visited.append(block[pointer])
    #             pointer += 1
    #         block_list.append(block)
    #     return block_list
    #
    # @staticmethod
    # def isomol():
    #     # TODO identify whether a molecule can be considered as an OMol
    #     pass
    #
    # @staticmethod
    # def check_hybid(symbol, nnbs):
    #     if symbol == 'H':
    #         return 'hydrogen'
    #     elif symbol in ['O', 'S', 'Se']:
    #         if nnbs == 2:
    #             return 'sp3'
    #         elif nnbs == 1:
    #             return 'sp2'
    #         elif nnbs == 0:
    #             return 'isolated'
    #     elif symbol in ['N', 'P', 'As']:
    #         if nnbs == 3:
    #             return 'sp3'
    #         elif nnbs == 2:
    #             return 'sp2'
    #         elif nnbs == 1:
    #             return 'sp'
    #         elif nnbs == 0:
    #             return 'isolated'
    #     elif symbol in ['C', 'Si', 'Ge']:
    #         if nnbs == 4:
    #             return 'sp3'
    #         elif nnbs == 3:
    #             return 'sp2'
    #         elif nnbs == 2:
    #             return 'sp'
    #         elif nnbs == 1:
    #             return 'ion'
    #         elif nnbs == 0:
    #             return 'isolated'
    #
