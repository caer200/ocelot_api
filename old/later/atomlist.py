import itertools
import warnings
from copy import deepcopy

import numpy as np
import rdkit.Chem as Chem
from pymatgen.core.structure import Molecule
from scipy.spatial.distance import pdist, squareform, cdist

from ocelot.routines.conformerparser import ACParser
from ocelot.routines.geometry import rotate_along_axis, rotation_matrix, coord_transform, norm
from ocelot.schema.conformation.atomsite import AtomSite
from ocelot.schema.graph.elementnode import ElementNode
from ocelot.schema.graph.molgraph import MolGraph
from ocelot.schema.graph.nodecollection import NodeCollection


class AtomList:
    def __init__(self, atomsites):
        """
        this is the parent obj for backbone, sidegroups, mol, ring, bond.

        duplicate AtomSite will be discarded

        :param list atomsites: a *list* of AtomSite
        """
        unique_sites = []
        for ms in atomsites:
            if ms not in unique_sites:
                unique_sites.append(ms)
        self.atomsites = unique_sites

    @property
    def atomic_numbers(self):
        return [s.Z for s in self]

    @property
    def idxlist(self):
        return [s.siteid for s in self]

    @property
    def i2siteid(self):
        d = {}
        for i in range(len(self)):
            d[i] = self[i].siteid
        return d

    @property
    def siteid2i(self):
        self.checkstatus('all assigned', 'unique ids')
        d = {}
        for i in range(len(self)):
            d[self.idxlist[i]] = i
        return d

    def __len__(self):
        return len(self.atomsites)

    def __contains__(self, atomsite):  # necessary?
        for s in self.atomsites:
            if s == atomsite:
                return True
        return False

    def __iter__(self):
        return self.atomsites.__iter__()

    def __getitem__(self, ind):
        return self.atomsites[ind]

    def __repr__(self):
        outs = [self.__class__.__name__ + ': ']
        for s in self.atomsites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @property
    def coordmat(self):
        """
        coordinates matrix

        :return: nx3 np array
        """
        coordmat = np.empty((len(self), 3))
        for i in range(len(self)):
            for j in range(3):
                coordmat[i][j] = self[i].coords[j]
        return coordmat

    def assign_siteid(self):
        for i in range(len(self)):
            self[i].siteid = i

    @staticmethod
    def get_nbrmap(bmat):
        """
        :param np.ndarray bmat: bool bond matrix, i is not bonded to itself
        :return: nbmap[i] is the index list of i's neighbors
        """
        ma = {}
        numsites = len(bmat)
        for i in range(numsites):
            nbs = []
            for j in range(numsites):
                if bmat[i][j] and j != i:
                    nbs.append(j)
            ma[i] = nbs
        return ma

    @staticmethod
    def get_distmat(coordmat):
        """
        distanct matrix

        :param coordmat: coordinates matrix nx3
        :return: distmat[i][j] is the euclid distance between coordmat[i] and coordmat[j]
        """
        return squareform(pdist(coordmat))

    def get_closest_sites(self, other):
        """
        other_sites -- self_border_site --- distmin --- other_border_site -- other_sites

        this can be used to see if a dimer is interesting or not

        :param AtomList other:
        :return: i, j, distmin
        self[i] is self_border_site,
        self[j] is other_border_site,
        """
        distmat = cdist(self.coordmat, other.coordmat)
        minid = np.unravel_index(np.argmin(distmat, axis=None), distmat.shape)
        i_self_border_site = distmat[minid[0]]
        j_other_border_site = distmat[minid[1]]
        distmin = distmat[minid]
        return i_self_border_site, j_other_border_site, distmin

    def get_site_by_coords(self, coords, tol=1e-5):
        """
        will only return one site

        :param coords:
        :param tol:
        :return:
        """
        for s in self.atomsites:
            if np.allclose(s.coords, coords, atol=tol):
                return s
        return None

    @staticmethod
    def get_bondmat(atomsites, distmat, co=1.3):
        """
        Bij = whether there is a bond between si and sj
        i is NOT bonded with itself

        :param co: coefficient for cutoff, default 1.3, based on covalent rad
        :return: bool matrix
        """
        numsites = len(atomsites)
        mat = np.ones((numsites, numsites), dtype=bool)
        for i in range(numsites):
            mat[i][i] = False
            for j in range(i + 1, numsites):
                if distmat[i][j] > (atomsites[i].element.atomic_radius + atomsites[j].element.atomic_radius) * co:
                    mat[i][j] = False
                    mat[j][i] = False
        return mat

    @property
    def elements(self):
        return [s.element.symbol for s in self]

    @property
    def geoc(self):
        """
        geometric center

        :return: 3x1 np array
        """
        v = np.zeros(3)
        for s in self.atomsites:
            v += s.coords
        return v / len(self.atomsites)

    @property
    def volume_slow(self, boxdensity=0.2):
        """
        http://wiki.bkslab.org/index.php/Calculate_volume_of_the_binding_site_and_molecules

        First, Lay a grid over the spheres.

        Count the number or points contained in the spheres (Ns).

        Count the number of points in the grid box (Ng).

        Calculate the volume of the grid box (Vb).

        don't use this as it's slow...

        :return: volume in \AA^3
        """
        mat = np.empty((len(self.atomsites), 4))
        for i in range(len(self.atomsites)):
            for j in range(3):
                mat[i][j] = self.atomsites[i].coords[j]
            mat[i][3] = self.atomsites[i].element.covrad

        box_min = np.floor(min(itertools.chain.from_iterable(mat))) - 2
        box_max = np.ceil(max(itertools.chain.from_iterable(mat))) + 2
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

    def intersection(self, other, copy=False):
        """
        :param AtomList other:
        :param bool copy: whether generate a list of deepcopied sites or not
        :return: a list of atomsites in self that belong to both AtomList
        """
        r = []
        if copy:
            for s in self.atomsites:
                if s in other:
                    r.append(deepcopy(s))
        else:
            for s in self.atomsites:
                if s in other:
                    r.append(s)
        return r

    def get_site_byid(self, siteid, multimatch=False):
        """
        :param int siteid:
        :param bool multimatch: if Ture, return a list, otherwise a msite
        :return: (a list of) msite obj
        """
        if multimatch:
            sites_matches = []
            for s in self.atomsites:
                if s.siteid == siteid:
                    sites_matches.append(s)
            if len(sites_matches) == 0:
                warnings.warn('W: cannot get site by id: {}'.format(siteid))
                return None
            return sites_matches
        else:
            for s in self.atomsites:
                if s.siteid == siteid:
                    return s

    def issubset(self, other):
        """
        based on atomsite.coords and atomsite.element, see __eq__ in AtomSite

        :param AtomList other:
        :rtype: bool
        """
        return len(self.atomsites) == len(self.intersection(other))

    def as_dict(self):
        """
        keys are

        atomsites, can, volume

        :return: a dict representation
        """
        d = {
            "atomsites": [s.as_dict() for s in self.atomsites],
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__
        }
        return d

    @classmethod
    def from_coordmat(cls, mat, names, ids=None):
        """
        :param mat: nx3 array
        :param names: a list of string as element_name
        :param ids:
        :return:
        """
        if ids is None:
            ids = [None] * len(mat)
        ss = []
        for i in range(len(names)):
            ss.append(AtomSite(names[i], mat[i], ids[i]))
        return cls(ss)

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        atomsites

        :param dict d: d['atomsites'] is a list of atomsites' dicts
        """
        return cls([AtomSite.from_dict(entry) for entry in d["atomsites"]])

    @classmethod
    def from_coordmat_and_elementlst_and_siteidlst(cls, mat, elementlst, ids=None):
        """
        :param np.ndarray mat: coords matrix nx3
        :param list elementlst: a list of string as element_name
        :param list ids: siteid list
        """
        if ids is None:
            ids = [None] * len(mat)
        if not len(mat) == len(elementlst) == len(ids):
            warnings.warn('W: length of the args are inconsistent')
            return None
        ss = []
        for i in range(len(elementlst)):
            ss.append(AtomSite(elementlst[i], mat[i], ids[i]))
        return cls(ss)

    @classmethod
    def from_file(cls, fname):
        """
        using pmg molecule

        :param str fname:
        """
        m = Molecule.from_file(fname)
        return cls([AtomSite.from_pymatgen_site(s) for s in m.sites])

    def to(self, ftype, fname):
        """
        using pmg molecule

        :param str ftype: 'xyz', 'gjf', etc.
        :param str fname: with extension
        """
        pymatgen_sites = [s.to_pymatgen_site() for s in self.atomsites]
        m = Molecule.from_sites(pymatgen_sites)
        m.to(ftype, fname)

    @property
    def pmgmol(self):
        """
        notice siteid info lost here

        :return: pymatgen molecule obj
        """
        pymatgen_sites = [s.to_pymatgen_site() for s in self.atomsites]
        return Molecule.from_sites(pymatgen_sites)

    def rotate_along(self, theta, end1, end2, unit='degree'):
        """
        for each AtomSite in the AtomList,
        rotate the vector defined by (site.coords - end1) along (end2 - end1) and add this vector to site.coords

        end2 ------ site

        \|

        \|

        \|

        end1

        notice the coords are changed in-place

        :param theta: angle of rotation
        :param end1: 3x1 float list/array
        :param end2: 3x1 float list/array
        :param str unit: degree/radian
        """
        end1 = np.array(end1)
        end2 = np.array(end2)
        for s in self.atomsites:
            v = s.coords - end1
            s.coords = end1 + rotate_along_axis(v, end2 - end1, theta, thetaunit=unit)

    def rotate_along_de_matrix(self, theta, end1, end2, unit='degree'):
        """
        rotation matrix for :func:`~AtomList.rotate_along`
        """
        end1 = np.array(end1)
        end2 = np.array(end2)
        axis = end2 - end1
        return rotation_matrix(axis, theta, thetaunit=unit)

    def rotate_along_with_matrix(self, matrix, end1):
        """
        a quicker version rotate_along if we konw the rotation matrix, end2 is included in the matrix
        """
        end1 = np.array(end1)
        for s in self.atomsites:
            v = s.coords - end1
            s.coords = end1 + np.dot(matrix, v)

    @classmethod
    def orient(cls, atomlist, pqo):
        """
        p, q, o are 3 orthonormal vectors in x, y, z basis

        basis transformation from xyz to pqo

        :param pqo: 3x3 array
        """
        pqo = np.array(pqo)
        p, q, o = pqo
        newosites = []
        for s in atomlist.atomsites:
            ns = deepcopy(s)
            ns.coords = coord_transform(p, q, o, s.coords - atomlist.geoc)
            newosites.append(ns)
        return cls(newosites)

    @property
    def status(self):
        """
        check status based on self.idxlist

        :return:
        """
        s = []

        if self.idxlist[0] == 0:
            s.append('0start')
        else:
            s.append('not0strat')

        if set(self.idxlist) == {None}:
            s.append('all none')
        elif None in set(self.idxlist):
            s.append('partially none')
        else:
            s.append('all assigned')

        if all(a + 1 == b for a, b in zip(self.idxlist, self.idxlist[1:])):
            s.append('continuous')
        else:
            s.append('not continuous')

        if len(set(self.idxlist)) < len(self.idxlist):
            s.append('duplicate ids')
        else:
            s.append('unique ids')

        return s

    def assign_siteid_default(self):
        for i in range(len(self)):
            self[i].siteid = i

    def get_bmat_based_on_siteid(self, co=1.3):
        self.checkstatus('all assigned', 'unique ids')
        distmat = self.get_distmat(self.coordmat)

        disti_to_siteid = {}
        for i in range(len(distmat)):
            disti_to_siteid[i] = self.atomsites[i].siteid
        siteid_to_disti = {v: k for k, v in disti_to_siteid.items()}

        mat = np.zeros((max(self.idxlist) + 10, max(self.idxlist) + 10), dtype=bool)
        for sidi in self.idxlist:
            i = siteid_to_disti[sidi]
            for sidj in self.idxlist:
                j = siteid_to_disti[sidj]
                distance = distmat[i][j]
                if distance <= (self[i].element.atomic_radius + self[j].element.atomic_radius) * co:
                    mat[sidi][sidj] = 1
                    mat[sidj][sidi] = 1
        return mat

    def to_nodecollection(self):
        self.checkstatus('all assigned', 'unique ids')
        sites = sorted(self.atomsites, key=lambda x: x.siteid)
        nodes = [ElementNode(s.element.symbol, s.siteid) for s in sites]
        return NodeCollection(nodes)

    def to_molgraph(self):
        """
        get a MolGraph, siteid --> nodename

        otherwise node index is assigned based on the order in self.atomsites
        """
        bmat_based_on_siteid = self.get_bmat_based_on_siteid(co=1.3)
        nc = self.to_nodecollection()
        return MolGraph.from_nodecollection(nc, bmat_based_on_siteid)

    def to_rdmol(self, charge=0, sani=True, charged_fragments=True, force_single=False, expliciths=True):
        """
        generate a rdmol obj with current conformer

        siteid --dict--> atomidx

        :param charge:
        :param sani:
        :param charged_fragments:
        :param force_single:
        :param expliciths:
        :return:
        """
        siteid2atomidx = self.siteid2i
        atomidx2siteid = self.i2siteid
        conf = Chem.Conformer(len(self))
        coordmat = self.coordmat
        for i in range(len(self)):
            conf.SetAtomPosition(i, (coordmat[i][0], coordmat[i][1], coordmat[i][2]))
        ac = self.get_bondmat(self.atomsites, self.get_distmat(self.coordmat), co=1.3)
        ap = ACParser(ac, charge, self.atomic_numbers, sani=sani)
        rdmol, smiles = ap.parse(charged_fragments=charged_fragments, force_single=force_single, expliciths=expliciths)
        rdmol.AddConformer(conf)
        return rdmol, smiles, siteid2atomidx, atomidx2siteid

    def checkstatus(self, *args):
        for a in args:
            if a in self.status:
                pass
            else:
                raise SiteidError('no <{}> in status!'.format(a))


class SiteidError(Exception):
    pass
