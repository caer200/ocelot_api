import itertools
import warnings
from copy import deepcopy
from typing import List

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import rdkit.Chem as Chem
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import Site
from pymatgen.io.xyz import XYZ
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolAlign
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from shapely.geometry import Polygon

from ocelot.routines.conformerparser import ACParser
from ocelot.routines.geometry import Fitter
from ocelot.routines.geometry import alpha_shape
from ocelot.routines.geometry import angle_btw
from ocelot.routines.geometry import coord_transform
from ocelot.routines.geometry import get_proj_point2plane
from ocelot.routines.geometry import norm
from ocelot.routines.geometry import rotate_along_axis
from ocelot.routines.geometry import rotation_matrix
from ocelot.routines.geometry import unify
from ocelot.routines.pbc import AtomicRadius
from ocelot.schema.graph import BackboneGraph
from ocelot.schema.graph import BasicGraph
from ocelot.schema.graph import FragmentGraph
from ocelot.schema.graph import MolGraph
from ocelot.schema.graph import SidechainGraph
from ocelot.schema.rdfunc import RdFunc

_coordination_rule = {
    'H': 1,
    'Li': 1,
    'Na': 1,
    'K': 1,
    'Be': 2,
    'Mg': 2,
    'Ca': 2,
    'B': 3,
    'Al': 3,
    'Ga': 3,
    'C': 4,
    'Si': 4,
    'Ge': 4,
    'Sn': 4,
    'N': 3,
    'P': 3,
    'As': 3,
    'O': 2,
    'S': 2,
    'Se': 2,
    'F': 1,
    'Cl': 1,
    'Br': 1,
    'I': 1,
}


class ConformerError(Exception):
    pass

class SiteidError(Exception):
    pass

_rdmol_sani = True
_rdmol_force_single = False
_rdmol_charged_fragments = False
_rdmol_expliciths = True


class SiteidOperation:

    def __getitem__(self, i: int):
        return self.sites[i]

    def __len__(self):
        return len(self.sites)

    def __iter__(self):
        return self.sites.__iter__()

    def __init__(self, sites: [Site]):
        self.sites = sites

        # allkeys = ['siteid']
        # for s in self.sites:
        #     allkeys += list(s.properties.keys())
        # allkeys = list(set(allkeys))
        #
        # for i in range(len(self)):
        #     k = 'siteid'
        #     # for k in allkeys:
        #     try:
        #         self[i].properties[k]
        #     except KeyError:
        #         self[i].properties[k] = None

    @property
    def sdict(self):
        d = {}
        for s in self.sites:
            d[s.properties['siteid']] = s
        return d

    @property
    def siteids(self):
        try:
            return [s.properties['siteid'] for s in self]
        except KeyError:
            raise SiteidError('at least one site does not have siteid field!')

    def get_site_byid(self, siteid):
        """
        :param int siteid:
        :return: (a list of) msite obj
        """
        return self.sdict[siteid]
        # if multimatch:
        #     sites_matches = []
        #     for s in self:
        #         if s.properties['siteid'] == siteid:
        #             sites_matches.append(s)
        #     if len(sites_matches) == 0:
        #         raise SiteidError('cannot get site by id: {}'.format(siteid))
        #     return sites_matches
        # else:
        #     for s in self:
        #         if s.properties['siteid'] == siteid:
        #             return s

    def get_sites_byids(self, siteids, copy=False):
        self.checkstatus('all assigned', 'unique ids')
        if not set(siteids).issubset(set(self.siteids)):
            raise SiteidError('siteids is not a subset of sitelist when init conformer')
        rs = []
        for i in siteids:
            s = self.sdict[i]
            # for s in self:
            # if s.properties['siteid'] in siteids:
            if copy:
                rs.append(deepcopy(s))
            else:
                rs.append(s)
        return rs

    @property
    def status(self):
        """
        check status based on conformer.siteids

        :return:
        """
        s = []

        # if self.siteids[0] == 0:
        #     s.append('0start')
        # else:
        #     s.append('not0strat')

        if set(self.siteids) == {None}:
            s.append('all none')
        elif None in set(self.siteids):
            s.append('partially none')
        else:
            s.append('all assigned')

        # try:
        #     cri = all(a + 1 == b for a, b in zip(self.siteids, self.siteids[1:]))
        #     if cri:
        #         s.append('continuous')
        #     else:
        #         s.append('not continuous')
        # except TypeError:
        #     pass

        if len(set(self.siteids)) < len(self.siteids):
            s.append('duplicate ids')
        else:
            s.append('unique ids')

        return s

    def assign_siteid(self, siteids):
        """
        :param siteids: default None means siteid = index, otherwise siteid = siteids[i]
        :return:
        """
        if siteids is False:
            pass
        elif siteids is None:
            for i in range(len(self)):
                self[i].properties['siteid'] = i
        elif len(siteids) == len(self) and all(isinstance(i, int) for i in siteids):
            for i in range(len(self)):
                self[i].properties['siteid'] = siteids[i]
        else:
            raise SiteidError('siteids are not legit!')

    def checkstatus(self, *args):
        for a in args:
            if a in self.status:
                pass
            else:
                raise SiteidError('no <{}> in status!'.format(a))


class BasicConformer(SiteidOperation):

    def __init__(self, sites, siteids=False):
        """
        :param sites:
        :param siteids: a list of siteids, siteids[index]=siteid of self[index]
        """
        super().__init__(sites)
        self.assign_siteid(siteids)
        self.checkstatus('all assigned', 'unique ids')

    @property
    def composition(self):
        return self.pmgmol.composition

    def __eq__(self, other):
        # using center should be enough for sites in a mol
        if self.pmgmol == other.pmgmol:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def compare(self, other):
        """
        get similarity based on rmsd, return inf if two different molecules
        """
        if len(self) != len(other):
            return np.inf
        if set(self.atomic_numbers) != set(other.atomic_numbers):
            return np.inf
        rdmol1, smiles1, siteid2atomidx, atomidx2siteid = self.to_rdmol()
        rdmol2, smiles2, siteid2atomidx, atomidx2siteid = other.to_rdmol()
        if smiles1 != smiles2:
            return np.inf
        rdmol1 = Chem.RemoveHs(rdmol1)
        rdmol2 = Chem.RemoveHs(rdmol2)
        return rdMolAlign.GetBestRMS(rdmol1, rdmol2)

    @property
    def pmgmol(self):
        return Molecule.from_sites(self.sites)

    def __repr__(self):
        outs = ["{}:".format(self.__class__.__name__)]
        for s in self:
            outs.append(s.__repr__() + '\t' + s.properties.__repr__())
        return "\n".join(outs)

    @property
    def cart_coords(self):
        coords = np.zeros((len(self), 3))
        for i in range(len(self)):
            coords[i] = self[i].coords
        return coords

    @property
    def index2siteid(self):
        d = {}
        for i in range(len(self)):
            d[i] = self[i].siteid
        return d

    @property
    def siteid2index(self):
        self.checkstatus('all assigned', 'unique ids')
        d = {}
        for i in range(len(self)):
            d[self.siteids[i]] = i
        return d

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

    @property
    def nbrmap(self):
        return self.get_nbrmap(self.bondmat)

    @staticmethod
    def get_distmat(coordmat):
        """
        distanct matrix

        :param coordmat: coordinates matrix nx3
        :return: distmat[i][j] is the euclid distance between coordmat[i] and coordmat[j]
        """
        return squareform(pdist(coordmat))

    @property
    def distmat(self):
        return self.get_distmat(self.cart_coords)

    def get_closest_sites(self, other):
        """
        other_sites -- self_border_site --- distmin --- other_border_site -- other_sites

        this can be used to see if a dimer is interesting or not

        :param AtomList other:
        :return: i, j, distmin
        conformer[i] is self_border_site,
        conformer[j] is other_border_site,
        """
        distmat = cdist(self.cart_coords, other.cart_coords)
        minid = np.unravel_index(np.argmin(distmat, axis=None), distmat.shape)
        i_self_border_site = distmat[minid[0]]
        j_other_border_site = distmat[minid[1]]
        distmin = distmat[minid]
        return i_self_border_site, j_other_border_site, self[i_self_border_site], other[j_other_border_site], distmin

    def get_site_by_coords(self, coords, tol=1e-5):
        """
        will only return one site

        :param coords:
        :param tol:
        :return:
        """
        for s in self:
            if np.allclose(s.coords, coords, atol=tol):
                return s
        return None

    @staticmethod
    def get_bondmat(sites, distmat, co=1.3):
        """
        Bij = whether there is a bond between si and sj, i is NOT bonded with itself

        if site is not a legit atom (e.g. bq in nics), it cannot have any bond

        :param sites:
        :param distmat:
        :param co: coefficient for cutoff, default 1.3, based on covalent rad
        :return: bool matrix
        """
        numsites = len(sites)
        mat = np.zeros((numsites, numsites), dtype=bool)
        for i in range(numsites):
            irad = AtomicRadius(sites[i])
            for j in range(i + 1, numsites):
                jrad = AtomicRadius(sites[j])
                cutoff = (irad + jrad) * co
                if 1e-5 < distmat[i][j] < cutoff:
                    mat[i][j] = True
                    mat[j][i] = True
        return mat

    @property
    def bondmat(self, co=1.3):
        return self.get_bondmat(self.sites, self.distmat, co)

    @property
    def symbols(self):
        return [s.species_string for s in self]

    @property
    def geoc(self):
        """
        geometric center

        :return: 3x1 np array
        """
        v = np.zeros(3)
        for s in self:
            v += s.coords
        return v / len(self)

    @property
    def volume_slow(self, boxdensity=0.2):
        """
        http://wiki.bkslab.org/index.php/Calculate_volume_of_the_binding_site_and_molecules

        First, Lay a grid over the spheres.

        Count the number or points contained in the spheres (Ns).

        Count the number of points in the grid box (Ng).

        Calculate the volume of the grid box (Vb).

        don't use this as it's slow...

        :return: volume in A^3
        """
        mat = np.empty((len(self), 4))
        for i in range(len(self)):
            for j in range(3):
                mat[i][j] = self[i].coords[j]
            mat[i][3] = self[i].element.covrad

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

    @property  # type: ignore
    def atomic_numbers(self):
        """List of atomic numbers."""
        return tuple(site.specie.Z for site in self)  # type: ignore

    # @property
    # def atomic_radii(self):
    #     return (Element(site.species_string).atomic_radius for site in self)
    # return tuple(atomic_radius Elesite.species_string for site in self)

    def intersection(self, other, copy=False):
        """
        :param other:
        :param bool copy: whether generate a list of deepcopied sites or not
        :return: a list of sites in conformer that belong to both
        """
        r = []
        if copy:
            for s in self:
                if s in other:
                    r.append(deepcopy(s))
        else:
            for s in self:
                if s in other:
                    r.append(s)
        return r

    def issubset(self, other):
        """
        :param other:
        :rtype: bool
        """
        return len(self) == len(self.intersection(other))

    def sort_by_siteid(self):
        self.checkstatus('all assigned')
        self.sites = sorted(self.sites, key=lambda x: x.properties['siteid'])

    def rotate_along(self, theta, end1, end2, unit='degree'):
        """
        for each site, rotate the vector defined by (site.coords - end1) along (end2 - end1) and add this vector to site.coords

        end1 - end2
        |
        site

        notice the coords are changed in-place

        :param theta: angle of rotation
        :param end1: 3x1 float list/array
        :param end2: 3x1 float list/array
        :param str unit: degree/radian
        """
        end1 = np.array(end1)
        end2 = np.array(end2)
        for s in self:
            v = s.coords - end1
            s.coords = end1 + rotate_along_axis(v, end2 - end1, theta, thetaunit=unit)

    @staticmethod
    def rotate_along_de_matrix(theta, end1, end2, unit='degree'):
        """
        rotation matrix :func:`~BasicConformer.rotate_along`
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
        for s in self:
            v = s.coords - end1
            s.coords = end1 + np.dot(matrix, v)

    def orient(self, pqo, origin='geoc'):
        """
        p, q, o are 3 orthonormal vectors in x, y, z basis

        basis transformation from xyz to pqo, notice sites are deepcopied

        :param pqo: 3x3 array
        :param origin: default geoc, otherwise a 3d coord
        """
        pqo = np.array(pqo)
        p, q, o = pqo
        newsites = []
        if origin == 'geoc':
            origin = self.geoc
        origin = np.array(origin)
        for s in self:
            ns = deepcopy(s)
            ns.coords = coord_transform(p, q, o, s.coords - origin)
            newsites.append(ns)
        return self.from_sites(newsites)

    def get_nbrmap_based_on_siteid(self, co=1.3):
        bmat = self.get_bondmat_based_on_siteid(co)
        ma = {}
        for i in self.siteids:
            nbs = []
            for j in self.siteids:
                if bmat[i][j] and j != i:
                    nbs.append(j)
            ma[i] = nbs
        return ma

    @property
    def can_rdmol(self):
        try:
            self.to_rdmol()
            return True
        except:
            return False

    def get_bondmat_based_on_siteid(self, co=1.3):
        # self.checkstatus('all assigned', 'unique ids')
        distmat = self.distmat
        siteid_to_disti = self.siteid2index
        mat = {}
        for ii in self.siteids:
            mat[ii] = {}
            for jj in self.siteids:
                mat[ii][jj] = False
            # print(self.siteids)
        # mat = np.zeros((max(self.siteids) + 10, max(self.siteids) + 10), dtype=bool)
        for isidi in range(len(self)):
            sidi = self.siteids[isidi]
            i = siteid_to_disti[sidi]
            irad = AtomicRadius(self[i])
            for isidj in range(isidi + 1, len(self)):
                sidj = self.siteids[isidj]
                j = siteid_to_disti[sidj]
                jrad = AtomicRadius(self[j])
                distance = distmat[i][j]
                cutoff = (irad + jrad) * co
                if 1e-5 < distance < cutoff:
                    mat[sidi][sidj] = True
                    mat[sidj][sidi] = True
        return mat

    def to_graph(self, nodename='siteid', graphtype='MolGraph', partition_scheme=None, joints: dict = None):
        """
        get a Graph, default siteid --> nodename

        otherwise nodename is assigned based on the order in self.sites
        """
        bondmat_by_siteid = self.get_bondmat_based_on_siteid(co=1.3)
        g = nx.Graph()
        if nodename == 'siteid':
            for i in range(len(self)):
                s = self[i]
                g.add_node(s.properties['siteid'], symbol=s.species_string)
                for j in range(i + 1, len(self)):
                    ss = self[j]
                    g.add_node(ss.properties['siteid'], symbol=ss.species_string)
                    if bondmat_by_siteid[s.properties['siteid']][ss.properties['siteid']]:
                        g.add_edge(s.properties['siteid'], ss.properties['siteid'])
        if graphtype == 'MolGraph':
            return MolGraph(g)
        elif graphtype == 'FragmentGraph':
            return FragmentGraph(g, joints=joints, partition_scheme=partition_scheme)
        elif graphtype == 'BackboneGraph':
            return BackboneGraph(g, joints=joints, partition_scheme=partition_scheme)
        elif graphtype == 'SidechainGraph':
            return SidechainGraph(g, joints=joints, partition_scheme=partition_scheme)
        return BasicGraph(g)

    def is_missing_hydrogen(self):
        try:
            rdmol, smiles, siteid2atomidx, atomidx2siteid = self.to_rdmol()  # hydrogens are explicitly added
        except:
            warnings.warn('to_rdmol failed for this conformer! we believe it misses hydrogen')
            return True
        # to_rdmol is quite sensitive to AC, if AC is wrong and charged_frag=False then nradical will be quite large
        nradical = Descriptors.NumRadicalElectrons(rdmol)
        if nradical != 0:
            return True
        rdmolh = Chem.AddHs(rdmol)
        if len(rdmol.GetAtoms()) != len(rdmolh.GetAtoms()):
            return True
        return False

    def to_rdmol(self, charge=0, sani=True, charged_fragments=None, force_single=False, expliciths=True):
        """
        generate a rdmol obj with current conformer

        siteid --dict--> atomidx in rdmol == index in conformer

        :param charge:
        :param sani:
        :param charged_fragments:
        :param force_single:
        :param expliciths:
        :return:
        """
        siteid2atomidx = self.siteid2index
        atomidx2siteid = self.index2siteid
        conf = Chem.Conformer(len(self))
        coordmat = self.cart_coords
        for i in range(len(self)):
            conf.SetAtomPosition(i, (coordmat[i][0], coordmat[i][1], coordmat[i][2]))

        if hasattr(self, 'joints'):  # LBYL
            apriori_radicals = {}
            for k in self.joints.keys():
                apriori_radicals[siteid2atomidx[k]] = len(self.joints[k])
        else:
            apriori_radicals = None
        ac = self.bondmat
        ap = ACParser(ac, charge, self.atomic_numbers, sani=sani, apriori_radicals=apriori_radicals)
        if charged_fragments is None:
            try:
                rdmol, smiles = ap.parse(charged_fragments=False, force_single=force_single, expliciths=expliciths)
            except Chem.rdchem.AtomValenceException:
                warnings.warn('AP parser cannot use radical scheme, trying to use charged frag')
                rdmol, smiles = ap.parse(charged_fragments=True, force_single=force_single, expliciths=expliciths)
        else:
            rdmol, smiles = ap.parse(charged_fragments=charged_fragments, force_single=force_single, expliciths=expliciths)
        rdmol.AddConformer(conf)
        return rdmol, smiles, siteid2atomidx, atomidx2siteid

    def to(self, fmt, filename):
        self.pmgmol.to(fmt, filename)

    @classmethod
    def from_sites(cls, sites, siteids=False):
        """
        we use this as the parent constructor
        """
        return cls(sites, siteids)

    @classmethod
    def from_pmgmol(cls, m: Molecule, siteids=False):
        return cls.from_sites(m.sites, siteids)

    @classmethod
    def from_siteids(cls, siteids, sites, copy=False):
        """
        build up a ring based on a list of siteid and a list of sites

        the sites should have been assigned siteids
        """
        rs = SiteidOperation(sites).get_sites_byids(siteids, copy)
        return cls.from_sites(rs, siteids=False)

    @classmethod
    def from_file(cls, filename):
        m = Molecule.from_file(filename)
        c = cls.from_pmgmol(m, siteids=None)
        warnings.warn(
            'conformer built from file, siteids are assigned by the order of entries in file: {}'.format(filename))
        return c

    @classmethod
    def from_str(cls, input_string, fmt):
        m = Molecule.from_str(input_string, fmt)
        conformer = cls.from_pmgmol(m, siteids=None)
        warnings.warn('conformer built from string, siteids are assigned by the order of entries in the string')
        return conformer

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "sites": [s.as_dict() for s in self],
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        sites = [Site.from_dict(sd) for sd in d['sites']]
        return cls.from_sites(sites, siteids=False)


class BondConformer(BasicConformer):

    def __init__(self, sites, siteids=False):
        super().__init__(sites, siteids)
        # self.checkstatus('all assigned', 'unique ids')
        if len(self) != 2:
            raise ConformerError('only use two sites to create a bond! you are using {}'.format(len(self)))
        self.a = self[0]
        self.b = self[1]

    @property
    def length(self):
        return norm(self.a.coords - self.b.coords)


class RingConformer(BasicConformer):
    def __init__(self, sites, siteids=False):
        super().__init__(sites, siteids)
        # self.checkstatus('all assigned', 'unique ids')
        self.n_member = len(self)
        if self.n_member < 3:
            raise ConformerError('you are initializing a ring with less than 3 sites!')
        normal, self.ptsmean, self.pl_error = Fitter.plane_fit(self.cart_coords)
        self.n1 = np.array(normal)
        self.n2 = np.array(-normal)

    @property
    def bonds_in_ring(self):
        """
        bond objects can be extracted from this ring
        e.g. for benzene the C-H bonds are NOT here

        :return: a list of bonds
        """
        bmat = self.bondmat
        bonds = []
        for ij in itertools.combinations(range(len(self)), 2):
            i, j = ij
            if bmat[i][j]:
                b = BondConformer.from_sites([self[i], self[j]])
                if b not in bonds:
                    bonds.append(b)
        return bonds

    def normal_along(self, refnormal, tol=45.0):
        """
        get the n1 or n2 that is along the direction of refnormal within a certain tol
        this is useful to identify 2 plane normals for a partially-bent structure

        :param refnormal: 3x1 array
        :param float tol: default 45 in degree
        :return: None if no normal along refnormal found
        """
        for n in [self.n1, self.n2]:
            if abs(angle_btw(n, refnormal)) < np.radians(tol):
                return n
        warnings.warn('W: no normal along this direction!')
        return None

    @staticmethod
    def get_avg_norms(rings):
        ref_ring = rings[0]
        avg_norm = np.zeros(3)
        for ring in rings:
            normal_along_avenorm = ring.normal_along(ref_ring.n1)
            if normal_along_avenorm is not None:
                avg_norm += ring.normal_along(ref_ring.n1)
        return unify(avg_norm), -unify(avg_norm)

    def iscoplane_with_norm(self, v2, tol=20.0, tolunit='degree'):
        """
        whether two rings are on the same plane with tol

        :param v2:
        :param tol: degree default 20
        :param tolunit: degree/radian
        :return: bool
        """
        angles = []
        if tolunit == 'degree':
            tol = np.radians(tol)
        for v1 in [self.n1, self.n2]:
            angles.append(abs(angle_btw(v1, v2)))
        if min(angles) < tol:
            return True
        return False

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

    def as_dict(self):
        d = super().as_dict()
        d['n_member'] = self.n_member
        # d['bonds'] = [b.as_dict() for b in self.bonds_in_ring]
        return d


def conformer_addh(c: BasicConformer, joints=None, original: BasicConformer = None):
    """
    all sites will be deep copied

    if we know nothing about the joints, we have to find the under-coordinated sites, then terminate them based on # of
    valence electrons

    if we know the joints as a list or set of siteids but not the original conformer, we just terminate them based on # of valence electrons

    if we know the joints as a dict but not the original conformer, we terminate them based on len(joints[siteid])

    if we know the joints and the orginal conformer, we can terminate them based on # of bonds broken during fragmenting, in this case joints is a dict

    :param c:
    :param joints:
    :param original:
    :return: d[siteid_of_the_joint] = a list of hydrogen sites
    """
    if joints is None:
        warnings.warn('adding h based on undercoordination')
        nbrmap = c.nbrmap
        hsites_by_joint_siteid = {}
        for i in range(len(c)):
            s = c[i]
            hsites = []
            try:
                nvalence = _coordination_rule[s.species_string]
            except KeyError:
                nvalence = 0
            if nvalence > len(nbrmap[i]):  # we got a joint
                nh_to_be_added = nvalence - len(nbrmap[i])
                if nh_to_be_added == 1 or nh_to_be_added == 3:
                    v_s_to_h = np.zeros(3)
                    for nb_idx in nbrmap[i]:
                        v_nb_to_s = s.coords - c[nb_idx].coords
                        v_s_to_h += v_nb_to_s
                    v_s_to_h = unify(v_s_to_h)
                    hcoords = s.coords + v_s_to_h * 1.1
                    hsite = Site('H', hcoords)
                    hsites.append(hsite)
                elif nh_to_be_added == 2:
                    if len(nbrmap[i]) == 2:
                        nb1 = c[nbrmap[i][0]]
                        nb2 = c[nbrmap[i][1]]
                    elif len(nbrmap[i]) == 1:
                        nb1 = c[nbrmap[i][0]]
                        nb2 = s.coords + np.array([1, 0, 0])
                        if unify(np.cross(nb1 - s.coords, nb2 - s.coords)) < 1e-5:
                            nb2 = s.coords + np.array([0, 1, 0])
                    else:
                        raise NotImplemented('currently there are only two possilities')
                    v_s_to_h = unify(np.cross(nb1.coords - s.coords, nb2.coords - s.coords))
                    hcoords = s.coords + v_s_to_h * 1.1
                    hsite = Site('H', hcoords)
                    hsites.append(hsite)
                    hcoords = s.coords + v_s_to_h * -1.1
                    hsite = Site('H', hcoords)
                    hsites.append(hsite)
                else:
                    raise NotImplementedError('currently there are only two possilities')
                hsites_by_joint_siteid[s.properties['siteid']] = hsites
                return hsites_by_joint_siteid

    else:
        if isinstance(joints, list) or isinstance(joints, set):
            if original is None:
                warnings.warn('adding h based on undercoordination at joints')
                nbrmap = c.nbrmap
                hsites_by_joint_siteid = {}
                for i in range(len(c)):
                    s = c[i]
                    if s.properties['siteid'] not in joints:
                        continue
                    hsites = []
                    try:
                        nvalence = _coordination_rule[s.species_string]
                    except KeyError:
                        nvalence = 0
                    if nvalence > len(nbrmap[i]):  # we got a joint
                        nh_to_be_added = nvalence - len(nbrmap[i])
                        if nh_to_be_added == 1 or nh_to_be_added == 3:
                            v_s_to_h = np.zeros(3)
                            for nb_idx in nbrmap[i]:
                                v_nb_to_s = s.coords - c[nb_idx].coords
                                v_s_to_h += v_nb_to_s
                            v_s_to_h = unify(v_s_to_h)
                            hcoords = s.coords + v_s_to_h * 1.1
                            hsite = Site('H', hcoords)
                            hsites.append(hsite)
                        elif nh_to_be_added == 2:
                            if len(nbrmap[i]) == 2:
                                nb1 = c[nbrmap[i][0]]
                                nb2 = c[nbrmap[i][1]]
                            elif len(nbrmap[i]) == 1:
                                nb1 = c[nbrmap[i][0]]
                                nb2 = s.coords + np.array([1, 0, 0])
                                if unify(np.cross(nb1 - s.coords, nb2 - s.coords)) < 1e-5:
                                    nb2 = s.coords + np.array([0, 1, 0])
                            else:
                                raise NotImplemented('currently there are only two possilities')
                            v_s_to_h = unify(np.cross(nb1.coords - s.coords, nb2.coords - s.coords))
                            hcoords = s.coords + v_s_to_h * 1.1
                            hsite = Site('H', hcoords)
                            hsites.append(hsite)
                            hcoords = s.coords + v_s_to_h * -1.1
                            hsite = Site('H', hcoords)
                            hsites.append(hsite)
                        else:
                            raise NotImplementedError('currently there are only two possilities')
                        hsites_by_joint_siteid[s.properties['siteid']] = hsites
                return hsites_by_joint_siteid
        elif isinstance(joints, dict):
            if original is None:
                warnings.warn(
                    'adding h based on fragmenting, but coords of H are calculated seperately as we do not know the original molecule')
                hsites_by_joint_siteid = {}
                nbrmap = c.nbrmap
                for i in range(len(c)):
                    s = c[i]
                    if s.properties['siteid'] not in joints.keys():
                        continue
                    hsites = []
                    nh_to_be_added = len(joints[s.properties['siteid']])
                    if nh_to_be_added == 1 or nh_to_be_added == 3:
                        v_s_to_h = np.zeros(3)
                        for nb_idx in nbrmap[i]:
                            v_nb_to_s = s.coords - c[nb_idx].coords
                            v_s_to_h += v_nb_to_s
                        v_s_to_h = unify(v_s_to_h)
                        hcoords = s.coords + v_s_to_h * 1.1
                        hsite = Site('H', hcoords)
                        hsites.append(hsite)
                    elif nh_to_be_added == 2:
                        if len(nbrmap[i]) == 2:
                            nb1 = c[nbrmap[i][0]]
                            nb2 = c[nbrmap[i][1]]
                        elif len(nbrmap[i]) == 1:
                            nb1 = c[nbrmap[i][0]]
                            nb2 = s.coords + np.array([1, 0, 0])
                            if unify(np.cross(nb1 - s.coords, nb2 - s.coords)) < 1e-5:
                                nb2 = s.coords + np.array([0, 1, 0])
                        else:
                            raise NotImplemented('currently there are only two possilities')
                        v_s_to_h = unify(np.cross(nb1.coords - s.coords, nb2.coords - s.coords))
                        hcoords = s.coords + v_s_to_h * 1.1
                        hsite = Site('H', hcoords)
                        hsites.append(hsite)
                        hcoords = s.coords + v_s_to_h * -1.1
                        hsite = Site('H', hcoords)
                        hsites.append(hsite)
                    else:
                        raise NotImplementedError('currently there are only two possilities')
                    hsites_by_joint_siteid[s.properties['siteid']] = hsites
                return hsites_by_joint_siteid

            else:
                warnings.warn('adding h based on fragmenting, we know the original molecule')
                hsites_by_joint_siteid = {}
                for i in range(len(c)):
                    s = c[i]
                    if s.properties['siteid'] not in joints.keys():
                        continue
                    hsites = []
                    for nb in joints[s.properties['siteid']]:
                        nb_site = original.get_site_byid(nb)
                        v_s_to_h = nb_site.coords - s.coords
                        v_s_to_h = unify(v_s_to_h)
                        hcoords = s.coords + v_s_to_h * 1.1
                        hsite = Site('H', hcoords)
                        hsites.append(hsite)
                    hsites_by_joint_siteid[s.properties['siteid']] = hsites
                return hsites_by_joint_siteid
        else:
            raise TypeError('joints can be a list, set or a dict!')


def conformer_addhmol(c: BasicConformer, joints=None, original: BasicConformer = None):
    hsite_dict = conformer_addh(c, joints, original)

    sites = deepcopy(c.sites)
    sites.sort(key=lambda x:len(x.properties.keys()), reverse=True)
    original_properties = sites[0].properties
    original_site_keys = original_properties.keys()
    hydrogen_assign_keys = ['occu', 'disg', 'iasym', 'imol', 'icell']
    warnings.warn("hydrogen site properties {} will be assigned based on {}".format(hydrogen_assign_keys, original_properties))

    hsites_to_be_added = []
    for k, v in hsite_dict.items():
        hsites_to_be_added += v
    for ihs in range(len(hsites_to_be_added)):
        for k in original_site_keys:
            if k in hydrogen_assign_keys:
                hsites_to_be_added[ihs].properties[k] = original_properties[k]
            elif k == 'siteid':
                hsites_to_be_added[ihs].properties[k] = -(ihs + 1)
            elif k == 'label':
                hsites_to_be_added[ihs].properties[k] = 'addedh'
            else:
                hsites_to_be_added[ihs].properties[k] = 'addh_null'
    sites += hsites_to_be_added
    return Molecule.from_sites(SiteidOperation(sites).sites)


from typing import Union


class FragConformer(BasicConformer):
    _conformer_properties = {}

    def __init__(
            self, sites,
            siteids=False,
            conformer_properties=None,
            graph: Union[FragmentGraph, BackboneGraph, SidechainGraph] = None,
    ):
        super().__init__(sites, siteids)
        if conformer_properties is None:
            self.conformer_properties = {}
        else:
            self.conformer_properties = conformer_properties
        if isinstance(graph, FragmentGraph):
            self.graph = graph
        else:
            raise ConformerError('Missing FragmentGraph!')
        self.joints = self.graph.joints
        self.rings = self.get_rings()

    def get_rings(self):
        rings = nx.minimum_cycle_basis(self.graph.graph)  # technically sssr
        rings = sorted(rings, key=lambda x: len(x))
        rcs = []
        for r in rings:
            ring = RingConformer.from_siteids(r, self.sites, copy=False)
            rcs.append(ring)
        return rcs

    def calculate_conformer_properties(self):
        pass

    @classmethod
    def from_sites(cls, sites, graph=None, siteids=False, conformer_properties=None):
        bc = BasicConformer.from_sites(sites, siteids)
        if graph is None:
            graph = bc.to_graph('siteid', 'FragmentGraph')
        return cls(sites, siteids, conformer_properties, graph)

    @classmethod
    def from_siteids(cls, siteids, sites, graph=None, copy=False, conformer_properties=None, ):
        rs = SiteidOperation(sites).get_sites_byids(siteids, copy)
        return cls.from_sites(rs, siteids=False, conformer_properties=conformer_properties, graph=graph)

    @classmethod
    def from_dict(cls, d: dict):
        graph = FragmentGraph.from_dict(d['graph'])
        prop = d['conformer_properties']
        sites = [Site.from_dict(sd) for sd in d['sites']]
        return cls.from_sites(sites, graph=graph, siteids=False, conformer_properties=prop)

    def as_dict(self):
        d = super().as_dict()
        d['conformer_properties'] = self.conformer_properties
        d['graph'] = self.graph.as_dict()
        return d

    @classmethod
    def from_pmgmol(cls, m: Molecule, siteids=False, conformer_properties=None, graph=None):
        return cls.from_sites(m.sites, siteids, conformer_properties, graph)


class BoneConformer(FragConformer):
    _bone_conformer_properties = {
        'pfit_vo': None,
        'pfit_vp': None,
        'pfit_vq': None,
        'lfiterror': None,
        'pfiterror': None,
    }

    def __init__(
            self, sites,
            siteids=False,
            conformer_properties=None,
            graph: BackboneGraph = None,
    ):
        if conformer_properties is None:
            conformer_properties = self._bone_conformer_properties
        super().__init__(sites, siteids, conformer_properties, graph)
        # self.checkstatus('all assigned', 'unique ids')

        if any(v is None for v in self.conformer_properties.values()) or any(
                k not in self.conformer_properties.keys() for k in self._bone_conformer_properties.keys()):
            self.calculate_conformer_properties()
            prop = self.conformer_properties
            self.pfit_vo = prop['pfit_vo']
            self.pfit_vp = prop['pfit_vp']
            self.pfit_vq = prop['pfit_vq']
            self.lfit_error = prop['lfit_error']
            self.pfit_error = prop['pfit_error']
            self.lfit_ptsmean = prop['lfit_ptsmean']
            self.pfit_ptsmean = prop['pfit_ptsmean']
            # for k, v in self.conformer_properties.items():
            #     self.__setattr__(k, v)

    def calculate_conformer_properties(self):
        if len(self.rings) == 1 or len(self.rings) == 0:
            warnings.warn('W: you are init backbone with one or no ring! fitting params will be set to defaults')
            site_to_geoc = [s.coords - self.geoc for s in self]
            lfit_vp = sorted(site_to_geoc, key=lambda x: norm(x))[-1]
            lfit_vp = unify(lfit_vp)
            lfit_error = 0
            lfit_ptsmean = self.geoc
        else:
            lfit_vp, lfit_ptsmean, lfit_error = Fitter.linear_fit([r.geoc for r in self.rings])
            lfit_vp = unify(lfit_vp)
        # this is the plane fit, vp taken from the projection of lfit_vp
        # the fitted plane is now used as a cart coord sys with origin at ptsmean
        pfit_vo, pfit_ptsmean, pfit_error = Fitter.plane_fit(self.cart_coords)
        pfit_vp = unify(get_proj_point2plane(pfit_ptsmean + lfit_vp, pfit_vo, pfit_ptsmean) - pfit_ptsmean)
        pfit_vq = unify(np.cross(pfit_vo, pfit_vp))
        # pfit_plane_params = get_plane_param(pfit_vo, pfit_ptsmean)
        self.conformer_properties = {
            'pfit_vo': pfit_vo,
            'pfit_vp': pfit_vp,
            'pfit_vq': pfit_vq,
            'lfit_error': lfit_error,
            'pfit_error': pfit_error,
            'lfit_ptsmean': lfit_ptsmean,
            'pfit_ptsmean': pfit_ptsmean
        }

    @property
    def lq(self):
        """
        maxdiff( (s.coord - ref) proj at vq )

        :return: long axis length
        """
        ref = self.rings[0].geoc
        projs = [np.dot(s.coords - ref, self.pfit_vo) for s in self]
        return max(projs) - min(projs)

    @property
    def lp(self):
        """
        maxdiff( (s.coord - ref) proj at vp )

        :return: short axis length
        """
        ref = self.rings[0].geoc
        projs = [np.dot(s.coords - ref, self.pfit_vo) for s in self]
        return max(projs) - min(projs)

    @classmethod
    def from_sites(cls, sites, graph=None, siteids=False, conformer_properties=None):
        bc = BasicConformer.from_sites(sites, siteids)
        if graph is None:
            graph = bc.to_graph('siteid', 'BackboneGraph')
        c = cls(sites, siteids, conformer_properties, graph)
        return c

    @classmethod
    def from_dict(cls, d: dict):
        graph = BackboneGraph.from_dict(d['graph'])
        prop = d['conformer_properties']
        sites = [Site.from_dict(sd) for sd in d['sites']]
        return cls.from_sites(sites, graph=graph, siteids=False, conformer_properties=prop)


class SidechainConformer(FragConformer):
    """
    this can be considered as a special backbone with just one joint
    """
    _sidechain_conformer_properties = {}


    def __init__(
            self, sites,
            siteids=False,
            conformer_properties=None,
            graph: SidechainGraph = None
    ):
        super().__init__(sites, siteids, conformer_properties, graph)
        if conformer_properties is None:
            self.conformer_properties = self._sidechain_conformer_properties

        if len(self.joints.keys()) != 1:
            raise ConformerError('sidechain init with no or >1 joints')
        self.sidechain_joint_siteid = list(self.joints.keys())[0]
        self.sidechain_joint = self.get_site_byid(self.sidechain_joint_siteid)

    @classmethod
    def from_sites(cls, sites, graph=None, siteids=False, conformer_properties=None):
        bc = BasicConformer.from_sites(sites, siteids)
        if graph is None:
            graph = bc.to_graph('siteid', 'SidechainGraph')
        c = cls(sites, siteids, conformer_properties, graph)
        return c

    @classmethod
    def from_dict(cls, d: dict):
        graph = SidechainGraph.from_dict(d['graph'])
        prop = d['conformer_properties']
        sites = [Site.from_dict(sd) for sd in d['sites']]
        return cls.from_sites(sites, graph=graph, siteids=False, conformer_properties=prop)

class MolConformer(BasicConformer):
    rings: List[RingConformer]
    _mol_conformer_properties = {
        # 'smiles': None,
    }

    def __init__(
            self,
            sites,
            siteids=False,
            prop=None,
            graph: MolGraph = None,
            rdmol: Chem.Mol = None,
            smiles: str = None,
            siteid2atomidx: dict = None,
            atomidx2siteid: dict = None,
            backbone: BoneConformer = None,
            sccs: [SidechainConformer] = None,
            backbone_graph: BackboneGraph = None,
            scgs: [SidechainGraph] = None,
            coplane_cutoff=30.0,
            chrombone: FragConformer = None,
            chromsccs: [FragConformer] = None,
            chrombone_graph: FragmentGraph = None,
            chromscgs: [FragmentGraph] = None,
    ):
        super().__init__(sites, siteids)
        # self.checkstatus('all assigned', 'unique ids')

        self.rdmol = rdmol
        self.smiles = smiles
        self.siteid2atomidx = siteid2atomidx
        self.atomidx2siteid = atomidx2siteid
        self.graph = graph
        self.rings = []
        for r in self.graph.rings:
            ring = RingConformer.from_siteids(r, self.sites, copy=False)
            self.rings.append(ring)

        if len(self.rings) < 2:
            self.is_solvent = True
        else:
            self.is_solvent = False

        self.backbone = backbone
        self.backbone_graph = backbone_graph
        self.sccs = sccs
        self.scgs = scgs
        self.chrombone = chrombone
        self.chromsccs = chromsccs
        self.chrombone_graph = chrombone_graph
        self.chromscgs = chromscgs
        self.coplane_cutoff = coplane_cutoff
        if prop is None:
            self.conformer_properties = self._mol_conformer_properties
        else:
            self.conformer_properties = prop

    def as_dict(self):
        d = super().as_dict()
        d['rdmol'] = RdFunc.mol_as_json(self.rdmol)
        d['smiles'] = self.smiles
        d['siteid2atomidx'] = self.siteid2atomidx
        d['atomidx2siteid'] = self.atomidx2siteid
        d['graph'] = self.graph.as_dict()
        d['rings'] = [r.as_dict() for r in self.rings]
        d['coplane_cutoff'] = self.coplane_cutoff
        d['conformer_properties'] = self.conformer_properties
        try:
            d['backbone'] = self.backbone.as_dict()
            d['backbone_graph'] = self.backbone_graph.as_dict()
            d['sccs'] = [sc.as_dict() for sc in self.sccs]
            d['scgs'] = [sg.as_dict() for sg in self.scgs]
        except AttributeError:
            d['backbone'] = None
            d['backbone_graph'] = None
            d['sccs'] = None
            d['scgs'] = None
        try:
            d['chrombone'] = self.chrombone.as_dict()
            d['chrombone_graph'] = self.chrombone_graph.as_dict()
            d['chromsccs'] = [sc.as_dict() for sc in self.chromsccs]
            d['chromscgs'] = [sg.as_dict() for sg in self.chromscgs]
        except AttributeError:
            d['chrombone'] = None
            d['chrombone_graph'] = None
            d['chromsccs'] = None
            d['chromscgs'] = None
        d['conformer_properties'] = self.conformer_properties
        return d

    @classmethod
    def from_dict(cls, d: dict):
        sites = [Site.from_dict(sd) for sd in d['sites']]
        siteids = False
        prop = d['conformer_properties']
        graph = MolGraph.from_dict(d['graph'])
        rdmol = RdFunc.mol_from_json(d['rdmol'])
        smiles = d['smiles']
        siteid2atomidx = d['siteid2atomidx']
        atomidx2siteid = d['atomidx2siteid']
        coplane_cutoff = d['coplane_cutoff']
        backbone = BoneConformer.from_dict(d['backbone'])
        sccs = [SidechainConformer.from_dict(scd) for scd in d['sccs']]
        backbone_graph = BackboneGraph.from_dict(d['backbone_graph'])
        scgs = [SidechainGraph.from_dict(scd) for scd in d['scgs']]
        chrombone = FragConformer.from_dict(d['chrombone'])
        chromsccs = [FragConformer.from_dict(fd) for fd in d['sccs']]
        chrombone_graph = FragmentGraph.from_dict(d['chrombone_graph'])
        chromscgs = [FragmentGraph.from_dict(fg) for fg in d['chromscgs']]

        mc = cls(
            sites,
            siteids,
            prop,
            graph,
            rdmol,
            smiles,
            siteid2atomidx,
            atomidx2siteid,
            backbone,
            sccs,
            backbone_graph,
            scgs,
            coplane_cutoff,
            chrombone,
            chromsccs,
            chrombone_graph,
            chromscgs
        )
        return mc

    def calculate_conformer_properties(self):
        pass

    @staticmethod
    def chrom_partition(mc: BasicConformer, rdmol, atomidx2siteid, molgraph: MolGraph, withhalogen=True):
        if withhalogen:
            cgs = RdFunc.get_conjugate_group(rdmol)
        else:
            cgs = RdFunc.get_conjugate_group_with_halogen(rdmol)
        try:
            chromol, aid_to_newid_chromol = cgs[0]
        except IndexError:
            raise ConformerError('rdkit cannot find a chrom here!')
        aids_in_chromol = list(aid_to_newid_chromol.keys())
        siteids_in_chromol = [atomidx2siteid[aid] for aid in aids_in_chromol]
        chromol_joints, other_joints, chromolsg, sg_components = MolGraph.get_joints_and_subgraph(
            siteids_in_chromol, molgraph.graph)
        cg, fgs = MolGraph.get_chrom_and_frags_from_nxgraph(chromolsg, sg_components)
        chromolc = FragConformer.from_siteids(
            siteids_in_chromol, mc.sites, graph=cg, copy=False
        )  # you need to make sure there is at least a ring here
        fragcs = []
        for fg in fgs:
            fragc = FragConformer.from_siteids(fg.graph.nodes, mc.sites, graph=fg, copy=False, )
            fragcs.append(fragc)
        return chromolc, fragcs, cg, fgs

    @staticmethod
    def geo_partition(bc: BasicConformer, molgraph: MolGraph, coplane_cutoff=30.0):
        try:
            lgfr = []
            for r in molgraph.lgfr:
                ring = RingConformer.from_siteids(r, bc.sites, False)
                lgfr.append(ring)
        except:
            raise ConformerError('cannot get lgfr!')
        avgn1, avgn2 = RingConformer.get_avg_norms(lgfr)

        def coplane_check(ring_siteids, tol=coplane_cutoff):
            ringconformer = RingConformer.from_siteids(ring_siteids, bc.sites, False)
            coplane = ringconformer.iscoplane_with_norm(avgn1, tol, 'degree')
            return coplane

        try:
            bg, scgs = molgraph.partition_to_bone_frags('lgcr', additional_criteria=coplane_check)
        except:
            raise ConformerError('cannot lgcr partition with coplane_cutoff: {}!'.format(coplane_cutoff))

        bone_conformer = BoneConformer.from_siteids(bg.graph.nodes, bc.sites, graph=bg, copy=False)

        sccs = []
        for scg in scgs:
            sc_joint_site_id = list(scg.graph.graph['joints'].keys())[0]
            bc_joint_site_id = scg.graph.graph['joints'][sc_joint_site_id][0]
            bc_joint_site = bc.get_site_byid(bc_joint_site_id)
            # print(bc_joint_site)
            v_sc_position = bc_joint_site.coords - bc.geoc
            sc_position_angle = angle_btw(v_sc_position, bone_conformer.pfit_vp, 'degree')
            scc = SidechainConformer.from_siteids(scg.graph.nodes, bc.sites, copy=False, graph=scg,
                                                  conformer_properties={'sc_position_angle': sc_position_angle,
                                                                        'bc_joint_siteid': bc_joint_site_id})
            sccs.append(scc)
        return bone_conformer, sccs, bg, scgs  # sccs <--> scgs bijection

    @classmethod
    def from_sites(cls, sites, siteids=False, prop=None, coplane_cutoff=30.0, withhalogen=True):
        bc = BasicConformer.from_sites(sites, siteids)
        rdmol, smiles, siteid2atomidx, atomidx2siteid = bc.to_rdmol()
        graph = bc.to_graph('siteid', 'MolGraph')
        try:
            backbone, sccs, backbone_graph, scgs = MolConformer.geo_partition(bc, graph, coplane_cutoff)
        except ConformerError:
            warnings.warn('geo_partition failed!')
            backbone, sccs, backbone_graph, scgs = [None] * 4
        try:
            chrombone, chromsccs, chrombone_graph, chromscgs = MolConformer.chrom_partition(bc, rdmol, atomidx2siteid,
                                                                                            graph, withhalogen)
        except ConformerError:
            warnings.warn('chrom_partition failed!')
            chrombone, chromsccs, chrombone_graph, chromscgs = [None] * 4
        mc = cls(
            bc.sites,
            siteids,
            prop,
            graph,
            rdmol,
            smiles,
            siteid2atomidx,
            atomidx2siteid,
            backbone,
            sccs,
            backbone_graph,
            scgs,
            coplane_cutoff,
            chrombone,
            chromsccs,
            chrombone_graph,
            chromscgs
        )
        return mc

    def to_addhmol(self, bonescheme='geo'):
        """
        add h and return pmg mol

        :return:
        """
        if bonescheme == 'geo':
            backbone_hmol = conformer_addhmol(self.backbone, joints=self.backbone_graph.joints, original=self)
            sccs = self.sccs
            scgs = self.scgs
        elif bonescheme == 'chrom':
            backbone_hmol = conformer_addhmol(self.chrombone, joints=self.chrombone_graph.joints, original=self)
            sccs = self.chromsccs
            scgs = self.chromscgs
        else:
            raise NotImplementedError('addh for {} not implemented'.format(bonescheme))

        schmols = []
        for i in range(len(sccs)):
            scc = sccs[i]
            scg = scgs[i]
            scch = conformer_addhmol(scc, scg.joints, self)
            schmols.append(scch)
        return backbone_hmol, schmols


class ConformerDimer:

    def __init__(self, conformer_ref: MolConformer, conformer_var: MolConformer, label=""):
        """
        basically 2 conformers

        :param conformer_ref:
        :param conformer_var:
        :param str label: mainly used to distinguish

        Attributes:
            vslip: slip vector in cart
            vslipnorm: normalized vslip
            pslip: projection of vslip along vp
            qslip: projection of vslip along vq
            oslip: projection of vslip along vo
            pangle: angle btw vps
            qangle: angle btw vps
            oangle: angle btw vps, always acute
            jmol: jmol draw arrow string in console
        """
        self.conformer_ref = conformer_ref
        self.conformer_var = conformer_var
        self.label = label

        ref_bone = self.conformer_ref.backbone
        var_bone = self.conformer_var.backbone

        self.vslip = var_bone.geoc - ref_bone.geoc
        self.vslipnorm = norm(self.vslip)
        ref_vp_fit = ref_bone.conformer_properties['pfit_vp']
        ref_vq_fit = ref_bone.conformer_properties['pfit_vq']
        ref_vo_fit = ref_bone.conformer_properties['pfit_vo']
        var_vp_fit = var_bone.conformer_properties['pfit_vp']
        var_vq_fit = var_bone.conformer_properties['pfit_vq']
        var_vo_fit = var_bone.conformer_properties['pfit_vo']

        self.pslip = abs(self.vslip @ ref_vp_fit)
        self.qslip = abs(self.vslip @ ref_vq_fit)
        self.oslip = abs(self.vslip @ ref_vo_fit)
        self.pangle = angle_btw(ref_vp_fit, var_vp_fit, output='degree')
        self.qangle = angle_btw(ref_vq_fit, var_vq_fit, output='degree')
        self.oangle = angle_btw(ref_vo_fit, var_vo_fit, output='degree')
        if self.oangle > 90:
            self.oangle = 180 - self.oangle

        self.jmol = "draw arrow {{ {:.4f}, {:.4f}, {:.4f} }} {{ {:.4f}, {:.4f}, {:.4f} }}".format(*ref_bone.geoc,
                                                                                                  *var_bone.geoc)

    def as_dict(self):
        """
        keys are

        vslipnorm, pslip, qslip, oslip, pangle, qangle, oangle, jmol, label, omol_ref, omol_var
        """
        d = {"@module": self.__class__.__module__, "@class": self.__class__.__name__, 'vslipnorm': self.vslipnorm,
             'pslip': self.pslip, 'qslip': self.qslip, 'oslip': self.oslip, 'pangle': self.pangle,
             'qangle': self.qangle, 'oangle': self.oangle, 'jmol': self.jmol, 'label': self.label,
             'ref': self.conformer_ref.as_dict(), 'var': self.conformer_var.as_dict()}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        omol_ref, omol_var, label
        """
        ref = MolConformer.from_dict(d['ref'])
        var = MolConformer.from_dict(d['var'])
        label = d['label']
        return cls(ref, var, label)

    def to_xyz(self, fn, center_label=False):
        sites = self.conformer_ref.sites + self.conformer_var.sites
        if center_label:
            sites += [Site('La', self.conformer_ref.geoc)]
            sites += [Site('La', self.conformer_var.geoc)]
        Molecule.from_sites(sites).to('xyz', fn)

    def to_xyzstring(self, center_label=False):
        sites = self.conformer_ref.sites + self.conformer_var.sites
        if center_label:
            sites += [Site('La', self.conformer_ref.geoc)]
            sites += [Site('La', self.conformer_var.geoc)]
        return str(XYZ(Molecule.from_sites(sites)))

    @property
    def is_identical(self):
        """
        whether two omols are identical, based on norm(vslip) < 1e-5
        """
        return np.linalg.norm(self.vslip) < 1e-5

    @property
    def is_not_identical(self):
        return not self.is_identical

    @property
    def is_bone_close(self, cutoff=6.5):
        """
        use to identify whether this dimer can have minimum wf overlap, ONLY consider bone distance

        this should be called is_not_faraway...

        :param cutoff: minbonedist less than which will be considered close
        :return: bool
        """
        return self.minbonedist < cutoff

    @property
    def is_close(self, cutoff=6.5):
        """
        use to identify whether this dimer can have minimum wf overlap

        this should be called is_not_faraway...

        :param cutoff:
        :return: bool
        """
        return self.mindist < cutoff

    # @property
    # def stack_type(self, cutoff=1.0):
    #     if self.oangle < 30.0 * cutoff:
    #         return 'face_to_face'
    #     else:
    #         return 'face_to_edge'

    @property
    def mindist(self):
        """
        :return: minimum dist between sites on different bones
        """
        distmat = cdist(self.conformer_ref.cart_coords, self.conformer_var.cart_coords)
        return np.min(distmat)

    @property
    def minbonedist(self):
        """
        :return: minimum dist between sites on different bones
        """
        distmat = cdist(self.conformer_ref.backbone.cart_coords, self.conformer_var.backbone.cart_coords)
        return np.min(distmat)

    def plt_bone_overlap(self, algo='convex', output='bone_overlap.eps'):
        """
        plot a 2d graph of how backbones overlap

        using concave or convex hull

        :param algo: concave/convex
        :param output: output filename
        """
        ref_bone = self.conformer_ref.backbone

        ref_p = ref_bone.conformer_properties['pfit_vp']
        ref_q = ref_bone.conformer_properties['pfit_vq']
        ref_o = ref_bone.conformer_properties['pfit_vo']
        origin = self.conformer_ref.backbone.geoc

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.conformer_ref.backbone.sites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.conformer_var.backbone.sites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            raise ConformerError('bone_overlap receive a wrong algo spec')
        x, y = var_hull.exterior.coords.xy
        points = np.array([x, y]).T
        fig, ax = plt.subplots(1)
        polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='r', facecolor='r')
        ax.add_patch(polygon_shape)
        x, y = ref_hull.exterior.coords.xy
        points = np.array([x, y]).T
        polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='b', facecolor='b', label='ref')
        ax.add_patch(polygon_shape)
        ax.axis('auto')
        fig.savefig(output)

    def mol_overlap(self, algo='concave'):
        """
        project var mol onto the plane of ref mol

        :param algo: concave/convex
        :return: area of the overlap, ref omol area, var omol area
        """
        ref_bone = self.conformer_ref.backbone

        ref_p = ref_bone.conformer_properties['pfit_vp']
        ref_q = ref_bone.conformer_properties['pfit_vq']
        ref_o = ref_bone.conformer_properties['pfit_vo']
        origin = self.conformer_ref.backbone.geoc

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.conformer_ref.sites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.conformer_var.sites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            raise ConformerError('overlap receive a wrong algo spec')
        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area))

    def bone_overlap(self, algo='concave'):
        """
        project var backbone onto the plane of ref backbone
        as there's the problem of alpha value in concave hull generation, maybe I should set default to convex, see hull_test

        :param algo: concave/convex
        :return: area of the overlap, ref omol area, var omol area
        """
        ref_bone = self.conformer_ref.backbone

        ref_p = ref_bone.conformer_properties['pfit_vp']
        ref_q = ref_bone.conformer_properties['pfit_vq']
        ref_o = ref_bone.conformer_properties['pfit_vo']
        origin = self.conformer_ref.backbone.geoc

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.conformer_ref.backbone.sites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.conformer_var.backbone.sites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            raise ConformerError('overlap receive a wrong algo spec')

        mindist2d = ref_hull.distance(var_hull)

        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area),
                mindist2d)

    @property
    def sites(self):
        return self.conformer_ref.sites + self.conformer_var.sites


class DimerCollection:
    def __init__(self, dimers: [ConformerDimer]):
        """
        just a list of dimers, they should share the same ref_mol
        """
        self.dimers = dimers

    def get_xyz_string(self, lalabel=False):
        """
        :return: xyz string to be written
        """
        sites = []

        d: ConformerDimer
        for d in self.dimers:
            sites += d.conformer_var.sites
        if lalabel:
            sites += [Site('La', ss.coords, properties=ss.properties) for ss in self.dimers[0].conformer_ref.sites]

        mol = Molecule.from_sites(sites)
        xyz = XYZ(mol)
        return str(xyz)

    def to_xyz(self, fn, lalabel=False):
        """
        :param lalabel:
        :param fn: xyz file name, with extension
        """
        with open(fn, 'w') as f:
            f.write(self.get_xyz_string(lalabel))
