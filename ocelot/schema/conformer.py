import itertools
import warnings
from copy import deepcopy

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import rdkit.Chem as Chem
from pymatgen.core.structure import Element
from pymatgen.core.structure import Molecule, Site
from pymatgen.io.xyz import XYZ
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from shapely.geometry import Polygon

from ocelot.routines.conformerparser import ACParser
from ocelot.routines.geometry import angle_btw, Fitter
from ocelot.routines.geometry import coord_transform
from ocelot.routines.geometry import get_proj_point2plane, alpha_shape
from ocelot.routines.geometry import norm
from ocelot.routines.geometry import rotate_along_axis
from ocelot.routines.geometry import rotation_matrix
from ocelot.routines.geometry import unify
from ocelot.schema.graph import BasicGraph, MolGraph

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


class ConformerOperationError(Exception):
    pass


class ConformerInitError(Exception):
    pass


class SiteidError(Exception):
    pass


class DisorderError(Exception):
    pass


class BondInitError(Exception):
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

    def __init__(self, sites):
        self.sites = sites

        allkeys = ['siteid']
        for s in self.sites:
            allkeys += list(s.properties.keys())
        allkeys = list(set(allkeys))

        for i in range(len(self)):
            for k in allkeys:
                try:
                    self[i].properties[k]
                except KeyError:
                    self[i].properties[k] = None

    @property
    def siteids(self):
        return [s.properties['siteid'] for s in self]

    def get_site_byid(self, siteid, multimatch=False):
        """
        :param int siteid:
        :param bool multimatch: if Ture, return a list, otherwise a site
        :return: (a list of) msite obj
        """
        if multimatch:
            sites_matches = []
            for s in self:
                if s.properties['siteid'] == siteid:
                    sites_matches.append(s)
            if len(sites_matches) == 0:
                raise SiteidError('cannot get site by id: {}'.format(siteid))
            return sites_matches
        else:
            for s in self:
                if s.properties['siteid'] == siteid:
                    return s

    def get_sites_byids(self, siteids, copy=False):
        self.checkstatus('all assigned', 'unique ids')
        if not set(siteids).issubset(set(self.siteids)):
            raise SiteidError('siteids is not a subset of sitelist when init conformer')
        rs = []
        for s in self:
            if s.properties['siteid'] in siteids:
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

        if self.siteids[0] == 0:
            s.append('0start')
        else:
            s.append('not0strat')

        if set(self.siteids) == {None}:
            s.append('all none')
        elif None in set(self.siteids):
            s.append('partially none')
        else:
            s.append('all assigned')

        try:
            cri = all(a + 1 == b for a, b in zip(self.siteids, self.siteids[1:]))
            if cri:
                s.append('continuous')
            else:
                s.append('not continuous')
        except TypeError:
            pass

        if len(set(self.siteids)) < len(self.siteids):
            s.append('duplicate ids')
        else:
            s.append('unique ids')

        return s

    def assign_siteid(self, siteids=None):
        """
        :param siteids: default None means siteid = index, otherwise siteid = siteids[i]
        :return:
        """
        if siteids is None:
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

    def __init__(self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None):
        """
        if siteids is None, siteid is set to site.properties['siteid'] or none if site.properties doesnot have a key as 'siteid'

        :param sites:
        :param charge:
        :param spin_multiplicity:
        :param validate_proximity:
        :param siteids: a list of siteids, siteids[index]=siteid of self[index]
        """
        super().__init__(sites)
        self.charge = charge
        self.spin_multiplicity = spin_multiplicity
        if not self.pmgmol.is_ordered:
            raise ConformerInitError('{} is disordered'.format(self.__class__.__name__))

        if validate_proximity:
            if not self.pmgmol.is_valid():
                raise ConformerInitError('sites too close')

        if siteids is not None:
            if len(siteids) == len(self) and all(isinstance(i, int) for i in siteids):
                for i in range(len(self)):
                    self[i].properties['siteid'] = siteids[i]
            else:
                raise SiteidError('siteids are not legit!')

    def __eq__(self, other):
        # using center should be enough for msites in a mol
        if self.pmgmol == other.pmgmol:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

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
            istring = sites[i].species_string
            irad = Element(istring).atomic_radius
            if irad is None:
                continue
            for j in range(i + 1, numsites):
                jstring = sites[j].species_string
                jrad = Element(jstring).atomic_radius
                if jrad is None:
                    continue
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
        self.checkstatus('all assigned', 'unique ids')
        distmat = self.distmat
        siteid_to_disti = self.siteid2index
        mat = np.zeros((max(self.siteids) + 10, max(self.siteids) + 10), dtype=bool)
        for sidi in self.siteids:
            i = siteid_to_disti[sidi]
            istring = self[i].species_string
            irad = Element(istring).atomic_radius
            if irad is None:
                continue
            for sidj in self.siteids:
                j = siteid_to_disti[sidj]
                jstring = self[j].species_string
                jrad = Element(jstring).atomic_radius
                if jrad is None:
                    continue
                distance = distmat[i][j]
                cutoff = (irad + jrad) * co
                if 1e-5 < distance < cutoff:
                    mat[sidi][sidj] = 1
                    mat[sidj][sidi] = 1
        return mat

    def to_graph(self, nodename='siteid', graphtype='molgraph'):
        """
        get a Graph, default siteid --> nodename

        otherwise nodename is assigned based on the order in self.sites
        """
        bondmat_by_siteid = self.get_bondmat_based_on_siteid(co=1.3)
        g = nx.Graph()
        if nodename == 'siteid':
            for s in self:
                g.add_node(s.properties['siteid'], symbol=s.species_string)
                for ss in self:
                    g.add_node(ss.properties['siteid'], symbol=ss.species_string)
                    if s == ss:
                        continue
                    else:
                        if bondmat_by_siteid[s.properties['siteid']][ss.properties['siteid']]:
                            g.add_edge(s.properties['siteid'], ss.properties['siteid'])
        if graphtype == 'molgraph':
            return MolGraph(g)
        return BasicGraph(g)

    def to_rdmol(self, charge=0, sani=True, charged_fragments=True, force_single=False, expliciths=True):
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
        ac = self.bondmat
        ap = ACParser(ac, charge, self.atomic_numbers, sani=sani)
        rdmol, smiles = ap.parse(charged_fragments=charged_fragments, force_single=force_single, expliciths=expliciths)
        rdmol.AddConformer(conf)
        return rdmol, smiles, siteid2atomidx, atomidx2siteid

    def to_file(self, fmt, filename):
        self.pmgmol.to(fmt, filename)

    @classmethod
    def from_sites(cls, sites, charge=0, spin_multiplicity=None,
                   validate_proximity=True, siteidcheck=('all assigned', 'unique ids'), siteids=None):
        """
        we use this as the parent constructor

        :param sites:
        :param charge:
        :param spin_multiplicity:
        :param validate_proximity:
        :param siteidcheck: if the check fails, assign siteid as the order in m
        :param siteids:
        :return:
        """
        so = SiteidOperation(sites)
        try:
            so.checkstatus(*siteidcheck)
        except SiteidError:
            warnings.warn('siteid check in {} failed, reassign siteid...'.format(cls.__class__.__name__))
            so.assign_siteid()
        # print(cls.__class__.__name__, so.sites)
        if siteids is None:
            try:
                idsfromsites = [s.properties['siteid'] for s in sites]
            except KeyError:
                idsfromsites = None
            siteids = idsfromsites
        return cls(so.sites, charge, spin_multiplicity, validate_proximity, siteids)

    @classmethod
    def from_pmgmol(cls, m: Molecule, validate_proximity=True, siteidcheck=('all assigned', 'unique ids'),
                    siteids=None):
        return cls.from_sites(m.sites, int(m.charge), m.spin_multiplicity, validate_proximity, siteidcheck, siteids)

    @classmethod
    def from_siteids(cls, siteids, sites, copy=False, charge=0, spin_multiplicity=None, validate_proximity=True,
                     siteidcheck=('all assigned', 'unique ids')):
        """
        build up a ring based on a list of siteid and a list of sites

        the sites should have been assigned siteids

        :param copy:
        :param charge:
        :param spin_multiplicity:
        :param validate_proximity:
        :param siteidcheck:
        :param siteids: a list of siteid
        :param sites: all sites (in a mol)
        """
        rs = SiteidOperation(sites).get_sites_byids(siteids, copy)
        return cls.from_sites(rs, charge, spin_multiplicity, validate_proximity, siteidcheck)

    @classmethod
    def from_file(cls, filename):

        m = Molecule.from_file(filename)
        c = cls.from_pmgmol(m)
        warnings.warn(
            'conformer built from file, siteids are assigned by the order of entries in file: {}'.format(filename))
        return c

    @classmethod
    def from_str(cls, input_string, fmt):
        m = Molecule.from_str(input_string, fmt)
        conformer = cls.from_pmgmol(m)
        warnings.warn('conformer built from string, siteids are assigned by the order of entries in the string')
        return conformer

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "sites": [s.as_dict() for s in self],
            "charge": self.charge,
            "spin_multiplicity": self.spin_multiplicity,
        }
        return d

    @classmethod
    def from_dict(cls, d: dict):
        sites = [Site.from_dict(sd) for sd in d['sites']]
        return cls.from_sites(sites, d['charge'], d['spin_multiplicity'])


def get_rings_from_conformer(conformer: BasicConformer, scheme='all'):
    """
    :param conformer:
    :param scheme: 'all', 'lgfr', 'lgcr'
    :return:a list of :class:`RingConformer`
    """
    conformer.checkstatus('all assigned', 'unique ids')
    molgraph = conformer.to_graph('siteid', 'molgraph')
    rings = []
    if scheme == 'all':
        for r in molgraph.rings:
            ring = RingConformer.from_siteids(r, conformer.sites, copy=False)
            rings.append(ring)
    elif scheme == 'lgfr':
        for r in molgraph.lgfr:
            ring = RingConformer.from_siteids(r, conformer.sites, copy=False)
            rings.append(ring)
    elif scheme == 'lgcr':
        for r in molgraph.lgcr:
            ring = RingConformer.from_siteids(r, conformer.sites, copy=False)
            rings.append(ring)
    else:
        raise NotImplementedError('I do not understand the scheme {}!'.format(scheme))
    return rings


class BondConformer(BasicConformer):

    def __init__(self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None):
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids)
        # self.checkstatus('all assigned', 'unique ids')
        if len(self) != 2:
            raise ConformerInitError('only use two sites to create a bond! you are using {}'.format(len(self)))
        self.a = self[0]
        self.b = self[1]

    @property
    def length(self):
        return norm(self.a.coords - self.b.coords)


class RingConformer(BasicConformer):
    def __init__(self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None):
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids)
        # self.checkstatus('all assigned', 'unique ids')
        self.n_member = len(self)
        if self.n_member < 3:
            raise ConformerInitError('you are initializing a ring with less than 3 sites!')
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
        d['bonds'] = [b.as_dict() for b in self.bonds_in_ring]
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
                        nb_site = original.get_site_byid(nb, multimatch=False)
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
    for k, v in hsite_dict.items():
        sites += v
    return Molecule.from_sites(SiteidOperation(sites).sites)


class FragConformer(BasicConformer):
    _conformer_properties = {}

    def __init__(
            self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None,
            conformer_properties=None,
            rings=None,
            joints=None,
    ):
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids)
        # self.checkstatus('all assigned', 'unique ids')
        if conformer_properties is None:
            self.conformer_properties = {}
        else:
            self.conformer_properties = conformer_properties

        self.joints = joints
        if self.joints is None:
            # raise ConformerInitError('BoneConformer init without joints dict defined')
            warnings.warn('{} init without joints'.format(self.__class__.__name__))

        if rings is None:
            self.rings = get_rings_from_conformer(self, 'all')
        else:
            self.rings = rings

    def calculate_conformer_properties(self):
        pass

    @classmethod
    def from_sites(
            cls, sites, charge=0, spin_multiplicity=None, validate_proximity=True,
            siteidcheck=('all assigned', 'unique ids'), siteids=None,
            conformer_properties=None,
            rings=None,
            joints=None,
    ):
        """

        :param conformer_properties:
        :param rings:
        :param joints:
        :param sites:
        :param charge:
        :param spin_multiplicity:
        :param validate_proximity:
        :param siteidcheck: if the check fails, assign siteid as the order in m
        :param siteids:
        :return:
        """
        if siteids is None:
            try:
                idsfromsites = [s.properties['siteid'] for s in sites]
            except KeyError:
                idsfromsites = None
            siteids = idsfromsites
        c = cls(sites, charge, spin_multiplicity, validate_proximity, siteids, conformer_properties, rings, joints)
        try:
            c.checkstatus(*siteidcheck)
        except SiteidError:
            warnings.warn('siteid check in {} failed, reassign siteid...'.format(cls.__class__.__name__))
            c.assign_siteid()
        return c

    @classmethod
    def from_rings(
            cls, rings, charge=0, spin_multiplicity=None, validate_proximity=True,
            siteidcheck=('all assigned', 'unique ids'), siteids=None,
            conformer_properties=None,
            joints=None,
    ):
        sites = []
        for r in rings:
            for rs in r:
                if rs not in sites:
                    sites.append(rs)
        c = cls(sites, charge, spin_multiplicity, validate_proximity, siteids, conformer_properties, rings, joints)
        try:
            c.checkstatus(*siteidcheck)
        except SiteidError:
            warnings.warn('siteid check in {} failed, reassign siteid...'.format(cls.__class__.__name__))
            c.assign_siteid()
        return c

    @classmethod
    def from_basic(
            cls, bc: BasicConformer, validate_proximity=True,
            conformer_properties=None,
            rings=None,
            joints=None,
    ):
        return cls.from_sites(bc.sites, bc.charge, bc.spin_multiplicity, validate_proximity, siteids=bc.siteids,
                              conformer_properties=conformer_properties, rings=rings, joints=joints)

    @classmethod
    def from_siteids(cls, siteids, sites, copy=False, charge=0, spin_multiplicity=None, validate_proximity=True,
                     siteidcheck=('all assigned', 'unique ids'),
                     conformer_properties=None,
                     rings=None,
                     joints=None,
                     ):
        rs = SiteidOperation(sites).get_sites_byids(siteids, copy)
        return cls.from_sites(
            rs, charge, spin_multiplicity, validate_proximity,
            ('all assigned', 'unique ids'), siteids=None,
            conformer_properties=conformer_properties,
            rings=rings,
            joints=joints,
        )

    @classmethod
    def from_dict(cls, d: dict):
        rings = [RingConformer.from_dict(rd) for rd in d['rings']]
        joints = d['joints']
        prop = d['conformer_properties']
        bc = BasicConformer.from_dict(d)
        return cls.from_basic(bc, True, prop, rings, joints)

    @classmethod
    def from_pmgmol(cls, m: Molecule, validate_proximity=True, siteidcheck=('all assigned', 'unique ids'), siteids=None,
                    conformer_properties=None, rings=None, joints=None, ):
        return cls.from_sites(m.sites, int(m.charge), m.spin_multiplicity, validate_proximity, siteidcheck, siteids,
                              conformer_properties, rings, joints)


class BoneConformer(FragConformer):
    _bone_conformer_properties = {
        'pfit_vo': None,
        'pfit_vp': None,
        'pfit_vq': None,
        'lfiterror': None,
        'pfiterror': None,
    }

    def __init__(
            self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None,
            conformer_properties=None,
            rings=None,
            joints=None,
    ):
        if conformer_properties is None:
            conformer_properties = self._bone_conformer_properties
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids, conformer_properties, rings,
                         joints)
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
        if len(self.rings) == 1:
            warnings.warn('W: you are init backbone with one ring! fitting params will be set to defaults')
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

    def as_dict(self):
        d = super().as_dict()
        d['rings'] = [r.as_dict() for r in self.rings]
        d['joints'] = self.joints
        d['conformer_properties'] = self.conformer_properties
        return d


class SidechainConformer(FragConformer):
    """
    this can be considered as a special backbone with just one joint
    """
    _sidechain_conformer_properties = {}

    def __init__(
            self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None,
            conformer_properties=None,
            rings=None,
            joints=None,
    ):
        if conformer_properties is None:
            conformer_properties = self._sidechain_conformer_properties
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids, conformer_properties, rings,
                         joints)

        try:
            if len(self.joints.keys()) != 1:
                raise ConformerInitError('sidechain init with no or >1 joints')
            self.sidechain_joint_siteid = list(self.joints.keys())[0]
            self.sidechain_joint = self.get_site_byid(self.sidechain_joint_siteid)
        except TypeError:
            warnings.warn('sidechain joints is not a dict!')
            self.sidechain_joint_siteid = None
            self.sidechain_joint = None


class MolConformer(BasicConformer):
    _mol_conformer_properties = {
        # 'smiles': None,
    }

    def __init__(
            self, sites, charge=0, spin_multiplicity=1, validate_proximity=True, siteids=None, prop=None
    ):
        super().__init__(sites, charge, spin_multiplicity, validate_proximity, siteids)
        # self.checkstatus('all assigned', 'unique ids')

        self.rings = get_rings_from_conformer(self, 'all')
        if len(self.rings) <= 1:
            self.is_solvent = True
        else:
            self.is_solvent = False

        self.backbone, self.sccs, self.backbone_graph, self.scgs = self.partition(coplane_cutoff=30.0)
        if prop is None:
            self.conformer_properties = self._mol_conformer_properties
        else:
            self.conformer_properties = prop

    def calculate_conformer_properties(self):
        pass

    def partition(self, coplane_cutoff=40.0, scheme=None):
        molgraph = self.to_graph('siteid', 'molgraph')
        lgfr = get_rings_from_conformer(self, 'lgfr')
        avgn1, avgn2 = RingConformer.get_avg_norms(lgfr)

        def coplane_check(ring_siteids, tol=coplane_cutoff):
            ringconformer = RingConformer.from_siteids(ring_siteids, self.sites, False)
            coplane = ringconformer.iscoplane_with_norm(avgn1, tol, 'degree')
            return coplane

        if scheme is None and isinstance(coplane_cutoff, float):

            bg, scgs = molgraph.partition_to_bone_frags('lgcr', additional_criteria=coplane_check)
        else:
            bg, scgs = molgraph.partition_to_bone_frags(scheme)

        bone_conformer = BoneConformer.from_siteids(
            bg.graph.nodes, self.sites, copy=False, joints=bg.graph.graph['joints'], rings=None
        )

        sccs = []
        for scg in scgs:
            sc_joint_site_id = list(scg.graph.graph['joints'].keys())[0]
            bc_joint_site_id = scg.graph.graph['joints'][sc_joint_site_id][0]
            bc_joint_site = self.get_site_byid(bc_joint_site_id)
            # print(bc_joint_site)
            v_sc_position = bc_joint_site.coords - self.geoc
            sc_position_angle = angle_btw(v_sc_position, bone_conformer.pfit_vp, 'degree')
            scc = SidechainConformer.from_siteids(scg.graph.nodes, self.sites, copy=False,
                                                  joints=scg.graph.graph['joints'],
                                                  rings=None,
                                                  conformer_properties={'sc_position_angle': sc_position_angle,
                                                                        'bc_joint_siteid': bc_joint_site_id})
            sccs.append(scc)
        return bone_conformer, sccs, bg, scgs  # sccs <--> scgs bijection

    def to_addhmol(self):
        """
        add h and return pmg mol

        :return:
        """
        backbone_hmol = conformer_addhmol(self.backbone, joints=self.backbone_graph.joints, original=self)
        schmols = []
        for i in range(len(self.sccs)):
            scc = self.sccs[i]
            scg = self.scgs[i]
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

    def plt_bone_overlap(self, algo='concave', output='bone_overlap.eps'):
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
            raise ConformerOperationError('bone_overlap receive a wrong algo spec')
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
            raise ConformerOperationError('overlap receive a wrong algo spec')
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
            raise ConformerOperationError('overlap receive a wrong algo spec')

        mindist2d = ref_hull.distance(var_hull)

        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area),
                mindist2d)

    @property
    def sites(self):
        return self.conformer_ref.sites + self.conformer_var.sites


class DimerCollection:
    def __init__(self, dimers):
        """
        just a list of dimers, they should share the same ref_mol
        """
        self.dimers = dimers

    def get_xyz_string(self, lalabel=False):
        """
        :return: xyz string to be written
        """
        sites = []

        for d in self.dimers:
            sites += d.sites
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
