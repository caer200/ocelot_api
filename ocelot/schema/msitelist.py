import numpy as np
import warnings
import itertools
from copy import deepcopy
from ocelot.routines.geometry import rotate_along_axis, rotation_matrix, coord_transform, norm
from scipy.spatial.distance import pdist, squareform, cdist
from pymatgen.core.structure import Molecule
from ocelot.routines.xyz2mol import xyz2mol  # https://github.com/jensengroup/xyz2mol
from ocelot.schema.msite import MSite
from ocelot.schema.element import Element
import rdkit.Chem as Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import rdkit.Chem.Descriptors3D as descriptors3d


class MSitelist:
    def __init__(self, msites):
        """
        this is the parent obj for backbone, sidegroups, omol, ring, bond.

        duplicate msite will be discarded

        :param list msites: a *list* of msites
        """
        unique_sites = []
        for ms in msites:
            if ms not in unique_sites:
                unique_sites.append(ms)
        self.msites = unique_sites

    def __len__(self):
        return len(self.msites)

    def __contains__(self, msite):  # TODO necessary?
        for s in self.msites:
            if np.allclose(msite.coords, s.coords) and s.element == msite.element:
                return True
        return False

    def __iter__(self):
        return self.msites.__iter__()

    def __getitem__(self, ind):
        return self.msites[ind]

    def __repr__(self):
        outs = [self.__class__.__name__ + ': ']
        for s in self.msites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    @property
    def rdkit_mol(self):
        """
        return a rdkit mol object using the code from

        https://github.com/jensengroup/xyz2mol

        #TODO lic issue
        :return:
        """
        atomic_number_list = [Element.atomic_numbers[name] for name in self.elementarray]
        mol = xyz2mol(atomicNumList=atomic_number_list, charge=0, xyz_coordinates=self.coordmat.tolist(),
                      charged_fragments=False, quick=True)
        return mol

    @property
    def shape_descriptors(self):
        """
        using methods from rdkit

        https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors3D.html

        Asphericity: 0 corresponds to spherical top molecules and 1 to linear molecules, For prolate (cigar-sh aped) molecules, ~ 0.25, whereas for oblate (disk-shaped) molecules, ~ 1.

        Eccentricity: 0 corresponds to spherical top molecules and 1 to linear molecules.

        RadiusOfGyration: a size descriptor based on the distribution of atomic masses in a molecule, a measure of molecular compactness for long-chain molecules and, specifically, small values are obtained when most o f the atoms are close to the center of mass

        SpherocityIndex: spherosity index varies from zero for flat molecules, such as benzene, to unity for totally spherical molecules

        :key: Asphericity, Eccentricity, InertialShapeFactor, NPR1, NPR2, PMI1, PMI2, PMI3, RadiusOfGyration, SpherocityIndex

        :return: dict
        """
        mol = self.rdkit_mol
        d = dict(
            Asphericity=descriptors3d.Asphericity(mol),
            Eccentricity=descriptors3d.Eccentricity(mol),
            InertialShapeFactor=descriptors3d.InertialShapeFactor(mol),
            NPR1=descriptors3d.NPR1(mol),
            NPR2=descriptors3d.NPR2(mol),
            PMI1=descriptors3d.PMI1(mol),
            PMI2=descriptors3d.PMI2(mol),
            PMI3=descriptors3d.PMI3(mol),
            RadiusOfGyration=descriptors3d.RadiusOfGyration(mol),
            SpherocityIndex=descriptors3d.SpherocityIndex(mol),
        )
        return d

    def fp_similarity(self, other, metric='Tanimoto'):
        """
        use RDK fingerprint similarity based on different metrics

        TODO add args to customize RDKfp, see https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint

        see Landrum2012 for more details

        :param str metric: "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet",
        :return:
        """
        mol = self.rdkit_mol
        other_mol = other.rdkit_mol
        fps = [FingerprintMols.FingerprintMol(mol), FingerprintMols.FingerprintMol(other_mol)]
        for func in DataStructs.similarityFunctions:
            if func[0] == metric:
                metric_function = func[1]
                return DataStructs.FingerprintSimilarity(fps[0], fps[1], metric=metric_function)
        return None

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

    @staticmethod
    def get_nbrmap(bmat):
        """
        :param np.ndarray bmat: bool bond matrix, i is not bonded to itself
        :return: nbmap[i] is the index list of i's neighbors
        """
        # TODO profile
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
        other_sites -- self_border_site ----- distmin ----- other_border_site --- other_sites

        :param MSitelist other:
        :return: self_border_site, other_border_site, distmin
        """
        distmat = cdist(self.coordmat, other.coordmat)
        minid = np.unravel_index(np.argmin(distmat, axis=None), distmat.shape)
        self_border_site = distmat[minid[0]]
        other_border_site = distmat[minid[1]]
        distmin = distmat[minid]
        return self_border_site, other_border_site, distmin

    def get_site_by_coords(self, coords):
        for s in self.msites:
            if np.allclose(s.coords, coords):
                return s
        return None

    @staticmethod
    def get_bondmat(msites, distmat, co=1.3):
        """
        Bij = whether there is a bond between si and sj
        i is NOT bonded with itself

        :param co: coefficient for cutoff, default 1.3, based on covalent rad
        :return: bool matrix
        """
        numsites = len(msites)
        mat = np.ones((numsites, numsites), dtype=bool)
        for i in range(numsites):
            mat[i][i] = False
            for j in range(i + 1, numsites):
                if distmat[i][j] > (msites[i].element.covrad + msites[j].element.covrad) * co:
                    mat[i][j] = False
                    mat[j][i] = False
        return mat

    @property
    def elementarray(self):
        """
        :return: a nx1 string np array
        """
        arr = []
        for i in range(len(self)):
            arr.append(self[i].element.name)
        return np.array(arr)

    @property
    def geoc(self):
        """
        geometric center

        :return: 3x1 np array
        """
        v = np.zeros(3)
        for s in self.msites:
            v += s.coords
        return v / len(self.msites)

    @property
    def canonical_smiles(self):
        """
        :return: smiles string, using rdkit Chem.MolToSmiles, isomericSmiles=False
        """
        return Chem.MolToSmiles(self.rdkit_mol, isomericSmiles=False)

    @property
    def volume_slow(self, boxdensity=0.2):
        """
        http://wiki.bkslab.org/index.php/Calculate_volume_of_the_binding_site_and_molecules

        First, Lay a grid over the spheres.

        Count the number or points contained in the spheres (Ns).

        Count the number of points in the grid box (Ng).

        Calculate the volume of the grid box (Vb).

        don't use this as it's slow, use `volume` as implemented in rdkit

        :return: volume in \AA^3
        """
        mat = np.empty((len(self.msites), 4))
        for i in range(len(self.msites)):
            for j in range(3):
                mat[i][j] = self.msites[i].coords[j]
            mat[i][3] = self.msites[i].element.covrad

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


    @property
    def volume(self):
        """
        AllChem.ComputeMolVolume in rdkit

        :return: volume in \AA^3
        """
        return AllChem.ComputeMolVolume(self.rdkit_mol, confId=-1, gridSpacing=0.2, boxMargin=2.0)

    def intersection(self, other, copy=False):
        """
        :param MSitelist other:
        :param bool copy: whether generate a list of deepcopied sites or not
        :return: a list of msites in self that belong to both Sitelists
        """
        r = []
        if copy:
            for s in self.msites:
                if s in other:
                    r.append(deepcopy(s))
        else:
            for s in self.msites:
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
            for s in self.msites:
                if s.siteid == siteid:
                    sites_matches.append(s)
            if len(sites_matches) == 0:
                warnings.warn('W: cannot get site by id: {}'.format(siteid))
                return None
            return sites_matches
        else:
            for s in self.msites:
                if s.siteid == siteid:
                    return s

    def issubset(self, other):
        """
        based on msite.coords and msite.element, see __eq__ in msite

        :param MSitelist other:
        :rtype: bool
        """
        return len(self.msites) == len(self.intersection(other))

    def as_dict(self):
        """
        keys are

        msites, can, volume

        :return: a dict representation
        """
        d = {
            "msites": [s.as_dict() for s in self.msites],
            "can": self.canonical_smiles,
            "volume": self.volume,
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__
        }
        return d

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        msites

        :param dict d: d['msites'] is a list of msites' dicts
        """
        return cls([MSite.from_dict(entry) for entry in d["msites"]])

    @classmethod
    def from_coordmat_and_elementlst_and_siteidlst(cls, mat, elementlst, ids=None):
        """
        :param np.ndarray mat: coords matrix nx3
        :param list elementlst: a list of string as element_name
        :param list ids: siteid list
        """
        if ids is None:
            ids = np.ones(len(elementlst))
            ids[:] = -1
        if not len(mat) == len(elementlst) == len(ids):
            warnings.warn('W: length of the args are inconsistent')
            return None
        ss = []
        for i in range(len(elementlst)):
            ss.append(MSite(elementlst[i], mat[i], ids[i]))
        return cls(ss)

    @classmethod
    def from_file(cls, fname):
        """
        using pmg molecule

        :param str fname:
        """
        m = Molecule.from_file(fname)
        return cls([MSite.from_pymatgen_site(s) for s in m.sites])

    def to(self, ftype, fname):
        """
        using pmg molecule

        :param str ftype: 'xyz', 'gjf', etc.
        :param str fname: with extension
        """
        pymatgen_sites = [s.to_pymatgen_site() for s in self.msites]
        m = Molecule.from_sites(pymatgen_sites)
        m.to(ftype, fname)

    @property
    def pmgmol(self):
        """
        notice siteid info lost here

        :return: pymatgen molecule obj
        """
        pymatgen_sites = [s.to_pymatgen_site() for s in self.msites]
        return Molecule.from_sites(pymatgen_sites)

    def rotate_along(self, theta, end1, end2, unit='degree'):
        """
        rotate the vectors defined by (site.coords - end1) along (end2 - end1)

        end2 ------ site

        \|

        \|

        \|

        end1

        notice the coords are changed in-situ

        :param theta: angle of rotation
        :param end1: 3x1 float list/array
        :param end2: 3x1 float list/array
        :param str unit: degree/radian
        """
        end1 = np.array(end1)
        end2 = np.array(end2)
        for s in self.msites:
            v = s.coords - end1
            s.coords = end1 + rotate_along_axis(v, end2 - end1, theta, thetaunit=unit)

    def rotate_along_de_matrix(self, theta, end1, end2, unit='degree'):
        """
        rotation matrix for :func:`~msitelist.MSitelist.rotate_along`
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
        for s in self.msites:
            v = s.coords - end1
            s.coords = end1 + np.dot(matrix, v)

    @classmethod
    def orient(cls, msitelist, pqo):
        """
        p, q, o are 3 orthonormal vectors in x, y, z basis

        basis transformation from xyz to pqo

        :param pqo: 3x3 array
        """
        pqo = np.array(pqo)
        p, q, o = pqo
        nmsites = []
        for s in msitelist.msites:
            ns = deepcopy(s)
            ns.coords = coord_transform(p, q, o, s.coords - msitelist.geoc)
            nmsites.append(ns)
        return cls(nmsites)
