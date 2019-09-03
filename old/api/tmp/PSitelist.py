import numpy as np
from pymatgen.core.structure import Structure
from api.schema.psite import PSite


class PSitelist:

    def __init__(self, psites, pbc, description=''):
        """
        basically Sitelist with PBC
        :param psites: a *list* of psites
        :param pbc:
        :param description:
        """
        self.psites = psites
        self.pbc = pbc
        self.description = description

    def __len__(self):
        return len(self.psites)

    def __contains__(self, psite):
        return psite in self.psites

    def __iter__(self):
        return self.psites.__iter__()

    def __getitem__(self, ind):
        return self.psites[ind]

    def __repr__(self):
        outs = ['PSitelist:']
        for s in self.psites:
            outs.append(s.__repr__())
        return '\n'.join(outs)

    def get_coordmat(self):
        """
        coordinates matrix
        """
        coordmat = np.empty((len(self.psites), 3))
        for i in range(len(self.psites)):
            for j in range(3):
                coordmat[i][j] = self.psites[i].coords[j]
        return coordmat

    def intersection(self, other):
        """
        :param other:
        :return: a list of msites in self that belong to both Sitelists, there is no copy!
        """
        r = []
        for s in self.psites:
            if s in other:
                r.append(s)
        return r

    def get_site(self, siteid):
        for s in self.psites:
            if s.siteid == siteid:
                return s
        return None

    def issubset(self, other):
        return len(self.psites) == len(self.intersection(other))

    @classmethod
    def from_coordmat(cls, mat, names, pbc, ids=None):
        if ids is None:
            ids = np.ones(len(names))
            ids[:] = -1
        ss = []
        for i in range(len(names)):
            ss.append(PSite(names[i], mat[i], pbc, ids[i]))
        return cls(ss, pbc)

    @classmethod
    def from_file(cls, fname):
        structure = Structure.from_file(fname)
        return cls([PSite.from_pymatgen_site(s) for s in structure.sites], structure.lattice.matrix)

    def to(self, ftype, fname):
        pymatgen_sites = [s.to_pymatgen_site() for s in self.psites]
        m = Structure.from_sites(pymatgen_sites)
        m.to(ftype, fname)


