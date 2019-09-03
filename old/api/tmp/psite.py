import sys
import numpy as np
from collections import Hashable
from monty.json import MSONable
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import PeriodicSite
from api.routines.geometry import frac2cart
from api.schema.msite import MSite
from api.schema.element import Element


class PSite(Hashable, MSONable):
    def __init__(self, element_name, coords, pbc, siteid=-1):
        """
        this is the basic unit with PBC
        :param element_name: string
        :param coords: 3*1 list, frac always
        :param pbc: 3*3 list
        :param siteid: -1, this should be set to non-negative only once when a Mol obj is initiated
        """
        self.element = Element(element_name)
        self._siteid = int(siteid)
        self._coords = np.empty(3)
        for i in range(3):
            self._coords[i] = float(coords[i])
        self._pbc = np.empty((3, 3))
        for i in range(3):
            for j in range(3):
                self._pbc[i][j] = float(pbc[i][j])

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
    def siteid(self):
        return self._siteid

    @siteid.setter
    def siteid(self, v):
        try:
            self._siteid = int(v)
        except ValueError:
            sys.exit("id must be int, but you tried to set it as: {}".format(v))

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, vs):
        for i in range(3):
            self._coords[i] = float(vs[i])

    @property
    def pbc(self):
        return self._pbc

    @pbc.setter
    def pbc(self, pbc):
        for i in range(3):
            for j in range(3):
                self._pbc[i][j] = float(pbc[i][j])

    def to_site(self):
        return MSite(self.element.name, frac2cart(self.coords, self.pbc), siteid=self.siteid)

    def __eq__(self, other):
        if self.element == other.element and np.allclose(
                self.coords, other.coords) and self.siteid == other.siteid and np.allclose(self.pbc, other.pbc):
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.element, self.siteid))

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "PSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element, *self.coords, self.siteid)

    def __str__(self):
        return "PSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element, *self.coords, self.siteid)

    def as_dict(self):
        """
        Json-serializable dict representation for Site.
        """
        d = {"element": self.element.name,
             "abc": self.coords.tolist(),
             "pbc": self.pbc.tolist(),
             "siteid": self.siteid,
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Create Site from dict representation
        """
        return cls(d["element"], d["abc"], d['pbc'], d["siteid"])

    @classmethod
    def from_pymatgen_site(cls, pymatgen_site):
        return cls(pymatgen_site.species_string, pymatgen_site.coords, pymatgen_site.lattice.matrix)

    def to_pymatgen_site(self):
        return PeriodicSite(self.element.name, self.coords, Lattice(self.pbc), properties={"siteid": self.siteid})
