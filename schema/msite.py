from schema.element import Element
import numpy as np
from collections import Hashable
from pymatgen.core.sites import Site
import sys
import warnings


class MSite(Hashable):
    def __init__(self, element_name, coords, siteid=-1):
        """
        :param element_name: string
        :param coords: 3*1 list, cart always
        :param siteid: -1, this should be set to non-negative only once when a Mol obj is initiated
        """
        self.element = Element(element_name)
        self._siteid = int(siteid)
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
    def siteid(self):
        return self._siteid

    @siteid.setter
    def siteid(self, v):
        try:
            self._siteid = int(v)
        except ValueError:
            sys.exit("E: id must be int, but you tried to set it as: {}".format(v))

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
        if self.element == other.element and np.allclose(self.coords, other.coords) and self.siteid == other.siteid:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash((self.element, self.x, self.y, self.z))

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "MSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element.name, *self.coords, self.siteid)

    def __str__(self):
        return "MSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element.name, *self.coords, self.siteid)

    def as_dict(self):
        """
        Json-serializable dict representation for Site.
        """
        d = {"element": self.element.name,
             "xyz": self.coords.tolist(),
             "siteid": self.siteid,
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        Create Site from dict representation
        """
        try:
            return cls(d["element"], d["xyz"], d["siteid"])
        except KeyError:
            warnings.warn('W: cannot init MSite from dict as no required keys')
            return None

    @classmethod
    def from_pymatgen_site(cls, pymatgen_site):
        return cls(pymatgen_site.species_string, pymatgen_site.coords)

    def to_pymatgen_site(self):
        # return Site(self.element.name, self.coords, properties={"siteid": self.siteid})
        return Site(str(self.element.name), self.coords)
