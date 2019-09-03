import sys
import numpy as np
from collections import Hashable
from monty.json import MSONable
from api.routines.geometry import cart2frac
from api.schema.element import Element
from api.schema.psite import PSite
<<<<<<< HEAD
from pymatgen.core.sites import Site
=======
>>>>>>> 88b3e96e2d37f403133c7972e5e02f735a27c716


class MSite(Hashable, MSONable):
    def __init__(self, element_name, coords, siteid=-1):
        """
        this is the basic structure unit
        :param element_name: string
        :param coords: 3*1 list, cart always
        :param siteid: -1, this should be set to non-negative only once when a Mol obj is initiated
        """
        self.element = Element(element_name)
        self._siteid = int(siteid)
        self._coords = np.empty(3)
        for i in range(3):
            self._coords[i] = float(coords[i])

    def to_psite(self, pbc):
        return PSite(self.element.name, cart2frac(self.coords, pbc), pbc, siteid=self.siteid)

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
        return hash(self.element)

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element, *self.coords, self.siteid)

    def __str__(self):
        return "Site: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element, *self.coords, self.siteid)

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
        return cls(d["element"], d["xyz"], d["siteid"])

    @classmethod
    def from_pymatgen_site(cls, pymatgen_site):
        return cls(pymatgen_site.species_string, pymatgen_site.coords)
<<<<<<< HEAD

    def to_pymatgen_site(self):
        return Site(self.element.name, self.coords, properties={"siteid": self.siteid})
=======
>>>>>>> 88b3e96e2d37f403133c7972e5e02f735a27c716
