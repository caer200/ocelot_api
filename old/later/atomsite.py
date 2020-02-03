import numpy as np
import sys
import warnings
from pymatgen.core.periodic_table import Element
from pymatgen.core.sites import Site


class AtomSite:
    def __init__(self, element_name, coords=(0, 0, 0), siteid=None):
        """
        :param str element_name: case sensitive
        :param np.ndarray coords: 3*1, cart always
        :param siteid: default None, this should be set as an int when init a molecule
        """
        self.element = Element(element_name)
        self.symbol = element_name
        self.Z = self.element.Z
        self._siteid = siteid
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
        # TODO profile
        return np.linalg.norm(self.coords - other.coords)

    def __eq__(self, other, tol=1e-5):
        """
        two osites are equal if they have the same element and are close in cart system

        this is useful for set or unique list

        :param AtomSite other:
        """
        if self.element.symbol == other.element.symbol and np.allclose(self.coords, other.coords, atol=tol):
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        rarely used, just in case
        """
        return hash((self.element.symbol, self.x, self.y, self.z, self.siteid))

    def __repr__(self):
        return "AtomSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element.symbol, *self.coords, self.siteid)

    def __str__(self):
        return "AtomSite: {} ({:.4f}, {:.4f}, {:.4f}) siteid: {}".format(
            self.element.symbol, *self.coords, self.siteid)

    def as_dict(self):
        """
        keys are

        element, xyz, siteid

        :return: json-serializable representation for AtomSite.
        :rtype: dict
        """
        d = {"element": self.element.symbol,
             "xyz": self.coords.tolist(),
             "siteid": self.siteid,
             "@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        return d

    @classmethod
    def from_dict(cls, d):
        """
        create AtomSite from dict representation, return None if failed

        keys are

        element, xyz, siteid

        :param dict d: representation of a AtomSite
        :return: AtomSite, or None if keys not found
        """
        try:
            return cls(d["element"], d["xyz"], d["siteid"])
        except KeyError:
            warnings.warn('W: cannot init AtomSite from dict as no required keys')
            return None

    @classmethod
    def from_pymatgen_site(cls, pymatgen_site):
        """
        create AtomSite from a pmg Site obj, the coords for that pmg site should be in cart

        if 'siteid' is not a key in pmg site property, siteid is default to None

        :param pmg_site pymatgen_site: the pmg site that will be converted
        :return: AtomSite
        """
        siteid = None
        if 'siteid' in pymatgen_site.properties.keys():
            siteid = pymatgen_site.properties['siteid']
        symbol = pymatgen_site.species.elements[0].symbol
        return cls(symbol, pymatgen_site.coords, siteid)

    def to_pymatgen_site(self):
        """
        :return: a pmg_site
        """
        return Site(self.element.symbol, self.coords, properties={"siteid": self.siteid})
