from ocelot.schema.element import Element
import numpy as np
from pymatgen.core.sites import Site
import sys
import warnings


class MSite:
    def __init__(self, element_name, coords, siteid=-1):
        """
        :param str element_name: case sensitive
        :param np.ndarray coords: 3*1, cart always
        :param int siteid: default -1, this should be set to non-negative only once when a OMol obj is initiated
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
        # TODO profile
        return np.linalg.norm(self.coords - other.coords)

    def __eq__(self, other):
        """
        two msites are equal if they have the same element and are close in cart system

        this is useful for set or unique list

        :param MSite other:
        """
        if self.element == other.element and np.allclose(self.coords, other.coords) and self.siteid == other.siteid:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        """
        rarely used, just in case
        """
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
        keys are

        element, xyz, siteid
        :return: json-serializable representation for MSite.
        :rtype: dict
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
        create MSite from dict representation, return None if failed

        keys are

        element, xyz, siteid
        :param dict d: representation of a MSite
        :return: MSite, or None if keys not found
        """
        try:
            return cls(d["element"], d["xyz"], d["siteid"])
        except KeyError:
            warnings.warn('W: cannot init MSite from dict as no required keys')
            return None

    @classmethod
    def from_pymatgen_site(cls, pymatgen_site):
        """
        create MSite from a pmg Site obj, the coords for that pmg site should be in cart and siteid is set to -1
        :param pmg_site pymatgen_site: the pmg site that will be converted
        :return: MSite
        """
        return cls(pymatgen_site.species_string, pymatgen_site.coords)

    def to_pymatgen_site(self):
        """
        notice the siteid info is lost in this conversion
        :return: a pmg_site
        """
        # return Site(self.element.name, self.coords, properties={"siteid": self.siteid})
        return Site(str(self.element.name), self.coords)
