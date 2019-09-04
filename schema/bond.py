import numpy as np
import warnings
from routines.geometry import norm
from schema.msitelist import MSitelist, MSite


class Bond(MSitelist):

    def __init__(self, sitea, siteb):
        super().__init__([sitea, siteb])
        self.omol_init = True
        for ms in self.msites:
            if ms.siteid == -1:
                warnings.warn('W: you are init a bond with sites not in an omol obj')
                self.omol_init = False
                break
        self.a = sitea
        self.b = siteb

    @property
    def order(self):  # self.omol_init = True
        if not self.omol_init:
            return None
        if self.a.element.valence == 1 or self.b.element.valence == 1:
            return 1
        return min([s.insaturation for s in self.msites])

    @property
    def center(self):
        return (self.a.coords + self.b.coords) * 0.5

    @property
    def length(self):
        return norm(self.a.coords - self.b.coords)

    @property
    def elements(self):
        return tuple(sorted([self.a.element.name, self.b.element.name]))

    def __eq__(self, other):
        # using center should be enough for msites in a mol
        if self.elements == other.elements:
            match = 0
            for ss in self.msites:
                for so in other.msites:
                    if np.allclose(ss.coords, so.coords):
                        match += 1
                        break
            if match == 2:
                return 1
        return 0

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            'site_a': self.a.as_dict(),
            'site_b': self.b.as_dict(),
            'order': self.order,
            'length': self.length,
            'center': self.center,
        }
        return d

    @classmethod
    def from_dict(cls, d):
        sa = MSite.from_dict(d['site_a'])
        sb = MSite.from_dict(d['site_b'])
        return cls(sa, sb)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)

    def __str__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)
