from ocelot.schema import Sitelist
from ocelot.routines import norm
import numpy as np
import sys


class Bond(Sitelist):

    def __init__(self, sites):
        super().__init__(sites)
        if len(self.sites) != 2:
            sys.exit('you are initializing a bond with <2 or >2 msites')

    @property
    def a(self):
        return self.sites[0]

    @property
    def b(self):
        return self.sites[1]

    @property
    def order(self):
        if self.a.element.valence == 1 or self.b.element.valence == 1:
            return 1
        return min([s.insaturation for s in self.sites])

    @property
    def center(self):
        return (self.a.coords + self.b.coords)*0.5

    @property
    def length(self):
        return norm(self.a.coords - self.b.coords)

    @property
    def elements(self):
        return {self.a.element, self.b.element}

    def __eq__(self, other):
        # using center should be enough for msites in a mol
        if self.elements == other.elements:
            match = 0
            for ss in self.sites:
                for so in other.sites:
                    if np.allclose(ss.coords, so.coords):
                        match += 1
                        break
            if match == 2:
                return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)

    def __str__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)
