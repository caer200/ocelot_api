import numpy as np
from ocelot.routines.geometry import norm
from ocelot.schema.conformation.atomlist import AtomList, AtomSite


class Bond(AtomList):

    def __init__(self, sitea, siteb):
        """
        try not to init it outside an omol

        :param sitea: first site, self.a
        :param siteb: second site, self.b
        """
        super().__init__([sitea, siteb])
        self.checkstatus('all assigned', 'unique ids')
        self.a = sitea
        self.b = siteb

    # def order(self, rdmol):
    #     bt = rdmol.GetBondBetweenAtoms(self.a.siteid, self.b.siteid).GetBondType()
    #     if bt == Chem.BondType.SINGLE:
    #         return 1
    #     elif bt == Chem.BondType.DOUBLE:
    #         return 2
    #     elif bt == Chem.BondType.TRIPLE:
    #         return 3
    #     elif bt == Chem.BondType.AROMATIC:
    #         return 1.5
    #     else:
    #         return None

    @property
    def center(self):
        """
        :return: cart coords of geo center
        """
        return (self.a.coords + self.b.coords) * 0.5

    @property
    def length(self):
        return norm(self.a.coords - self.b.coords)

    @property
    def elements(self):
        """
        :return: 2x1 string tuple
        """
        return tuple(sorted([self.a.element.symbol, self.b.element.symbol]))

    def __eq__(self, other):
        # using center should be enough for msites in a mol
        if self.elements == other.elements:
            match = 0
            for ss in self.atomsites:
                for so in other.atomsites:
                    if np.allclose(ss.coords, so.coords):
                        match += 1
                        break
            if match == 2:
                return 1
        return 0

    def as_dict(self):
        """
        keys are

        site_a, site_b, order, length, center
        :return: a dict
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            'site_a': self.a.as_dict(),
            'site_b': self.b.as_dict(),
            'length': self.length,
            'center': self.center,
        }
        return d

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        site_a, site_b
        """
        sa = AtomSite.from_dict(d['site_a'])
        sb = AtomSite.from_dict(d['site_b'])
        return cls(sa, sb)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __repr__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)

    def __str__(self):
        return "Bond: {} \n Center: ({:.4f}, {:.4f}, {:.4f}) \n Length: {}".format(
            self.elements, *self.center, self.length)
