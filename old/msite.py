import numpy as np
from constants import atom_electrons
# from ops import edist


class MSite:
    def __init__(self, symbol, coords, mid=-1, properties=None):
        """
        :param mid: id in a molecule
        :param symbol:
        :param coords:
        :param properties:
        """
        self.symbol = symbol
        if properties is None:
            self.properties = {'number': atom_electrons[self.symbol]}
        self.mid = mid
        self.coords = np.empty(3)
        for i in range(3):
            self.coords[i] = float(coords[i])
        self.x, self.y, self.z = self.coords

    def distance(self, other):

        return np.linalg.norm(self.coords - other.coords)

    def __eq__(self, other):
        if self.symbol == other.symbol and np.allclose(self.coords, other.coords):
            if set(self.properties.keys()) == set(other.properties.keys()):
                if all([self.properties[k] == other.properties[k] for k in self.properties.keys()]):
                    return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.symbol)

    def __contains__(self, item):
        for i in item:
            if self.__eq__(i):
                return True
        return False

    def __repr__(self):
        return "MolSite: {} ({:.4f}, {:.4f}, {:.4f}) {} {}".format(
            self.symbol, *self.coords, self.properties, self.mid)

    def __str__(self):
        return "MolSite: {} ({:.4f}, {:.4f}, {:.4f}) {} {}".format(
            self.symbol, *self.coords, self.properties, self.mid)
