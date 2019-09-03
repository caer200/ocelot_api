from Routines import *
from Schema import *
import itertools
import sympy


class FPBC_optimizer:

    def __init__(self, pbc, mol):
        self.mol = mol
        self.pbc = np.zeros((len(pbc), 3))
        for i in range(len(pbc)):
            for j in range(3):
                self.pbc = pbc[i][j]

    def gen_transv(self, coelim=1, linear_dependent=False):
        vs = []
        if linear_dependent:
            coes = itertools.product(range(coelim+1), range(coelim+1), range(coelim+1))
            coes = np.array(list(coes))
            # https://stackoverflow.com/questions/28816627/ does not work for floats?
            _, inds = sympy.Matrix(coes).T.rref()
            coes = coes[inds]
            for coe in coes:
                vs.append(coe[0]*self.pbc[0]+coe[1]*self.pbc[1]+coe[2]*self.pbc[2])
        return vs

    def nbchecker(self, pbc, criterion='fixed'):
        if criterion == 'fixed':
            SiteList.translate_copy(sl, )

