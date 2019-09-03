from copy import deepcopy
from scipy.spatial.distance import pdist, squareform, cdist
from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.core.structure import Structure
from api.schema.element import Element
import numpy as np
import itertools

class PBCmol:

    def __init__(self, pstructure, name="", comment=None):
        self.pstructure = pstructure
        if comment is None:
            self.comment = dict()
        else:
            self.comment = comment
        self.name = name

        self.coordmat = self.pstructure.coords

        self.shortest_vs, d2 = pbc_shortest_vectors(self.pstructure.lattice, self.coordmat, self.coordmat, return_d2=True)
        self.distmat = np.sqrt(d2)

        self.bondmat = self.get_bondmat()
        self.nbrmap = self.get_nbrmap()

    def get_bondmat(self, co=1.3):
        """
        Bij = whether there is a bond between si and sj
        i is NOT bonded with itself
        :param co: coefficient for cutoff
        :return:
        """
        numsites = len(self.pstructure.sites)
        mat = np.zeros((numsites, numsites), dtype=bool)
        for i in range(numsites):
            for j in range(i + 1, numsites):
                if self.distmat[i][j] < (Element(self.pstructure.sites[i].species_string).covrad + Element(
                        self.pstructure.sites[j].species_string).covrad) * co:
                    mat[i][j] = True
                    mat[j][i] = True
        return mat

    def get_nbrmap(self):
        # TODO profile
        bmat = self.bondmat
        ma = {}
        numsites = len(bmat)
        for i in range(numsites):
            nbs = []
            for j in range(numsites):
                if bmat[i][j] and j != i:
                    nbs.append(j)
            ma[i] = nbs
        return ma

    # def to_mols(self):
    #     """
    #     unwrap the poscar and return a list of mols sorted by # of msites
    #     :return:
    #     """
    #
    @classmethod
    def unwrap(cls, pstructure):
        """
        return a new unwrap structure
        :return:
        """
        psites = deepcopy(pstructure.psites)
        nsites = len(pstructure.psites)
        pindices = range(nsites)
        visited = []
        block_list = []
        # unwrap = []
        # unwrap_block_list = []
        while len(visited) != nsites:
            # initialization
            unvisited = [idx for idx in pindices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            # unwrap.append(ini_idx)
            # unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx]._coords)]
            # unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords)]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    if i in pstructure.nbrmap[block[pointer]]:
                        block.append(i)
                        psites[i].coords = psites[i].coords + pstructure.shortest_transvs[]

                    # if distance < cutoff:
                    #     block.append(i)
                    #     psites[i] = PeriodicSite(psites[i].species_string, psites[i]._fcoords + fctrans,
                    #                              pstructure.lattice)
                    #     unwrap_block.append(Site(psites[i].species_string, psites[i]._coords))
                    #     # psites[i] = PeriodicSite(psites[i].species_string, psites[i].fcoords + fctrans, pstructure.lattice)
                    #     # unwrap_block.append(Site(psites[i].species_string, psites[i].coords))
                    #     unwrap.append(psites[i])
                visited.append(block[pointer])
                pointer += 1
            # unwrap_block_list.append(unwrap_block)
            block_list.append(block)
        # mols = [IMolecule.from_sites(i) for i in unwrap_block_list]
        # # for i in range(len(mols)):
        # #     mols[i].to(fmt='xyz', filename='mol_' + str(i) + '.xyz')
        # unwrap = sorted(unwrap, key=lambda x: x.species_string)
        # unwrap_str = IStructure.from_sites(unwrap)
        # unwrap_str.to(fmt='poscar', filename='unwrap.poscar')


