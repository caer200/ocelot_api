from copy import deepcopy
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.util.coord import pbc_shortest_vectors
from api.schema.element import Element
import numpy as np


class Unwrapper:

    def __init__(self, pstructure, co=1.3):
        self.pstructure = deepcopy(pstructure)
        self.co = co
        self.sites = self.pstructure.sites
        self.nsites = len(self.sites)
        self.lattice = self.pstructure.lattice
        self.coordmat = self.pstructure.coords
        self.shortest_vs, d2 = pbc_shortest_vectors(self.lattice, self.coordmat, self.coordmat, return_d2=True)
        self.distmat = np.sqrt(d2)

        self.bondmat = np.zeros((self.nsites, self.nsites), dtype=bool)
        for i in range(self.nsites):
            for j in range(i+1, self.nsites):
                ielement = self.sites[i].species_string
                jelement = self.sites[j].species_string
                cutoff = (Element(ielement).covrad + Element(jelement).covrad) * self.co
                if self.distmat[i][j] < cutoff:
                    self.bondmat[i][j] = True
                    self.bondmat[j][i] = True
        self.nbrmap = dict()
        for i in range(self.nsites):
            nbs = []
            for j in range(self.nsites):
                if self.bondmat[i][j] and j != i:
                    nbs.append(j)
            self.nbrmap[i] = nbs

    @staticmethod
    def get_dist_and_trans(fc1, fc2, lat):
        """
        get the shortest distance and corresponding translation vector between two frac coords
        :param fc1:
        :param fc2:
        :param lat:
        :return:
        """
        v, d2 = pbc_shortest_vectors(lat, fc1, fc2, return_d2=True)
        fc = lat.get_fractional_coords(v[0][0]) + fc1 - fc2
        return np.sqrt(d2[0, 0]), fc

    def unwrap(self):
        pindices = range(self.nsites)
        visited = []
        block_list = []
        # group = []
        # group_block_list = []
        while len(visited) != self.nsites:
            # initialization
            unvisited = [idx for idx in pindices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            # group.append(ini_idx)
            # group_block = [Site(self.msites[ini_idx].species_string, self.msites[ini_idx]._coords)]
            # unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords)]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    if i in self.nbrmap[block[pointer]]:
                        block.append(i)
                        self.sites[]._fcoords = self.sites[i]._fcoords + self.shortest_vs[block[pointer]][i]

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

























        pindices = range(len(self.sites))
        visited = []
        block_list = []
        unwrap = []
        unwrap_block_list = []
        while len(visited) != len(self.sites):
            # initialization
            unvisited = [idx for idx in pindices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            unwrap.append(self.sites[ini_idx])
            unwrap_block = [Site(self.sites[ini_idx].species_string, self.sites[ini_idx]._coords)]
            # unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords)]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    distance, fctrans = get_dist_and_trans(psites[block[pointer]]._fcoords, psites[i]._fcoords,
                                                           pstructure.lattice)
                    # distance, fctrans = get_dist_and_trans(psites[block[pointer]].fcoords, psites[i].fcoords,
                    #                                        pstructure.lattice)
                    cutoff = get_cutoff(psites[block[pointer]].species_string, psites[i].species_string)
                    if distance < cutoff:
                        block.append(i)
                        psites[i] = PeriodicSite(psites[i].species_string, psites[i]._fcoords + fctrans,
                                                 pstructure.lattice)
                        unwrap_block.append(Site(psites[i].species_string, psites[i]._coords))
                        # psites[i] = PeriodicSite(psites[i].species_string, psites[i].fcoords + fctrans, pstructure.lattice)
                        # unwrap_block.append(Site(psites[i].species_string, psites[i].coords))
                        unwrap.append(psites[i])
                visited.append(block[pointer])
                pointer += 1
            unwrap_block_list.append(unwrap_block)
            block_list.append(block)
        mols = [IMolecule.from_sites(i) for i in unwrap_block_list]
        # for i in range(len(mols)):
        #     mols[i].to(fmt='xyz', filename='mol_' + str(i) + '.xyz')
        unwrap = sorted(unwrap, key=lambda x: x.species_string)
        unwrap_str = IStructure.from_sites(unwrap)
        unwrap_str.to(fmt='poscar', filename='unwrap.poscar')
        return mols, unwrap_str
