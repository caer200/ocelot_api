import sys
from copy import deepcopy

import numpy as np
from scipy.spatial.distance import pdist, squareform

from api.routines.loop import Loopsearcher
from api.schema.Backbone import Backbone
from api.schema.Ring import Ring
from api.schema.Sidechain import Sidechain
from api.schema.msite import MSite
from api.schema.msitelist import Sitelist


class Mol(Sitelist):

    def __init__(self, sites, name=None, charge=0, multiplicity=1):
        super().__init__(deepcopy(sites))

        self.name = name
        self.comment = {}
        self.charge = charge
        self.multiplicity = multiplicity

        # these are originally matrices
        self.coordmat = self.get_coordmat()
        self.distmat = self.get_distmat()
        self.bondmat = self.get_bondmat()
        self.nbrmap = self.get_nbrmap()

        # set mol site attributes 
        for i in range(len(self.sites)):
            self.sites[i].id = i
            self.sites[i].nbs = [self.sites[j] for j in self.nbrmap[i]]  # do not change connection table in a mol
            self.sites[i].num_nbs = len(self.sites[i].nbs)
            if isinstance(self.sites[i].element.valence, int):
                self.sites[i].insaturation = self.sites[i].element.valence - self.sites[i].num_nbs
                if self.sites[i].insaturation < 0:
                    sys.exit('when init mol site {} has {} nbs but only {} valence electrons'.format(
                        self.sites[i].id, len(self.sites[i].nbs), self.sites[i].element.valence))
            else:
                self.sites[i].insaturation = None

        rs = Loopsearcher(self.nbrmap)
        self.rings = []
        for ring_size in range(3, 9):
            idxlsts = rs.alex_method(ring_size)
            self.rings += [Ring.from_idxlst(idxlst, self) for idxlst in idxlsts]
        self.backbone = Backbone.from_mol(self)
        self.sidechains = self.get_sidechains()

    def get_distmat(self):
        return squareform(pdist(self.coordmat))

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

    # def get_shortest_path(s1, s2, path=None):
        # backtracking as in
        # https://www.python.org/doc/essays/graphs/
        # turns out to be problematic if the starting node has two directions to loop over
        # if path == None:
        #     path = [s1]
        # else:
        #     path += [s1]
        # if s1 == s2:
        #     return path
        # shortest = None
        # for node in s1.nbs:
        #     if node not in path:
        #         newpath = self.get_shortest_path(node, s2, path)
        #         if newpath:
        #             if not shortest or len(newpath) < len(shortest):
        #                 shortest = newpath

    @staticmethod
    def get_shortest_path(s1id, s2id, nbrmap):

        # dijsktra algo as in http://benalexkeen.com/implementing-djikstras-shortest-path-algorithm-with-python/
        # looks like percolation to me
        shortest_paths = {s1id: (None, 0)}
        current_node = s1id
        visited = set()

        while current_node != s2id:
            visited.add(current_node)
            destinations = nbrmap[current_node]
            weight_to_current_node = shortest_paths[current_node][1]

            for next_node in destinations:
                weight = 1 + weight_to_current_node
                if next_node not in shortest_paths:
                    shortest_paths[next_node] = (current_node, weight)
                else:
                    current_shortest_weight = shortest_paths[next_node][1]
                    if current_shortest_weight > weight:
                        shortest_paths[next_node] = (current_node, weight)

            next_destinations = {node: shortest_paths[node] for node in shortest_paths if node not in visited}
            if not next_destinations:
                return "Route Not Possible"
                # return None
            # next node is the destination with the lowest weight
            current_node = min(next_destinations, key=lambda k: next_destinations[k][1])

        # Work back through destinations in shortest path
        path = []
        while current_node is not None:
            path.append(current_node)
            next_node = shortest_paths[current_node][0]
            current_node = next_node
        # Reverse path
        path = path[::-1]
        return path

    def get_bondmat(self, co=1.3):
        """
        Bij = whether there is a bond between si and sj
        i is NOT bonded with itself
        :param co: coefficient for cutoff
        :return:
        """
        numsites = len(self.sites)
        mat = np.ones((numsites, numsites), dtype=bool)
        for i in range(numsites):
            mat[i][i] = False
            for j in range(i + 1, numsites):
                if self.distmat[i][j] > (self.sites[i].element.covrad + self.sites[j].element.covrad) * co:
                    mat[i][j] = False
                    mat[j][i] = False
        return mat

    def frlst(self):
        # fused rings, [[r1, r2, r3], [r4], [r5, r6], ...]
        indices = range(len(self.rings))
        block_list = []
        visited = []
        while len(visited) != len(self.rings):
            unvisited = [idx for idx in indices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in indices if idx not in block and idx not in visited]
                for i in outside:
                    if self.rings[block[pointer]].isfused_with(self.rings[i]):
                        block.append(i)
                visited.append(block[pointer])
                pointer += 1
            block_list.append([self.rings[j] for j in block])
        return block_list

    def largest_fr(self):
        """
        reture the largest fused ring system
        :return: a list of ring objs
        """
        frlst = self.frlst()
        if len(frlst) < 0:
            sys.exit('no ring identified while you want to grep the largest fused ring system!')
        frlst.sort(key=lambda x: len(x), reverse=True)
        return frlst[0]

    @classmethod
    def from_xyz(cls, xyz, name=None, charge=0, multiplicity=1):
        with open(xyz, 'r') as f:
            lines = f.readlines()[2:]
        sites = []
        for line in lines:
            items = line.strip().split()
            if len(items) == 4:
                s = MSite(items[0], items[1:])
                sites.append(s)
        return cls(sites, name, charge, multiplicity)

    @classmethod
    def from_gjf(cls, gjf, name=None):
        gjfname = gjf[:-4].split('/')[-1]
        with open(gjf, 'r') as f:
            lines = f.readlines()
        idx = 0
        hashline = ''
        for i in range(len(lines)):
            if lines[i].strip().startswith('#'):
                hashline = lines[i].strip()
                idx = i
                break
        if hashline == '':
            sys.exit('no hashline found in this gjf file!')

        charge, multiplicity = [int(s) for s in lines[idx + 4].strip().split()]
        sites = []
        for i in range(idx + 5, len(lines)):
            items = lines[i].strip().split()
            if len(items) == 4:
                sites.append(MSite(items[0], items[1:]))
        if name is None:
            mol = cls(sites, gjfname, charge, multiplicity)
        else:
            mol = cls(sites, name, charge, multiplicity)
        mol.comment['hashline'] = hashline
        return mol

    def get_sidechains(self):
        bs = self.backbone.sites
        bsids = []
        nbsids = []
        for s in self.sites:
            if s in bs:
                bsids.append(s.id)
            else:
                nbsids.append(s.id)
        side_chains = []
        side_chains_idx = []
        for sid in bsids:
            for nid in self.nbrmap[sid]:
                if nid not in bsids:
                    side_chains_idx.append([sid, nid])
        scid = 0

        for side_chain_ids in side_chains_idx:
            bj, sj = side_chain_ids
            for sid in nbsids:
                path = self.get_shortest_path(sid, bj, self.nbrmap)
                if len(set(path).intersection(set(bsids))) == 1:
                    if sid not in (bj, sj):
                        side_chain_ids.append(sid)
            side_chains.append(Sidechain([self.sites[ii] for ii in side_chain_ids], self, id=scid))
            scid += 1
        return side_chains
