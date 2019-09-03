import warnings
import sys
from copy import deepcopy
from pymatgen.core.structure import Molecule
from api.routines.loop import Loopsearcher
from api.schema.msitelist import MSitelist
from api.schema.msite import MSite
from api.schema.ring import Ring
from api.schema.backbone import Backbone
from api.schema.sidechain import Sidechain


class OMol(MSitelist):
    """
    siteid is set only once when init this obj
    """

    def __init__(self, msites):
        # mss = [MSite(ms.element.name, deepcopy(ms.coords)) for ms in msites]
        # super().__init__(mss)

        super().__init__(msites)
        # these are original matrices
        self.distmat = self.get_distmat(self.coordmat)
        self.bondmat = self.get_bondmat(self.msites, self.distmat)
        self.nbrmap = self.get_nbrmap(self.bondmat)

        # set mol site attributes
        for i in range(len(self)):
            self.msites[i].siteid = i
            self.msites[i].nbs = [self.msites[j] for j in self.nbrmap[i]]  # do not change connection table in a mol
            self.msites[i].nbs_idx = self.nbrmap[i]
            self.msites[i].num_nbs = len(self.msites[i].nbs)
            if isinstance(self.msites[i].element.valence, int):
                self.msites[i].insaturation = self.msites[i].element.valence - self.msites[i].num_nbs
                if self.msites[i].insaturation < 0:
                    warnings.warn('when init mol site {} has {} nbs but only {} valence electrons'.format(
                        self.msites[i].siteid, len(self.msites[i].nbs), self.msites[i].element.valence))
            else:
                self.msites[i].insaturation = None

        rs = Loopsearcher(self.nbrmap)
        self.rings = []
        for ring_size in range(3, 9):
            idxlsts = rs.alex_method(ring_size)
            self.rings += [Ring.from_idxlst(idxlst, self.msites) for idxlst in idxlsts]

        self.fused_rings_list = self.get_fused_rings_list()  # [[r1, r2, r3], [r5, r6], [r4]...]

        if len(self.rings) <= 1:  # we think this is a solvent
            self.largest_fused_ring = None
            self.is_solvent = True
            self.backbone = None
            self.sidechains = None
        else:
            self.largest_fused_ring = self.fused_rings_list[0]
            self.is_solvent = False
            self.backbone = Backbone.from_omol(self)
            self.sidechains = self.get_sidechains()

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "msites": [s.as_dict() for s in self.msites],
            "can": self.canonical_smiles,
            "volume": self.volume,
            "is_solvent": self.is_solvent,
        }
        if not self.is_solvent:
            d['backbone'] = self.backbone.as_dict()
            d['sidechains'] = [sidechain.as_dict() for sidechain in self.sidechains]
            d['rings'] = [r.as_dict() for r in self.rings]
            d['largest_fused_ring'] = [r.as_dict() for r in self.largest_fused_ring]
        return d

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

    def get_fused_rings_list(self):
        # sorted fused rings, [[r1, r2, r3], [r5, r6], [r4]...]
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
        block_list.sort(key=lambda x: len(x), reverse=True)
        return block_list

    def get_sidechains(self):
        bs = self.backbone.msites
        bsids = []
        nbsids = []
        for s in self.msites:
            if s in bs:
                bsids.append(s.siteid)
            else:
                nbsids.append(s.siteid)
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
            side_chain_obj = Sidechain.from_omol([self.msites[ii] for ii in side_chain_ids], scid=scid, omol=self)
            side_chains.append(side_chain_obj)
            scid += 1
        return side_chains

    @classmethod
    def from_pymatgen_mol(cls, mol):
        msites = [MSite.from_pymatgen_site(s) for s in mol.sites]
        return cls(msites)

    def to_pymatgen_mol(self):
        sites = [MSite.to_pymatgen_site(s) for s in self.msites]
        return Molecule.from_sites(sites)
