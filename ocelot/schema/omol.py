import warnings
from rdkit import Chem
# import sys
# from copy import deepcopy
import sys
from pymatgen.core.structure import Molecule
from ocelot.routines.loop import Loopsearcher
from ocelot.schema.msitelist import MSitelist
from ocelot.schema.msite import MSite
from ocelot.schema.ring import Ring
from ocelot.schema.backbone import Backbone
from ocelot.schema.sidechain import Sidechain


class OMol(MSitelist):
    def __init__(self, msites, keepsiteid=False):
        """
        siteid is set only once when init this obj

        :param msites: a list of msites
        :param bool keepsiteid: rarely used, default Faluse

        :var distmat: distance matrix, not a property!
        :var bondmat: bond bool matrix, not a property!
        :var nrnbrmap: index map as a dict, not a property!
        :var rings: a list of Ring objs
        :var fused_rings_list: [[r1, r2, r3], [r5, r6], [r4]...], rn is a Ring obj
        :var is_solvent: True if there is only one ring in this omol
        :var largest_fused_ring: None if solvent, otherwise a list of Ring objs
        :var backbone: Backbone obj, None if solvent
        :var sidechains: a list of sc objs, None if solvent
        """
        # mss = [MSite(ms.element.name, deepcopy(ms.coords)) for ms in msites]
        # super().__init__(mss)

        super().__init__(msites)
        # these are original matrices
        self.distmat = self.get_distmat(self.coordmat)
        self.bondmat = self.get_bondmat(self.msites, self.distmat)
        self.nbrmap = self.get_nbrmap(self.bondmat)

        # set mol site attributes
        for i in range(len(self)):
            if keepsiteid:
                if self.msites[i].siteid == -1:
                    sys.exit('E: you want to keep site id when init an omol but there is at least one site with site_id as -1')
            else:
                self.msites[i].siteid = i
            self.msites[i].nbs = [self.msites[j] for j in self.nbrmap[i]]  # do not change connection table in a mol
            self.msites[i].nbs_idx = self.nbrmap[i]
            self.msites[i].num_nbs = len(self.msites[i].nbs)
            if isinstance(self.msites[i].element.valence, int):
                self.msites[i].insaturation = self.msites[i].element.valence - self.msites[i].num_nbs
                if self.msites[i].insaturation < 0:
                    warnings.warn('W: when init mol site {} has {} nbs but only {} valence electrons'.format(
                        self.msites[i].siteid, len(self.msites[i].nbs), self.msites[i].element.valence))
            else:
                self.msites[i].insaturation = None

        rs = Loopsearcher(self.nbrmap)
        idxlsts = rs.sssr_alex(3, 12)
        self.rings = [Ring.from_idxlst(idxlst, self.msites) for idxlst in idxlsts]

        # # somehow this doesnt work, GetRingInfo().AtomRings() returns an empty tuple
        # self.rings = []
        # # you need to make sure the order in msites is identical to siteid
        # for idxlst in self.rdkit_mol.GetRingInfo().AtomRings():
        #     self.rings.append(Ring.from_idxlst(idxlst, self.msites))

        self.fused_rings_list = self.get_fused_rings_list()  # [[r1, r2, r3], [r5, r6], [r4]...]

        if len(self.rings) <= 1:  # we think this omol is a solvent if it only has one or zero ring
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
        """
        keys are

        msites, can, volume, is_solvent

        if not solvent, additional:

        backbone, sidechains, rings, largest_fused_ring
        """
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

    @classmethod
    def from_dict(cls, d):
        msites = [MSite.from_dict(msd) for msd in  d['msites']]
        return cls(msites, keepsiteid=True)

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
        """
        dijsktra algo as in http://benalexkeen.com/implementing-djikstras-shortest-path-algorithm-with-python/

        looks like percolation to me

        :param s1id: siteid of site1
        :param s2id: siteid of site2
        :param nbrmap: nb map of this omol
        :return: "Route Not Possible" if disconnected, a list of idx if connected
        """

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
        """
        sorted (number of rings!) fused rings, [[r1, r2, r3], [r5, r6], [r4]...]

        TODO what if [[r1, r2, r3], [r4, r5, r6]] ?

        :return: a list of list of rings
        """
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
        """
        :return: a list of sc objs, notice sc.scid will be set here
        """
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
