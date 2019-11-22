from ocelot.schema.omol import OMol
from ocelot.schema.element import nvalence_electron
from ocelot.routines.loop import Loopsearcher
# from copy import deepcopy

"""
identify whether a molecule is an OMOL based on connection map (graph)

OMOL consists of a backbone and a set of side chains
"""

def getnval(atom):
    symbol = atom.GetSymbol()
    valence_elect = nvalence_electron(symbol)
    return valence_elect

class Atomlist:
    def __init__(self, atoms):
        self.atoms = atoms
        self.idx = [atom.GetIdx() for atom in self.atoms]
        self.nvelect = 0
        for a in atoms:
            self.nvelect += getnval(a)

    def __repr__(self):
        outs = [self.__class__.__name__ + ': ']
        for s in self.atoms:
            outs.append(str(s.GetIdx()))
        return ' '.join(outs)

    def interscet(self, other):
        """
        get idx of shared sites

        :param Ring other:
        :return: a list of siteid
        """
        return list(set(self.idx) & set(other.idx))

    def isfused_with(self, other):
        """
        'fused' means the two lists (rings) share at least 2 sites

        :param Atomlist other:
        :return: bool
        """
        if len(self.interscet(other)) > 1:
            return 1
        return 0

    def __getitem__(self, ind):
        return self.atoms[ind]

    @property
    def unsat(self):
        """
        notice the unsat is calculated based on atom.GetNeighbors.
        usually atom is generated from mol.GetAtoms()
        this means the unsat is calculated based on mol, not this Atomlist

        :return:
        """
        n_unsat = 0
        for atom in self.atoms:
            if getnval(atom) - len(atom.GetNeighbors()) > 0:
                n_unsat += 1
        return n_unsat / len(self.atoms)

    @classmethod
    def from_idxlst(cls, idxlst, allatoms):
        atoms = []
        for i in idxlst:
            atoms.append(allatoms[i])
        return cls(atoms)


class SidechainRd(Atomlist):
    def __init__(self, atoms, bone_joint, scid, rankmap, hasring):
        super().__init__(atoms)
        # self.bone_joint = deepcopy(bone_joint)
        self.bone_joint = bone_joint
        self.side_joint = self.atoms[0]
        self.scid = scid
        self.rankmap = rankmap
        self.hasring = hasring

    @classmethod
    def from_omol(cls, atoms, scid, omolrd):
        """
        this method assigns several attributes to msites on the side_chain and bone_joint

        site.rank, again bone_joint has rank=1

        """
        hasring = False  # inner loop
        sc_atoms = atoms[1:]
        bone_joint = atoms[0]
        bone_joint.rank = 1
        for i in range(len(sc_atoms)):
            sc_atoms[i].path_to_bj = OMol.get_shortest_path(sc_atoms[i].GetIdx(), bone_joint.GetIdx(), omolrd.nbmap)
            sc_atoms[i].rank = len(sc_atoms[i].path_to_bj)  # bone_joint has rank 1, side_joint has rank 2
        maxrank = max([s.rank for s in sc_atoms])
        rankmap = [None, None]  # rankmap[2] is [sidejoint], rankmap[3] is [msites with rank==3]
        for rank in range(2, maxrank+1):
            rankmap.append([s for s in sc_atoms if s.rank == rank])

        return cls(sc_atoms, bone_joint, scid, rankmap, hasring=hasring)


class BackboneRd(Atomlist):
    def __init__(self, atoms):
        super().__init__(atoms)

    @classmethod
    def from_omol(cls, omolrd):
        """
        just largest_fused_ring in omolrd

        :param omolrd:
        :return:
        """
        fused_rings_list = omolrd.fused_rings_list  # [[r1, r2, r3], [r5, r6], [r4]...]
        atoms = []
        for ring in fused_rings_list[0]:
            atoms += ring.atoms
        return cls(atoms)


class RingRd(Atomlist):
    def __init__(self, atoms):
        super().__init__(atoms)


class OmolRd(Atomlist):

    def __init__(self, rdmol):
        self.rdmol = rdmol
        atoms = self.rdmol.GetAtoms()
        super().__init__(atoms)
        self.natoms = len(self.atoms)

        self.rings = []
        # you need to make sure the order in msites is identical to siteid
        # for idxlst in self.rdmol.GetRingInfo().AtomRings():
        #     self.rings.append(RingRd.from_idxlst(idxlst, self.atoms))
        rs = Loopsearcher(self.nbmap)
        idxlsts = rs.sssr_alex(3, 10)
        self.rings = [RingRd.from_idxlst(idxlst, self.atoms) for idxlst in idxlsts]
        self.fused_rings_list = self.get_fused_rings_list()  # [[r1, r2, r3], [r5, r6], [r4]...]
        if len(self.rings) > 0 and len(self.fused_rings_list) > 0:
            self.backbone, self.scs = self.omol_partition()
        else:
            self.backbone = None
            self.scs = []

    @property
    def nbmap(self):
        map = {}
        for i in range(self.natoms):
            atom = self.atoms[i]
            nbs = [nb.GetIdx() for nb in atom.GetNeighbors()]
            map[i] = nbs
        return map

    def get_fused_rings_list(self):
        """
        sorted (number of rings!) fused rings, [[r1, r2, r3], [r5, r6], [r4]...]

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

    def omol_partition(self):
        backbone = BackboneRd.from_omol(self)
        bsids = backbone.idx
        nbsids = [i for i in self.idx if i not in bsids]

        side_chains = []
        side_chains_idx = []
        for sid in bsids:
            for nid in self.nbmap[sid]:
                if nid not in bsids:
                    side_chains_idx.append([sid, nid])
        scid = 0
        for side_chain_ids in side_chains_idx:
            bj, sj = side_chain_ids
            for sid in nbsids:
                path = OMol.get_shortest_path(sid, bj, self.nbmap)
                if len(set(path).intersection(set(bsids))) == 1:
                    if sid not in (bj, sj):
                        side_chain_ids.append(sid)
            side_chain_obj = SidechainRd.from_omol([self.atoms[ii] for ii in side_chain_ids], scid=scid, omolrd=self)
            side_chains.append(side_chain_obj)
            scid += 1
        return backbone, side_chains

"""
it is possible to use sssr algo within the rdkit

m = Chem.MolFromSmiles('CCCCCC')
atoms = m.GetAtoms()
for ringidx in m.GetRingInfo().AtomRings():
    for i in ringidx:
        print(atoms[i].GetSymbol(), atoms[i].GetIdx(), ringidx)

the problem is somehow the rdkit molecule object has to be built from smiles or other files where connectivity 
is well-defined, otherwise GetRingInfo() is an empty object, e.g. if you get the molecule from xyz2mol by Jensen.
This behavior may be related to the fact that in the AC2BO function (xyz2mol.py) the atomic_valence is not fully defined
for every element, so bonds are not defined in this case.
"""
