import sys

from ocelot.routines import angle_btw
from ocelot.schema import Sitelist


class Sidechain(Sitelist):

    def __init__(self, sites, mol, id=-1):
        """
        the first site in msites is bone_joint, the obj has msites started from side_joint
        :param sites:
        :param mol:
        :param id:
        """
        if len(sites) < 2:
            sys.exit('Sidechain obj must be init with at least 2 msites!')
        super().__init__(sites[1:])
        self.id = id
        self.mol = mol
        self.bone_joint = sites[0]
        self.side_joint = self.sites[0]

        for i in range(len(self.sites)):
            self.sites[i].path_to_bj = self.mol.get_shortest_path(self.sites[i].id, self.bone_joint.id, self.mol.nbrmap)
            self.sites[i].rank = len(self.sites[i].path_to_bj)  # bone_joint has rank 1, side_joint has rank 2
        self.bone_joint.rank = 1
        for s in self.sites:
            lower_rank_nbs = []
            eq_rank_nbs = []
            higher_rank_nbs = []
            for nb in s.nbs:
                if nb.rank < s.rank:
                    lower_rank_nbs.append(nb)
                elif nb.rank == s.rank:
                    eq_rank_nbs.append(nb)
                else:
                    higher_rank_nbs.append(nb)
            s.lower_rank_nbs = lower_rank_nbs
            s.eq_rank_nbs = eq_rank_nbs
            s.higher_rank_nbs = higher_rank_nbs
            if len(s.lower_rank_nbs) > 1:
                sys.exit('inner loop in this side chain, id=' + str(self.id))
            # TODO add more hybrid identification rules
            if s.element.name in ['C', 'Si', 'Ge', 'Sn']:
                if len(s.nbs) == 4:
                    s.hybrid = 'sp3'
                elif len(s.nbs) == 3:
                    s.hybrid = 'sp2'
                elif len(s.nbs) == 2:
                    s.hybrid = 'sp'
                else:
                    s.hybrid = None
            else:
                s.hybrid = None

        self.rankmap = [None, None]
        for rank in range(2, self.maxrank):
            self.rankmap.append(sorted([s for s in self.sites if s.rank == rank], key=lambda x: len(x.higher_rank_nbs)))

    @property
    def ishydrogen(self):
        return len(self.sites) == 1 and self.sites[0].symbol == 'H'

    @property
    def angle_vp(self):
        v1 = self.bone_joint.coords - self.mol.backbone.geoc
        v2 = self.mol.backbone.vp
        avp = angle_btw(v1, v2, output='degree')
        # if avp > 90.0:
        #     self.avp = 180 - self.avp
        return avp

    @property
    def maxrank(self):
        return max([s.rank for s in self.sites])



    # def v_tree(self):
    #     tree = {}
    #     for rank in range(2, self.maxrank):
    #         for s in self.rankmap[rank]:
    #             if rank == 2:
    #                 vref = s.coords - self.bone_joint.coords
    #                 v = [np.dot(vref, self.mol.backbone.vp), np.dot(vref, self.mol.backbone.vq), np.dot(vref, self.mol.backbone.vo)]
    #                 stem = 0  # 'fixed'
    #             else:
    #     for rank in self.msites:
    #         if s.rank == 2:
    #         else:
    #             nbroot = s.lower_rank_nbs[0]
    #             vref = s.coords - nbroot.coords
    #             v = [np.dot(vref, self.mol.backbone.vp), np.dot(vref, self.mol.backbone.vq), np.dot(vref, self.mol.backbone.vo)]
    #             if nbroot.insaturation > 0:
    #                 stem = 0  # 'fixed'
    #             else:
    #                 stem = 1  # 'flex'
    #         tree[s.id] = [v, stem]


    # TODO add similarity and morphology descriptors
