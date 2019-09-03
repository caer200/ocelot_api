import warnings
from api.routines.geometry import angle_btw
from api.schema.msitelist import MSitelist


class Sidechain(MSitelist):

    def __init__(self, sc_msites, bone_joint, scid, rankmap, angle_vp, hasring=False):
        if len(sc_msites) < 1:
            warnings.warn('Sidechain obj must be init with at least 1 sc msite!')
        super().__init__(sc_msites)
        self.scid = scid
        self.bone_joint = bone_joint
        self.side_joint = self.msites[0]
        self.rankmap = rankmap
        self.angle_vp = angle_vp
        self.hasring = hasring  # inner loop

        ranksorted_sites = sorted(self.msites, key=lambda x: x.rank)
        self.branch_msite = None
        self.branch_msite_rank = None
        self.umbrella_sites = []
        self.umbrella = None

        for ss in ranksorted_sites:
            if len(ss.higher_rank_nbs) > 1 and len(ss.lower_rank_nbs) == 1:
                self.branch_msite = ss
                self.branch_msite_rank = ss.rank
                self.umbrella_sites = [s for s in self.msites if s.rank >= self.branch_msite_rank]
                self.umbrella = MSitelist(self.umbrella_sites)
                break

    def as_dict(self):
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "can": self.canonical_smiles,
            "volume": self.volume,
            "msites": [s.as_dict() for s in self.msites],
            "max_rank": self.maxrank,
            "branch_rank": self.branch_msite_rank,
            "is_hydrogen": self.ishydrogen,
            "has_ring": self.hasring,
            "scid": self.scid,
            "angle_vp": self.angle_vp
        }
        return d

    @property
    def ishydrogen(self):
        return len(self.msites) == 1 and self.msites[0].element.name == 'H'

    @classmethod
    def from_omol(cls, msties, scid, omol):
        """

        :param msties: including bone-joint site
        :param scid:
        :param omol:
        :return:
        """
        hasring = False  # inner loop
        sc_msites = msties[1:]
        bone_joint = msties[0]
        for i in range(len(sc_msites)):
            sc_msites[i].path_to_bj = omol.get_shortest_path(sc_msites[i].siteid, bone_joint.siteid, omol.nbrmap)
            sc_msites[i].rank = len(sc_msites[i].path_to_bj)  # bone_joint has rank 1, side_joint has rank 2
            bone_joint.rank = 1
        for s in sc_msites:
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
                hasring = True
                warnings.warn('inner loop in this side chain, scid={}'.format(scid))
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
        maxrank = max([s.rank for s in sc_msites])
        rankmap = [None, None]  # rankmap[2] is [sidejoint], rankmap[3] is [msites with rank==3]
        for rank in range(2, maxrank):
            rankmap.append(sorted([s for s in sc_msites if s.rank == rank], key=lambda x: len(x.higher_rank_nbs)))

        v1 = bone_joint.coords - omol.backbone.geoc
        v2 = omol.backbone.vp_fit
        avp = angle_btw(v1, v2, output='degree')
        if avp > 90.0:
            avp = 180 - avp
        return cls(sc_msites, bone_joint, scid, rankmap, angle_vp=avp, hasring=hasring)

    @property
    def maxrank(self):
        return max([s.rank for s in self.msites])
