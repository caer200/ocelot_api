import warnings
from ocelot.routines.geometry import angle_btw
from ocelot.schema.msitelist import MSitelist
from ocelot.schema.msite import MSite


class Sidechain(MSitelist):

    def __init__(self, sc_msites, bone_joint, scid, rankmap, angle_vp, hasring=False):
        """

        bone_joint(backbone) -- side_joint(sidechain)

        :param sc_msites: msites on the sidechain, start from side_joint
        :param bone_joint: the msite obj of bone_joint
        :param scid: id of this sidechain in the omol
        :param rankmap: rankmap is a list, always rankmap[0] = rankmap[1] = None, note bone_joint has rank as 1 but is not included here, rankmap[i] with i > 1 is a list of sites with rank i, rankmap[2] is side_joint msite
        :param angle_vp: the angle between v(bone_geoc --> bone_jopint) and v()
        :param hasring: boolean, inner loop within the side chain  TODO DAG?

        :var umbrella: a msitelist, considering a sidechain like -=-TIPS, the umbrella means the msites with rank that is equal or larger than Si
        :var branch_msite: Si in -=-TIPS

        """
        if len(sc_msites) < 1:
            warnings.warn('W: Sidechain obj must be init with at least 1 sc msite!')
        super().__init__(sc_msites)
        self.omol_init = True
        for ms in sc_msites:
            if ms.siteid == -1:
                warnings.warn('W: you are init a ring with sites not in an omol obj')
                self.omol_init = False
                break
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

    @property
    def site_ids(self):
        return [ms.siteid for ms in self.msites]

    def is_identical_with(self, other):
        return set(self.site_ids) == set(other.site_ids)

    def is_intersection_with(self, other):
        return bool(set(self.site_ids).intersection(set(other.site_ids)))

    def as_dict(self):
        """
        keys are

        can, volume, msites, max_rank, branch_rank, is_hydrogen, has_ring, scid, angle_vp, bone_joint, rankmap
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "can": self.canonical_smiles,
            "volume": self.volume,
            "msites": [s.as_dict() for s in self.msites],
            "max_rank": self.maxrank,
            "shape_descriptor": self.shape_descriptors,
            "branch_rank": self.branch_msite_rank,
            "is_hydrogen": self.ishydrogen,
            "has_ring": self.hasring,
            "scid": self.scid,
            "angle_vp": self.angle_vp,
            "bone_joint": self.bone_joint.as_dict(),
            "rankmap": [None, None] + [[s.as_dict() for s in sitelist] for sitelist in self.rankmap[2:]]
        }
        if self.umbrella is None:
            d['umbrella'] = None
        else:
            d['umbrella'] = self.umbrella.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        msites, bone_joint, scid, rankmap, angle_vp, hasring
        """
        sc_msites = [MSite.from_dict(sitedict) for sitedict in d['msites']]
        bone_joint = MSite.from_dict(d['bone_joint'])
        scid = d['scid']
        rankmap = d['rankmap']
        for i in range(2, len(rankmap)):
            rankmap[i] = [MSite.from_dict(sdict) for sdict in rankmap[i]]
        angle_vp = d['angle_vp']
        hasring = d['hasring']
        return cls(sc_msites, bone_joint, scid, rankmap, angle_vp, hasring)

    @property
    def ishydrogen(self):
        """
        :return: bool
        """
        return len(self.msites) == 1 and self.msites[0].element.name == 'H'

    @classmethod
    def from_omol(cls, msties, scid, omol):
        """
        this method assigns several attributes to msites on the side_chain and bone_joint

        site.rank, again bone_joint has rank=1

        site.lower_rank_nbs, a list of msites

        site.higher_rank_nbs, a list of msites

        site.eq_rank_nbs, a list of msites

        :param msties: including the bone-joint site
        :param scid: sidechain_id
        :param omol: parent omol obj
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
                try:  # I added this due to situation like RAZTAN.cif where one sc has two bone_joints
                    rank_dummy = nb.rank
                except AttributeError:
                    continue
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
                warnings.warn('W: inner loop in this side chain, scid={}'.format(scid))
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
        rankmap = [None, [bone_joint]]  # rankmap[2] is [sidejoints], rankmap[3] is [msites with rank==3]
        for rank in range(2, maxrank+1):
            rankmap.append(sorted([s for s in sc_msites if s.rank == rank], key=lambda x: len(x.higher_rank_nbs)))

        v1 = bone_joint.coords - omol.backbone.geoc
        v2 = omol.backbone.vp_fit
        avp = angle_btw(v1, v2, output='degree')
        if avp > 90.0:
            avp = 180 - avp
        return cls(sc_msites, bone_joint, scid, rankmap, angle_vp=avp, hasring=hasring)

    @property
    def maxrank(self):
        """
        :return: int, max rank of this sidechain after assigning rank to each site
        """
        return max([s.rank for s in self.msites])
