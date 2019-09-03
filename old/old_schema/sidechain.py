import sys
import warnings
from api.routines.geometry import angle_btw, norm
from api.schema.msitelist import MSitelist


class Sidechain(MSitelist):

    def __init__(self, sites, omol, scid=-1):
        """
        the first site in msites is bone_joint, the obj has msites started from side_joint
        :param sites:
        :param omol:
        :param scid:
        """
        if len(sites) < 2:
            sys.exit('Sidechain obj must be init with at least 2 msites!')
        super().__init__(sites[1:])
        self.scid = scid
        self.omol = omol
        self.bone_joint = sites[0]
        self.side_joint = self.sites[0]
        self.hasring = False  # inner loop

        for i in range(len(self.sites)):
            self.sites[i].path_to_bj = self.omol.get_shortest_path(self.sites[i].siteid, self.bone_joint.siteid,
                                                                   self.omol.nbrmap)
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
                self.hasring = True
                warnings.warn('inner loop in this side chain, id={}'.format(self.scid))
                # sys.exit('inner loop in this side chain, id=' + str(self.scid))
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

        self.rankmap = [None, None]  # rankmap[2] is [sidejoint], rankmap[3] is [msites with rank==3]
        for rank in range(2, self.maxrank):
            self.rankmap.append(sorted([s for s in self.sites if s.rank == rank], key=lambda x: len(x.higher_rank_nbs)))

    def get_sphere_diameter(self):
        ranksorted_sites = sorted(self.sites, key=lambda x: x.rank)
        for i in range(len(ranksorted_sites)):
            if len(ranksorted_sites[i].higher_rank_nbs) == 3:
                node = ranksorted_sites[i]
                break
        try:
            r = max([norm(s.coords - node.coords) for s in self.sites if s.rank > node.rank])
        except:
            r = 0
        return 2*r


    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d['scif'] = self.scid
        d['can'] = self.canonical_smiles
        d['hasring'] = self.hasring
        d['is_hydrogen'] = self.ishydrogen
        d['max_rank'] = self.maxrank
        d['volume'] = self.volume
        d['msites'] = [s.as_dict() for s in self.sites]
        return d


    @property
    def ishydrogen(self):
        return len(self.sites) == 1 and self.sites[0].symbol == 'H'

    @property
    def angle_vp(self):
        v1 = self.bone_joint.coords - self.omol.backbone.geoc
        v2 = self.omol.backbone.vp
        avp = angle_btw(v1, v2, output='degree')
        # if avp > 90.0:
        #     self.avp = 180 - self.avp
        return avp

    @property
    def maxrank(self):
        return max([s.rank for s in self.sites])

    # TODO add similarity and morphology descriptors
