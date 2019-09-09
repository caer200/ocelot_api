import numpy as np
from copy import deepcopy
import warnings
from routines.geometry import angle_btw, Fitter, unify, get_proj_point2plane
from schema.msitelist import MSitelist
from schema.msite import MSite
from schema.ring import Ring

_coplane_cutoff = 25.0  # in degrees


class Backbone(MSitelist):

    def __init__(self, msites, lfit_linearity, pfit_vp, pfit_vq, pfit_vo, pfit_error, backbone_rings):
        """
        backbone within a omol object, again build it from omol

        :param msites: a list of msites
        :param lfit_linearity: error from linear fit
        :param pfit_vp: long axis vector from plane_fit
        :param pfit_vq: short axis vector
        :param pfit_vo: normal vector
        :param pfit_error: plane_fit error
        :param backbone_rings: a list of individual rings
        """
        for ms in msites:
            if ms.siteid == -1:
                warnings.warn('you are init a backbone with sites not in an omol obj')
        super().__init__(msites)
        self.siteids = [s.siteid for s in self.msites]  # TODO this could be in the parent obj
        self.backbone_rings = backbone_rings

        self.vo_fit = pfit_vo
        self.vp_fit = pfit_vp
        self.vq_fit = pfit_vq
        self.lfit_linearity = lfit_linearity
        self.pfit_error = pfit_error
        self.backbone_rings = backbone_rings

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        msites, linearity, vp, vq, vo, plane_fit_error, backbone_rings
        :param dict d:
        """
        msites = [MSite.from_dict(sdict) for sdict in d['msites']]
        lfit_linearity = d['linearity']
        pfit_vp = d['vp']
        pfit_vq = d['vq']
        pfit_vo = d['vo']
        pfit_error = d['plane_fit_error']
        backbone_rings = [Ring.from_dict(rdict) for rdict in d['backbone_rings']]
        return cls(msites, lfit_linearity, pfit_vp, pfit_vq, pfit_vo, pfit_error, backbone_rings)

    def as_dict(self):
        """
        keys are

        msites, linearity, vp, vq, vo, plane_fit_error, backbone_rings, n_backbone_rings, p_length, q_length, can, volume
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "p_length": self.lp,
            "q_length": self.lq,
            "backbone_rings": [r.as_dict() for r in self.backbone_rings],
            "n_backbone_rings": len(self.backbone_rings),
            "linearity": self.lfit_linearity,
            "plane_fit_error": self.pfit_error,
            "can": self.canonical_smiles,
            "volume": self.volume,
            "msites": [s.as_dict() for s in self.msites],
            "vp": self.vp_fit,
            "vq": self.vq_fit,
            "vo": self.vo_fit,
        }
        return d

    @property
    def lq(self):
        """
        maxdiff( (s.coord - ref) proj at vq )

        :return: long axis length
        """
        ref = self.backbone_rings[0].geoc
        projs = [np.dot(s.coords - ref, self.vq_fit) for s in self.msites]
        return max(projs) - min(projs)

    @property
    def lp(self):
        """
        maxdiff( (s.coord - ref) proj at vp )

        :return: short axis length
        """
        ref = self.backbone_rings[0].geoc
        projs = [np.dot(s.coords - ref, self.vp_fit) for s in self.msites]
        return max(projs) - min(projs)

    @classmethod
    def from_omol(cls, omol):
        """
        backbone is built from the largest fused ring system with connected rings sharing the same plane but
        without appending non-ring insaturated sites

        e.g.

        TIPS-fusedring1-fusedring2=o will have a backbone as fusedring1-fusedring2
        if the angle < _coplane_cutoff (25.0 degrees default)

        otherwise it gives fusedring2 since it contains more INDIVIDUAL rings than fusedring1
        """
        fused_rings_list = omol.fused_rings_list  # [[r1, r2, r3], [r5, r6], [r4]...]
        largest_fused_ring = fused_rings_list[0]
        other_fused_rings = fused_rings_list[1:]

        s_in_largest_fused_ring = []
        for r in largest_fused_ring:
            s_in_largest_fused_ring += r.msites
        # MSitelist(s_in_largest_fused_ring).to_xyz('d.xyz')
        vo, ptsmean, pfit_error = Fitter.plane_fit([s.coords for s in s_in_largest_fused_ring])

        backbone_rings = [r for r in largest_fused_ring]
        for fr in other_fused_rings:
            s_in_fr = []
            fr_nbs = []
            for r in fr:
                s_in_fr += r.msites
                for s in r:
                    fr_nbs += s.nbs_idx
            backbone_rings_site_idx = []
            for r in backbone_rings:
                backbone_rings_site_idx += [s.siteid for s in r]

            if len(set(fr_nbs).intersection(set(backbone_rings_site_idx))) > 0:
                fr_vo, fr_ptsmean, fr_pfit_error = Fitter.plane_fit([s.coords for s in s_in_fr])
                if angle_btw(fr_vo, vo, output="degree") < _coplane_cutoff or \
                        angle_btw(-fr_vo, vo, output="degree") < _coplane_cutoff:
                    backbone_rings += fr

        lfit_vp, lfit_ptsmean, lfit_linearity = Fitter.linear_fit([r.geoc for r in backbone_rings])
        lfit_vp = unify(lfit_vp)

        s_in_backbone = []
        for r in backbone_rings:
            s_in_backbone += r.msites

        # this is the plane fit, vp taken from the projection of lfit_vp
        # the fitted plane is now used as a cart coord sys with origin at ptsmean
        pfit_vo, pfit_ptsmean, pfit_error = Fitter.plane_fit([s.coords for s in s_in_backbone])
        pfit_vp = unify(get_proj_point2plane(pfit_ptsmean + lfit_vp, pfit_vo, pfit_ptsmean) - pfit_ptsmean)
        pfit_vq = unify(np.cross(pfit_vo, pfit_vp))
        # pfit_plane_params = get_plane_param(pfit_vo, pfit_ptsmean)
        return cls(s_in_backbone, lfit_linearity, pfit_vp, pfit_vq, pfit_vo, pfit_error, backbone_rings)

    def terminate(self):
        """
        basically add-H

        only terminate those have nbs different from what they had in omol
        (+1 or +2, otherwise do nothing--deepcopy only),

        this means your cif should be legit

        the H added wiill have siteid as -10

        :return: a list of msites
        """
        terminated_sites = deepcopy(self.msites)
        for siteid in self.siteids:  # siteid is the id to be terminated
            origincoords = self.get_site_byid(siteid).coords
            omol_nbs_ids = self.get_site_byid(siteid).nbs_idx
            backbon_nbs_ids = [sid for sid in omol_nbs_ids if sid in self.siteids]
            if len(omol_nbs_ids) == len(backbon_nbs_ids) + 1:
                vsh = np.zeros(3)
                for sid in backbon_nbs_ids:
                    nbcoords = self.get_site_byid(sid).coords
                    vsh += nbcoords - origincoords
                coords = -1.1 * unify(vsh) + origincoords
                hsite = MSite('H', coords, siteid=-10)
                terminated_sites.append(hsite)
            elif len(omol_nbs_ids) == len(backbon_nbs_ids) + 2:
                # TODO right now we add hydrogens vertical to the plane, it should be sp3 like
                nb1_id, nb2_id = backbon_nbs_ids[:2]
                nb1_coords = self.get_site_byid(nb1_id).coords
                nb2_coords = self.get_site_byid(nb2_id).coords
                vsh = np.cross(nb1_coords - origincoords, nb2_coords - origincoords)
                coords_1 = -1.1 * unify(vsh) + origincoords
                coords_2 = 1.1 * unify(vsh) + origincoords
                terminated_sites.append(MSite('H', coords_1, siteid=-10))
                terminated_sites.append(MSite('H', coords_2, siteid=-10))
        return terminated_sites  # TODO maybe it's better to return a backbone obj

