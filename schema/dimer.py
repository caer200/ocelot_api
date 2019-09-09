import numpy as np
from scipy.spatial.distance import cdist
from pymatgen.core.sites import Site
import sys
from pymatgen.core.structure import Molecule
from routines.geometry import norm, angle_btw, get_proj_point2plane, coord_transform, alpha_shape
from shapely.geometry import Polygon
from schema.omol import OMol
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class Dimer:

    def __init__(self, omol_ref, omol_var, label=""):
        """
        basically 2 omols

        :param omol_ref: first omol obj
        :param omol_var: second omol obj
        :param str label: mainly used to distinguish
        :var vslip: slip vector in cart
        :var vslipnorm: normalized vslip
        :var pslip: projection of vslip along vp
        :var qslip: projection of vslip along vq
        :var oslip: projection of vslip along vo
        :var pangle: angle btw vps
        :var qangle: angle btw vps
        :var oangle: angle btw vps, always acute
        :var jmol: jmol draw arrow string in console
        """
        self.omol_ref = omol_ref
        self.omol_var = omol_var
        self.label = label

        ref_bone = self.omol_ref.backbone
        var_bone = self.omol_var.backbone

        self.vslip = var_bone.geoc - ref_bone.geoc
        self.vslipnorm = norm(self.vslip)
        self.pslip = abs(self.vslip @ ref_bone.vp_fit)
        self.qslip = abs(self.vslip @ ref_bone.vq_fit)
        self.oslip = abs(self.vslip @ ref_bone.vo_fit)
        self.pangle = angle_btw(ref_bone.vp_fit, var_bone.vp_fit, output='degree')
        self.qangle = angle_btw(ref_bone.vq_fit, var_bone.vq_fit, output='degree')
        self.oangle = angle_btw(ref_bone.vo_fit, var_bone.vo_fit, output='degree')
        if self.oangle > 90:
            self.oangle = 180 - self.oangle

        self.jmol = "draw arrow {{ {:.4f}, {:.4f}, {:.4f} }} {{ {:.4f}, {:.4f}, {:.4f} }}".format(*ref_bone.geoc,
                                                                                                  *var_bone.geoc)

    def as_dict(self):
        """
        keys are

        vslipnorm, pslip, qslip, oslip, pangle, qangle, oangle, jmol, label, omol_ref, omol_var
        """
        d = {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
        }
        d['vslipnorm'] = self.vslipnorm
        d['pslip'] = self.pslip
        d['qslip'] = self.qslip
        d['oslip'] = self.oslip
        d['pangle'] = self.pangle
        d['qangle'] = self.qangle
        d['oangle'] = self.oangle
        d['jmol'] = self.jmol
        d['label'] = self.label
        d['omol_ref'] = self.omol_ref.as_dict()
        d['omol_var'] = self.omol_var.as_dict()
        return d

    @classmethod
    def from_dict(cls, d):
        """
        keys are

        omol_ref, omol_var, label
        """
        omol_ref = OMol.from_dict(d['omol_ref'])
        omol_var = OMol.from_dict(d['omol_var'])
        label = d['label']
        return cls(omol_ref, omol_var, label)

    def to_xyz(self, fn):
        sites = [s.to_pymatgen_site() for s in self.omol_ref.sites]
        sites += [s.to_pymatgen_site() for s in self.omol_var.sites]
        sites += [Site('La', self.omol_ref.geoc)]
        Molecule.from_sites(sites).to('xyz', fn)

    @property
    def is_identical(self):
        """
        whether two omols are identical, based on norm(vslip) < 1e-5
        """
        return np.linalg.norm(self.vslip) < 1e-5

    @property
    def is_not_identical(self):
        return not self.is_identical

    @property
    def is_close(self, cutoff=5.5):
        """
        use to identify whether this dimer can have minimum wf overlap

        this should be called is_not_faraway...

        :param cutoff: minbonedist less than which will be considered close
        :return: bool
        """
        return self.minbonedist < cutoff

    # @property
    # def stack_type(self, cutoff=1.0):
    #     if self.oangle < 30.0 * cutoff:
    #         return 'face_to_face'
    #     else:
    #         return 'face_to_edge'

    @property
    def minbonedist(self):
        """
        :return: minimum dist between sites on different bones
        """
        distmat = cdist(self.omol_ref.backbone.coordmat, self.omol_var.backbone.coordmat)
        return np.min(distmat)

    def plt_bone_overlap(self, algo='concave', output='bone_overlap.eps'):
        """
        plot a 2d graph of how backbones overlap

        using concave or convex hull

        :param algo: concave/convex
        :param output: output filename
        """
        ref_o = self.omol_ref.backbone.vo_fit
        ref_p = self.omol_ref.backbone.vp_fit
        ref_q = self.omol_ref.backbone.vq_fit
        origin = self.omol_ref.backbone.ptsmean

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.omol_ref.backbone.sites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.omol_var.backbone.sites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            sys.exit('bone_overlap receive a wrong algo spec')
        x, y = var_hull.exterior.coords.xy
        points = np.array([x, y]).T
        fig, ax = plt.subplots(1)
        polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='r', facecolor='r')
        ax.add_patch(polygon_shape)
        x, y = ref_hull.exterior.coords.xy
        points = np.array([x, y]).T
        polygon_shape = patches.Polygon(points, linewidth=1, edgecolor='b', facecolor='b', label='ref')
        ax.add_patch(polygon_shape)
        ax.axis('auto')
        fig.savefig(output)

    def mol_overlap(self, algo='concave'):
        """
        project var mol onto the plane of ref mol

        :param algo: concave/convex
        :return: area of the overlap, ref omol area, var omol area
        """
        ref_o = self.omol_ref.backbone.vo_fit
        ref_p = self.omol_ref.backbone.vp_fit
        ref_q = self.omol_ref.backbone.vq_fit
        origin = self.omol_ref.backbone.ptsmean

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.omol_ref.sites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.omol_var.sites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            sys.exit('bone_overlap receive a wrong algo spec')
        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area))

    def bone_overlap(self, algo='concave'):
        """
        project var backbone onto the plane of ref backbone

        :param algo: concave/convex
        :return: area of the overlap, ref omol area, var omol area
        """
        ref_o = self.omol_ref.backbone.vo_fit
        ref_p = self.omol_ref.backbone.vp_fit
        ref_q = self.omol_ref.backbone.vq_fit
        # origin = self.omol_ref.backbone.ptsmean
        origin = self.omol_ref.backbone.geoc

        ref_proj_pts = [get_proj_point2plane(rs.coords, ref_o, origin) for rs in self.omol_ref.backbone.msites]
        var_proj_pts = [get_proj_point2plane(vs.coords, ref_o, origin) for vs in self.omol_var.backbone.msites]

        # convert to 2d points on the plane
        ref_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in ref_proj_pts]
        var_2dpts = [coord_transform(ref_p, ref_q, ref_o, p3d)[:2] for p3d in var_proj_pts]

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        elif algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull
        else:
            sys.exit('E: bone_overlap receive a wrong algo spec')

        mindist2d = ref_hull.distance(var_hull)

        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area),
                mindist2d)
