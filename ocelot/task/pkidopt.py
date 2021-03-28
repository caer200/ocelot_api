from collections import OrderedDict
from copy import deepcopy

from pymatgen.core.structure import Lattice
from pymatgen.core.structure import Site
from scipy.spatial.distance import cdist
from scipy.spatial.distance import np
from shapely.geometry import Polygon

from ocelot.routines.disparser import DisParser
from ocelot.routines.geometry import alpha_shape
from ocelot.routines.geometry import angle_btw
from ocelot.routines.geometry import coord_transform
from ocelot.routines.geometry import get_proj_point2plane
from ocelot.routines.geometry import norm
from ocelot.schema.configuration import MolConformer
from ocelot.schema.configuration import Molecule
from ocelot.schema.configuration import Structure
from ocelot.schema.conformer import BasicConformer, BoneConformer, RingConformer

"""
one type of molecule in the unit cell
"""


class ModelBone:

    def __init__(self, pts: np.ndarray, vp, vq, vo):
        self.vp = vp
        self.vq = vq
        self.vo = vo
        self.pts = pts

    @property
    def geoc(self):
        return np.mean(self.pts, axis=0)

    def __len__(self):
        return len(self.pts)

    def __getitem__(self, item: int):
        return self.pts[item]

    def __iter__(self):
        return self.pts.__iter__()

    def mindist(self, other):
        distmat = cdist(self.pts, other.pts)
        return np.min(distmat)

    def to_xyz(self, fn):
        sites = []
        for pt in self:
            sites.append(Site('C', pt))
        Molecule.from_sites(sites).to('xyz', fn)


class ModelBoneDimer:

    def __init__(self, ref: ModelBone, var: ModelBone, label=''):
        self.ref = ref
        self.var = var
        self.label = label

        self.pslip = abs(self.vslip @ ref.vp)
        self.qslip = abs(self.vslip @ ref.vq)
        self.oslip = abs(self.vslip @ ref.vo)
        self.pangle = angle_btw(ref.vp, var.vp, output='degree')
        self.qangle = angle_btw(ref.vq, var.vq, output='degree')
        self.oangle = angle_btw(ref.vo, var.vo, output='degree')
        if self.oangle > 90:
            self.oangle = 180 - self.oangle

    @property
    def vslip(self):
        return self.var.geoc - self.ref.geoc

    @property
    def vslipnorm(self):
        return norm(self.vslip)

    @property
    def mindist(self):

        return self.ref.mindist(self.var)

    def __repr__(self):
        return 'from {} to {}'.format(self.ref.geoc, self.var.geoc)

    def to_xyz(self, fn):
        sites = []
        for pt in self.ref.pts:
            sites.append(Site('C', pt))
        for pt in self.var.pts:
            sites.append(Site('C', pt))
        Molecule.from_sites(sites).to('xyz', fn)

    @property
    def is_identical(self):
        """
        whether two omols are identical, based on norm(vslip) < 1e-5
        """
        return self.vslipnorm < 1e-5

    @property
    def is_close(self, cutoff=6.5):
        """
        use to identify whether this dimer can have minimum wf overlap, ONLY consider bone distance

        this should be called is_not_faraway...

        :param cutoff: minbonedist less than which will be considered close
        :return: bool
        """
        return self.mindist < cutoff

    def overlap(self, algo='concave'):
        """
        project var backbone onto the plane of ref backbone
        as there's the problem of alpha value in concave hull generation, maybe I should set default to convex, see hull_test

        :param algo: concave/convex
        :return: area of the overlap, ref omol area, var omol area
        """
        origin = self.ref.geoc
        ref_p = self.ref.vp
        ref_q = self.ref.vq
        ref_o = self.ref.vo

        ref = self.ref
        var = self.var

        ref_2dpts = []
        var_2dpts = []
        for i in range(len(ref)):
            ref_proj = get_proj_point2plane(ref[i], ref_o, origin)
            var_proj = get_proj_point2plane(var[i], ref_o, origin)
            ref_2d = coord_transform(ref_p, ref_q, ref_o, ref_proj)[:2]
            var_2d = coord_transform(ref_p, ref_q, ref_o, var_proj)[:2]
            ref_2dpts.append(ref_2d)
            var_2dpts.append(var_2d)

        if algo == 'concave':
            ref_hull, ref_edge_points = alpha_shape(ref_2dpts)
            var_hull, var_edge_points = alpha_shape(var_2dpts)
        else:  # algo == 'convex':
            ref_hull = Polygon(ref_2dpts).convex_hull
            var_hull = Polygon(var_2dpts).convex_hull

        mindist2d = ref_hull.distance(var_hull)

        return (ref_hull.intersection(var_hull).area,
                float(ref_hull.area),
                float(var_hull.area),
                mindist2d)


class PackingIdentifier:
    def __init__(self, bone_structure: Structure, model_bones: [ModelBone]):
        self.structure = bone_structure
        self.latt = self.structure.lattice
        self.model_bones = model_bones
        self.z = len(self.model_bones)

    def identify_heuristic(self):
        """
        return a dictionary, keys are

        n_close_azm_and_parallel,
        n_close_azm_and_notparallel,
        n_close_vertical_and_parallel,
        n_close_vertical_and_notparallel,
        n_parallel_and_overlap,
        n_notparallel_and_overlap, packing

        these keys are defined based on:
        close_vertical: d.oslip within [1.5, 4]
        parallel: d.oangle < 15 deg
        close_azm: d.oslip <= 1.5
        overlap: overlap > 1e-5

        variables used:
        1. is_not_identical: slip vector norm >= 1e-5
        2. is_close: minbonedist < 5.5
        3. overlap: project sites onto the plane defined by ovector of ref_mol and ref_mol.geoc,
            then get overlap of concave/convex hulls
        4. mindist2d: min distance between two hulls
        5. d.oslip: vertical slip
        6. d.oangle: angle between o_vector


        workflow:
        1. from bone config get all bone dimers, maxfold=2
        2. for each molecule (a) in the unit cell
            1. get first 27*z dimers, sorted by vslipnorm, whose ref_mol is (a)
            2. if identical or not_close, do nothing
            3. get values for the keys
            4. apply classification based on the keys
        :return:
        """
        bc_dimers, bc_dimers_transv_fcs = get_modeldimers_array(self.model_bones, self.latt, maxfold=2)

        report = OrderedDict()
        for i in range(self.z):
            n_close_azm_and_parallel = 0  # min distance between backbone proj<2, 1.5>oslip
            n_close_azm_and_notparallel = 0  # min distance between backbone proj<2, 1.5>oslip
            n_close_vertical_and_parallel = 0  # min distance between backbone proj<2, 4>oslip>1.5
            n_close_vertical_and_notparallel = 0
            n_parallel_and_overlap = 0
            n_notparallel_and_overlap = 0
            dimers_ref_i = bc_dimers[i].flatten()

            # dimers_ref_i = sorted(dimers_ref_i, key=lambda x: x.vslipnorm)[:27 * self.z]  # 3x3x3
            # 10.1039/c0ce00044b says 14 is enough
            dimers_ref_i = sorted(dimers_ref_i, key=lambda x: x.vslipnorm)[:16]

            sites = []
            bcd: ModelBoneDimer
            for bcd in dimers_ref_i:
                for pt in bcd.var.pts:
                    sites.append(Site('C', pt))
            for pt in dimers_ref_i[0].ref:
                sites.append(Site('N', pt))
            # Molecule.from_sites(sites).to('xyz', 'dccollection-{}.xyz'.format(i))

            dimers_ref_i = [d for d in dimers_ref_i if not d.is_identical and d.is_close]

            d: ModelBoneDimer
            for d in dimers_ref_i:
                # if d.is_not_identical and d.is_close:
                overlap, refboneproj, varboneproj, mindist2d = d.overlap(algo='convex')
                # print(d.minbonedist)
                # print(mindist2d)
                # debug += 1
                # d.to_xyz('dimer_{}.xyz'.format(debug))
                # overlap, refboneproj, varboneproj = d.mol_overlap()
                if 1e-5 < mindist2d < 2.5:  # exclude overlap as dist=0.0 that case
                    if 4 > d.oslip > 1.5:
                        if d.oangle < 15.0:
                            n_close_vertical_and_parallel += 1
                        else:
                            n_close_vertical_and_notparallel += 1
                    elif d.oslip <= 1.5:
                        if d.oangle < 15.0:
                            n_close_azm_and_parallel += 1
                        else:
                            n_close_azm_and_notparallel += 1
                if overlap > 1e-5:
                    if d.oangle < 15.0:
                        n_parallel_and_overlap += 1
                    else:
                        n_notparallel_and_overlap += 1

            if n_parallel_and_overlap > 2:
                packing = 'brickwork'

            elif n_close_vertical_and_parallel > 2 or n_close_azm_and_parallel > 2 or n_close_vertical_and_parallel + n_close_azm_and_parallel + n_parallel_and_overlap > 2:
                packing = 'edge_brickwork'

            elif n_parallel_and_overlap == 2:
                packing = 'slipped_stack'

            elif n_close_vertical_and_parallel == 2 or n_parallel_and_overlap + n_close_vertical_and_parallel + n_close_azm_and_parallel == 2:
                packing = 'edge_slipped_stack'

            elif n_notparallel_and_overlap >= 1:
                packing = 'herringbone'

            elif n_close_vertical_and_notparallel >= 1 or n_close_azm_and_notparallel >= 1:
                packing = 'edge_herringbone'

            elif n_parallel_and_overlap == 1 or n_close_vertical_and_parallel == 1:
                packing = 'dimerized'

            elif n_notparallel_and_overlap == n_parallel_and_overlap == n_close_azm_and_notparallel == n_close_azm_and_parallel == n_close_vertical_and_parallel == n_close_vertical_and_notparallel == 0:
                packing = 'sparse'

            else:
                packing = 'other'

            report[i] = OrderedDict(
                n_close_azm_and_parallel=n_close_azm_and_parallel,  # min distance between backbone proj<2, 1.5>oslip
                n_close_azm_and_notparallel=n_close_azm_and_notparallel,
                # min distance between backbone proj<2, 1.5>oslip
                n_close_vertical_and_parallel=n_close_vertical_and_parallel,
                # min distance between backbone proj<2, 4>oslip>1.5
                n_close_vertical_and_notparallel=n_close_vertical_and_notparallel,
                n_parallel_and_overlap=n_parallel_and_overlap,
                n_notparallel_and_overlap=n_notparallel_and_overlap,
                packing=packing,
            )

        return report


def partition_basicconformer(conformer: BasicConformer):
    molgraph = conformer.to_graph('siteid', 'MolGraph')
    rings = []
    for r in molgraph.lgfr:
        ring = RingConformer.from_siteids(r, conformer.sites, copy=False)
        rings.append(ring)
    avgn1, avgn2 = RingConformer.get_avg_norms(rings)

    def coplane_check(ring_siteids, tol=30):
        ringconformer = RingConformer.from_siteids(ring_siteids, conformer.sites, False)
        coplane = ringconformer.iscoplane_with_norm(avgn1, tol, 'degree')
        return coplane
    bg, scgs = molgraph.partition_to_bone_frags('lgcr', additional_criteria=coplane_check)
    bone_conformer = BoneConformer.from_siteids(bg.graph.nodes, conformer.sites, graph=bg, copy=False)
    return bone_conformer

def get_bone_structure_cheaper(config: Structure, mols: [Molecule]):
    """
    cheaper method to get bone structure, relies on siteid prop assigned during to_configs method in DP
    """
    conf_ids = [s.properties['siteid'] for s in config]
    backboneid_in_config = []
    model_bones = []
    for mol in mols:
        clean_siteids = [s.properties['siteid'] for s in mol if s.properties['siteid'] in conf_ids]
        basic_conformer = BasicConformer.from_siteids(clean_siteids, mol.sites)
        b = partition_basicconformer(basic_conformer)
        # clean_mc = MolConformer.from_siteids(clean_siteids, mol.sites)
        # b = clean_mc.backbone
        model_bone = ModelBone(b.cart_coords, b.pfit_vp, b.pfit_vq, b.pfit_vo)
        # bone_ids = clean_mc.backbone.siteids
        bone_ids = b.siteids
        backboneid_in_config += bone_ids
        model_bones.append(model_bone)
    bc = Structure.from_sites([s for s in config if s.properties['siteid'] in backboneid_in_config])
    return bc, model_bones


def get_modeldimers_array(bones: [ModelBone], lattice: Lattice, maxfold=2):
    """
    :param maxfold: translation vector in fc can be [h, h, h] where maxfold <= h <= maxfold
    :return: dimers_array, z x z x n array, dimers[i][j][k] is the dimer of omols[i], omols[j] with translation vector as transv_fcs[k]
             transv_fcs
    """
    z = len(bones)
    transv_1d = list(range(-maxfold, maxfold + 1))
    transv_fcs = np.array(np.meshgrid(transv_1d, transv_1d, transv_1d)).T.reshape(-1, 3)
    # symmetry dimers[i][j][transv_fcs[k]] = dimers[j][i][-transv_fcs[k]]
    dimers = np.empty((z, z, len(transv_fcs)), dtype=ModelBoneDimer)
    used_transv_fcs = transv_fcs

    for i in range(z):
        ref = deepcopy(bones[i])
        for j in range(z):
            var = deepcopy(bones[j])
            for k in range(len(used_transv_fcs)):
                transv = lattice.get_cartesian_coords(used_transv_fcs[k])
                var_k: ModelBone = deepcopy(var)
                for h in range(len(var_k)):
                    var_k.pts[h] += transv
                dimer_ijk = ModelBoneDimer(ref, var_k, label="{}_{}_{}".format(i, j, k))
                dimers[i][j][k] = dimer_ijk
    return dimers, transv_fcs


class PkidError(Exception):pass

def pkid_ciffile(ciffile):
    dp = DisParser.from_ciffile(ciffile)
    pstructure, unwrap_str, mols, conf_infos= dp.to_configs(write_files=False)  # writes conf_x.cif
    if len(set([m.composition for m in mols])) != 1:  # TODO ideally one should check at the graph level
        raise PkidError("not one type of molecule!")
    config = conf_infos[0][0]
    # mols still contains disordered sites!
    import time
    ts1 = time.time()
    bs, model_bones = get_bone_structure_cheaper(config, mols)
    ts2 = time.time()
    print('boneconfig', ts2 - ts1)
    pid = PackingIdentifier(bs, model_bones)
    pdata = pid.identify_heuristic()
    ts3 = time.time()
    print('idheuristic', ts3-ts2)
    return pdata

