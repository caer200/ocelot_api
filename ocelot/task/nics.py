import warnings
from copy import deepcopy

import numpy as np
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import Site
from pymatgen.io.gaussian import GaussianInput
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

from ocelot.routines.geometry import angle_btw
from ocelot.routines.geometry import genlinepts
from ocelot.schema.conformer import MolConformer
from ocelot.schema.conformer import RingConformer
from ocelot.schema.conformer import _coordination_rule

"""
given a non-solvent omol, generate input for nics scan, parse output
the omol should be reoriented such that z axis is perpendicular to the backbone plane
"""


class NicsError(Exception):
    pass


class NICSjob:
    lgfr: [RingConformer]

    def __init__(self, omol: MolConformer, nmrscheme='GIAO'):
        self.omol = omol
        self.nmrscheme = nmrscheme
        if self.omol.is_solvent:
            NicsError('NICS does not work for solvent-like molecule')
        if self.omol.is_solvent:
            warnings.warn('W: trying to do NICS job for a solvent-like molecule')

        bone = self.omol.backbone
        lgfr = bone.rings
        geocs = [r.geoc for r in lgfr]  # nx3
        dmat = squareform(pdist(np.array(geocs)))
        maxi, maxj = np.unravel_index(dmat.argmax(), dmat.shape)  # idx for rings at two ends
        fr = sorted(lgfr, key=lambda r: np.linalg.norm(r.geoc - lgfr[maxi].geoc))
        self.lgfr = fr

    def plot_path(self):
        """
        # TODO plot the chemical structure of self.omol and show a path projection there
        :return: 
        """
        pass

    def nics_sigma_structure(self, normal_idx=0):
        """
        add hydrogen to all sites on lfr(largest fused rings), and sites that are connected to lfr with double bonds
        notice this could be different from backbone, as backbone contains other fr connected to lfr with single bond 
        (eg. BDT 3mer)


              :param normal_idx: 0, 1
        :return: a MSitelist objects
        """
        lgfr = self.lgfr

        if normal_idx == 0:
            ref_normal = lgfr[0].n1
        else:
            ref_normal = lgfr[0].n2
        normals = [r.normal_along(ref_normal) for r in lgfr]
        if any([n is None for n in normals]):
            warnings.warn('W: nics_sigma_structure failed as largest_fused_ring maybe too twisted')
            return None

        frsites = []
        for r in lgfr:
            for s in r.sites:
                if s not in frsites:
                    frsites.append(s)

        sigma_sites = deepcopy(self.omol.sites)

        terminated_sites_id = []
        added_hydrogens = []
        nbrmap = self.omol.get_nbrmap_based_on_siteid(1.3)

        for i in range(len(lgfr)):
            ringconf = lgfr[i]
            for s in ringconf.sites:
                nbs = nbrmap[s.properties['siteid']]
                unsaturation = _coordination_rule[s.species_string] - len(nbs)
                if unsaturation == 1 and s.properties['siteid'] not in terminated_sites_id:
                    terminated_sites_id.append(s.properties['siteid'])
                    hs = Site('H', 1.0 * normals[i] + s.coords, properties={'siteid': -10})
                    added_hydrogens.append(hs)
                for nb in nbs:
                    nbsite = self.omol.get_site_byid(nb)
                    nb_unsaturation = _coordination_rule[nbsite.species_string] - len(
                        nbrmap[nbsite.properties['siteid']])
                    if nb not in terminated_sites_id and nb_unsaturation == 1 and nbsite not in frsites and all(
                            nbsite not in ring.sites for ring in self.omol.rings):
                        terminated_sites_id.append(nb)
                        hs = Site('H', 1.0 * normals[i] + nbsite.coords, properties={'siteid': -10})
                        added_hydrogens.append(hs)
        return sigma_sites + added_hydrogens

    def nics_line_scan_path(self, step_size, nrings, height=1.7, normaldirection=0):
        """
        pts, n x 3 array, cart coords for path

        ring_idx, n x 1 array, ring_idx[i] represents the ring idx of the ring in which pts[i] reside

        xnumbers, the path travelled

        xticks, seg ends position on path

        :param step_size: in \AA
        :param nrings: only look at this # of rings in lfr
        :param height: default 1.7 is suggested by Stanger, Chemistry–A European Journal 20.19 (2014): 5673-5688
        :param normaldirection: 0 or 1, as it's possible to have two different paths for bent molecules
        :return: pts, ring_idx, xnumbers, xticks
        """
        fr = self.lgfr[:nrings]
        available_normals = [fr[0].n1, fr[0].n2]
        normals = [r.normal_along(available_normals[normaldirection]) for r in fr]

        pts = []
        ring_idx = []
        xnumbers = []
        xticks = [0]
        cross_point = np.zeros(3)
        for i in range(len(fr) - 1):
            ring = fr[i]
            nb_ring = fr[i + 1]
            segend1 = ring.geoc
            segend2 = nb_ring.geoc
            # nb_ring center -- cross_pt -- ring center -- prev_cross_pt -- prev ring center
            if i == 0:
                bond_centers = sorted([b.geoc for b in ring.bonds_in_ring],
                                      key=lambda bc: angle_btw(segend2 - segend1, bc - segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]],
                                                  key=lambda c: np.linalg.norm(c - 0.5 * (segend2 + segend1)))
                subpath = genlinepts(start_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(ring.geoc - start_point)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)

                subpath = genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(ring.geoc - cross_point)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)
            elif i != 0 and i != len(fr) - 2:
                prev_cross_point = cross_point
                bond_centers = sorted([b.geoc for b in ring.bonds_in_ring],
                                      key=lambda bc: angle_btw(segend2 - segend1, bc - segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]],
                                                  key=lambda c: np.linalg.norm(c - 0.5 * (segend2 + segend1)))
                subpath = genlinepts(prev_cross_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)

                subpath = genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)
            elif i == len(fr) - 2:
                prev_cross_point = cross_point
                bond_centers = sorted([b.geoc for b in ring.bonds_in_ring],
                                      key=lambda bc: angle_btw(segend2 - segend1, bc - segend1))
                cross_point, start_point = sorted([bond_centers[0], bond_centers[-1]],
                                                  key=lambda c: np.linalg.norm(c - 0.5 * (segend2 + segend1)))
                subpath = genlinepts(prev_cross_point, ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(prev_cross_point - ring.geoc)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)

                subpath = genlinepts(ring.geoc, cross_point, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - ring.geoc)]
                subpath = [p + normals[i] * height for p in subpath]
                pts += subpath
                ring_idx += [i] * len(subpath)

                subpath = genlinepts(cross_point, nb_ring.geoc, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(cross_point - nb_ring.geoc)]
                subpath = [p + normals[i + 1] * height for p in subpath]
                pts += subpath
                ring_idx += [i + 1] * len(subpath)

                bond_centers = sorted([b.geoc for b in nb_ring.bonds_in_ring],
                                      key=lambda bc: angle_btw(segend1 - segend2, bc - segend2))
                cross_point, end_point = sorted([bond_centers[0], bond_centers[-1]],
                                                key=lambda c: np.linalg.norm(c - 0.5 * (segend2 + segend1)))
                subpath = genlinepts(nb_ring.geoc, end_point, step_size)
                xnumbers += [np.linalg.norm(p - subpath[0]) + xticks[-1] for p in subpath]
                xticks += [xticks[-1] + np.linalg.norm(end_point - nb_ring.geoc)]
                subpath = [p + normals[i + 1] * height for p in subpath]
                pts += subpath
                ring_idx += [i + 1] * len(subpath)
        return pts, ring_idx, xnumbers, xticks

    @staticmethod
    def add_bqlines(stringinput, pts):
        i_lastnonempty_line = 0
        ls = stringinput.split('\n')
        for i in range(len(ls)):
            line = ls[i]
            if len(line.strip()) != 0:
                i_lastnonempty_line = i
        ls = ls[:i_lastnonempty_line + 1]
        for pt in pts:
            bqline = '{}\t{}\t{}\t{}'.format('bq', *pt)
            ls.append(bqline)
        return '\n'.join(ls) + '\n \n '

    def gen_input(
            self, step_size, nrings, maxnbq, height=1.7, normaldirection=0, charge=None, spin_multiplicity=None,
            title=None, functional='HF', basis_set='6-31G(d)', route_parameters=None, input_parameters=None,
            link0_parameters=None, dieze_tag='#P', gen_basis=None
    ):
        """
        this will return two lists of gauss input string, xticks, xnumbers, pt_idx(ring idx) for plotting

        pts are cart coords for bqs

        return None if failed

        one problem is if the number of ghost atoms cannot be too large, so the param maxnbq is introduced as
        the max number of ghost atoms that is allowed to show in one input file
        for other param docs see nics_line_scan_path and pymatgen.io.gassian
        """
        pts, pt_idx, xnumbers, xticks = self.nics_line_scan_path(step_size, nrings, height, normaldirection)
        sigma_mol_msites = self.nics_sigma_structure(normal_idx=1 - normaldirection)
        sigma_pmgmol = Molecule.from_sites(sigma_mol_msites)
        total_pmgmol = self.omol.pmgmol
        chunksize = maxnbq

        total_ginobj = GaussianInput(total_pmgmol, charge, spin_multiplicity, title, functional, basis_set,
                                     route_parameters, input_parameters, link0_parameters, dieze_tag, gen_basis)
        sigma_ginobj = GaussianInput(sigma_pmgmol, charge, spin_multiplicity, title, functional, basis_set,
                                     route_parameters, input_parameters, link0_parameters, dieze_tag, gen_basis)

        total_outs_list = []
        sigma_outs_list = []
        if len(pts) > maxnbq:
            pts_sep_list = [pts[i * chunksize:(i + 1) * chunksize] for i in
                            range((len(pts) + chunksize - 1) // chunksize)]
            for idx in range(len(pts_sep_list)):
                total_outs = total_ginobj.to_string(cart_coords=True)
                total_outs = self.add_bqlines(total_outs, pts_sep_list[idx])
                sigma_outs = sigma_ginobj.to_string(cart_coords=True)
                sigma_outs = self.add_bqlines(sigma_outs, pts_sep_list[idx])

                total_outs_list.append(total_outs)
                sigma_outs_list.append(sigma_outs)
        else:  # is this redundant?
            total_outs = total_ginobj.to_string(cart_coords=True)
            total_outs = self.add_bqlines(total_outs, pts)
            sigma_outs = sigma_ginobj.to_string(cart_coords=True)
            sigma_outs = self.add_bqlines(sigma_outs, pts)
            total_outs_list.append(total_outs)
            sigma_outs_list.append(sigma_outs)
        return sigma_outs_list, total_outs_list, xticks, xnumbers, pt_idx, pts

    @staticmethod
    def get_zzmstensor(logstring, tensor_kw='GIAO Magnetic shielding tensor'):
        """
        :param logstring: string of gauss log file
        :param tensor_kw: default 'GIAO Magnetic shielding tensor'
        :return: a list of floats, zz components of magnetic shield tensor from log string
        """
        lines = logstring.split('\n')
        ns = 0
        for i in range(len(lines)):
            if 'Charge' in lines[i] and 'Multiplicity' in lines[i]:
                cursor = i + 1
                while len(lines[cursor].strip().split()) == 4:
                    ns += 1
                    cursor += 1
                break

        tensor_zzs = []
        for i in range(len(lines)):
            if tensor_kw in lines[i]:
                for j in range(i + 1, i + 1 + ns * 5):
                    if 'BQ' in lines[j].upper():
                        tensor_zzs.append(-float(lines[j + 3].strip().split()[-1]))
                break
        return tensor_zzs

    @staticmethod
    def plt_pizz_xypath(tensor_zzs_sigma, tensor_zzs_total, xnumbers, xticks, fn):
        import matplotlib.pyplot as plt

        plt.switch_backend('agg')
        plt.rc('font', family='Times New Roman')
        plt.rcParams["font.size"] = 14
        vc = ['gray', 'black']  # gray vline is ring edge, black is ring center
        i = 1
        for xv in xticks:
            plt.axvline(xv, c=vc[i])
            i = 1 - i
        ys_total = np.array(tensor_zzs_total)
        ys_sigma = np.array(tensor_zzs_sigma)
        ys_pi = ys_total - ys_sigma
        plt.scatter(xnumbers, ys_pi, c='b', label='pi')
        plt.scatter(xnumbers, ys_total, c='k', label='total')
        plt.scatter(xnumbers, ys_sigma, c='r', label='sigma')
        plt.xlim([xnumbers[0], xnumbers[-1]])
        plt.xlabel(r'scan path ($\rm{\AA}$)')
        plt.ylabel('NICSzz (ppm)')
        plt.savefig(fn, dpi=800)
        plt.tight_layout()
        plt.clf()
