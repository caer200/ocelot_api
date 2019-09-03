import sys
import math
import copy
import itertools
import warnings
import numpy as np
import sitesop
import pymatop
import mathop
from pymatgen.core.structure import IStructure, IMolecule, Site, Molecule
from scipy.spatial.distance import pdist, squareform
# TODO remove pymatgen dependence
# TODO currently only support quasi-linear backbone


class Dimer:

    def __init__(self, slip, ref_organicmol, oth_organicmol, label):
        self.label = label
        self.slip = slip
        self.ref_organicmol = ref_organicmol
        self.oth_organicmol = oth_organicmol


class Backbone:

    def __init__(self, mol):
        self.mol = copy.deepcopy(mol)
        self.com = self.mol.center_of_mass
        self.dis_list = sorted([i.distance_from_point(self.com) for i in self.mol.sites])
        self.shortest = self.dis_list[0]
        self.longest = self.dis_list[-1]
        self.total_electrons = self.mol.nelectrons
        self.electron_density = self.total_electrons / len(self.mol.sites)
        self.site_coords = [i.coords for i in self.mol.sites]
        self.composition = self.mol.composition

        # self.uv_p, self.uv_q, self.uv_o, self.lp, self.lq = self.parse_rect_shape(self.mol)

        self.comat = sitesop.coordmat(self.mol.sites)
        self.distmat = squareform(pdist(self.comat))
        self.centers = self.ring_centers()
        self.center_line_v, self.center_line_b, self.linearity = mathop.linear_fit_3d(self.centers)
        self.uv_p = mathop.unify(self.center_line_v)
        self.uv_o = np.zeros(3)
        for i in range(len(self.mol.sites)):
            vdummy = self.comat[i] - self.comat[i+1]
            if mathop.angle_btw(self.uv_p, vdummy, output='degree') > 1:
                self.uv_o = mathop.unify(np.cross(self.uv_p, vdummy))
                break
        self.uv_q = mathop.unify(np.cross(self.uv_o, self.uv_p))
        ref = self.centers[0]
        projs = [np.dot(s.coords - ref, self.uv_q) for s in self.mol.sites]
        self.lq = max(projs) - min(projs)
        mdistance = max(self.distmat.flatten())
        mi, mj = np.where(self.distmat == mdistance)[0][0], np.where(self.distmat == mdistance)[1][0]
        v = self.comat[mi] - self.comat[mj]
        self.lp = abs(np.dot(v, self.center_line_v))
        self.mol.to('xyz', 'backbone.xyz')

    def ring_centers(self):
        centers = []
        ss = self.mol.sites
        nrnbrtable = sitesop.nrnbrmap(ss)
        loops4 = mathop.loop_finding(nrnbrtable, loopsize=4)
        loops5 = mathop.loop_finding(nrnbrtable, loopsize=5)
        loops6 = mathop.loop_finding(nrnbrtable, loopsize=6)
        if len(loops4 + loops5 + loops6) == 0:
            warnings.warn('no ring found!')
            sys.exit(1)
        for loop_idx_list in loops4 + loops5 + loops6:
            loop_coords = [ss[i].coords for i in loop_idx_list]
            center = np.average(loop_coords, axis=0)
            centers.append(center)
        return np.array(centers)

    @staticmethod
    def parse_rect_shape(mol):
        """
        assuming a rectangle shape of the mol
        p0                p3
        |
        lq
        |
        p1   --lp--       p2
        :param mol:
        :return:[long_axis_direction, short_axis_direction]
        """
        site_dis_list = sorted([[s, s.distance_from_point(mol.center_of_mass)] for s in mol.sites],
                               key=lambda x: x[1], reverse=True)
        vertices = [item[0] for item in site_dis_list[:4]]
        p0 = vertices[0]
        dis_v0 = sorted([[s, p0.distance(s)] for s in vertices], key=lambda x: x[1], reverse=True)
        p2 = dis_v0[0][0]
        p3 = dis_v0[1][0]
        p1 = dis_v0[2][0]
        v_p = mathop.unify(p3.coords - p0.coords)
        v_dummy = mathop.unify(p1.coords - p0.coords)
        v_norm = np.cross(v_p, v_dummy)
        v_q = np.cross(v_p, v_norm)
        return [v_p, v_q, v_norm, np.linalg.norm(p2.coords - p1.coords), np.linalg.norm(p0.coords - p1.coords)]


class SideChain:

    def __init__(self, mol):
        self.mol = copy.deepcopy(mol)
        self.com = self.mol.center_of_mass
        self.dis_list = sorted([i.distance_from_point(self.com) for i in self.mol.sites])
        self.dis_list_mean = np.mean(self.dis_list)
        self.dis_list_std = np.std(self.dis_list)
        self.shortest = self.dis_list[0]
        self.longest = self.dis_list[-1]
        self.total_electrons = self.mol.nelectrons
        self.electron_density = self.total_electrons / len(self.mol.sites)
        self.site_coords = [i.coords for i in self.mol.sites]
        self.composition = self.mol.composition
        self.sum_v = np.sum([i.coords - self.com for i in self.mol.sites], axis=0)
        self.sum_v_norm = np.linalg.norm(self.sum_v)


class OrganicMolecule:

    def __init__(self, mol):
        self.mol = copy.deepcopy(mol)
        self.composition = self.mol.composition
        self.bone, self.side_chain_mols = OrganicMolecule.partition(self.mol)
        self.side_chains = [SideChain(i) for i in self.side_chain_mols]

    @staticmethod
    def partition(mol):
        insaturate = []
        heteros = []
        for s in mol.sites:
            hybrid = pymatop.identify_cnhybrid(s, mol)
            if hybrid == 'csp2' or hybrid == 'nsp2':
                insaturate.append(s)
            elif s.species_string != 'C' and s.species_string != 'H':
                heteros.append(s)
        backbone = []
        for s in insaturate + heteros:
            nbs = pymatop.get_nbs_in_mol(s, mol)
            if len([nb for nb in nbs if nb in insaturate or nb in heteros]) > 1:
                backbone.append(s)
        side_groups = pymatop.group_close_sites([s for s in mol.sites if s not in backbone])
        return Backbone(Molecule.from_sites(backbone)), [Molecule.from_sites(sgm) for sgm in side_groups]


class Slip:

    def __init__(self, array):
        self.long_axis_slip = float(array[0])
        self.short_axis_slip = float(array[1])
        self.o_slip = float(array[2])
        self.long_axis_slip_ratio = float(array[3])
        self.short_axis_slip_ratio = float(array[4])
        self.h_angle = float(array[5])
        self.p_angle = float(array[6])
        self.jmol_command = array[7]


class Configuration:

    def __init__(self, filename):

        self.filename = filename
        self.structure = IStructure.from_file(self.filename)
        self.mols, self.unwrap_structure = pymatop.unwrapper(self.structure)
        self.boxmols = [self.molinbox(s) for s in self.mols]
        self.organic_mols = [OrganicMolecule(i) for i in self.mols]
        self.dimers, self.slipvis_structure = self.dimmer_analysis(self.organic_mols, self.unwrap_structure.lattice,
                                                                   maxfold=3)
        # self.sidechain = self.organic_mols[0].side_chains[0]  # assume identical
        self.slipoutstring = self.write_slip()

    @staticmethod
    def molinbox(mol):
        mol = mol.get_centered_molecule()
        coords = np.array([[j.x, j.y, j.z] for j in mol.sites], dtype='float')
        ba = max(coords[:, 0]) - min(coords[:, 0])
        bb = max(coords[:, 1]) - min(coords[:, 1])
        bc = max(coords[:, 2]) - min(coords[:, 2])
        box_structure = mol.get_boxed_structure(ba + 10, bb + 10, bc + 10, random_rotation=True, min_dist=10)
        return box_structure

    @staticmethod
    def bone_slip_data(bones, lattice, transv_fcs, cutoff):
        if cutoff == 'med':
            p_slip_ration_max = 1.5
            q_slip_ration_max = 1.5
            o_slip_max = 5.0
            o_slip_min = 0.7
        elif cutoff == 'tight':
            p_slip_ration_max = 1.0
            q_slip_ration_max = 1.0
            o_slip_max = 4.0
            o_slip_min = 0.9
        else:
            p_slip_ration_max = 5.5
            q_slip_ration_max = 5.5
            o_slip_max = 10.0
            o_slip_min = 0.1
        leni, lenk = len(bones), len(transv_fcs)
        slip_data = np.zeros((len(list(itertools.combinations_with_replacement(range(leni), 2))) * lenk, 19))
        slip_data[:] = np.NaN
        pointer = 0
        for i in range(leni):
            for j in range(i, leni):
                for k in range(lenk):
                    vslip = bones[j].com + lattice.get_cartesian_coords(transv_fcs[k]) - bones[i].com
                    if abs(np.dot(vslip, bones[i].uv_p) / bones[i].lp) < p_slip_ration_max \
                            and abs(np.dot(vslip, bones[i].uv_q) / bones[i].lq) < q_slip_ration_max \
                            and o_slip_min < abs(np.dot(vslip, bones[i].uv_o)) < o_slip_max:
                        slip_data[pointer][0] = i  # identifier for ref mol
                        slip_data[pointer][1] = j  # identifier for oth mol
                        slip_data[pointer][2] = k  # translation for oth mol
                        slip_data[pointer][3] = vslip[0]
                        slip_data[pointer][4] = vslip[1]
                        slip_data[pointer][5] = vslip[2]
                        slip_data[pointer][6] = abs(np.dot(vslip, bones[i].uv_o))  # o slip
                        slip_data[pointer][7] = abs(np.dot(vslip, bones[i].uv_p))  # p slip
                        slip_data[pointer][8] = abs(np.dot(vslip, bones[i].uv_q))  # q slip
                        slip_data[pointer][9] = abs(np.dot(vslip, bones[i].uv_p)) / bones[i].lp  # p slip ratio
                        slip_data[pointer][10] = abs(np.dot(vslip, bones[i].uv_q)) / bones[i].lq  # q slip ratio
                        slip_data[pointer][11] = math.degrees(mathop.angle_btw(bones[i].uv_o,
                                                                               bones[j].uv_o))  # h angle
                        slip_data[pointer][12] = math.degrees(mathop.angle_btw(bones[i].uv_p,
                                                                               bones[j].uv_p))  # p angle
                        slip_data[pointer][13] = bones[i].com[0]  # ref com
                        slip_data[pointer][14] = bones[i].com[1]  # ref com
                        slip_data[pointer][15] = bones[i].com[2]  # ref com
                        oth_img_com = bones[j].com + lattice.get_cartesian_coords(transv_fcs[k])
                        slip_data[pointer][16] = oth_img_com[0]  # oth com
                        slip_data[pointer][17] = oth_img_com[1]  # oth com
                        slip_data[pointer][18] = oth_img_com[2]  # oth com
                        pointer += 1
        slip_data = [i for i in slip_data[~np.isnan(slip_data).all(1)] if (i[3], i[4], i[5]) != (0., 0., 0.)]
        return sorted(slip_data, key=lambda x: abs(x[6]))

    @staticmethod
    def dimmer_analysis(organicmolsincell, lattice, maxfold):
        bones = [i.bone for i in organicmolsincell]
        transv_1d = range(-maxfold, maxfold + 1)
        transv_fcs = tuple(itertools.product(transv_1d, transv_1d, transv_1d))
        slip_data = Configuration.bone_slip_data(bones, lattice, transv_fcs, cutoff='tight')
        if len(slip_data) == 0:
            warnings.warn('slip analysis cutoff switch to med')
            slip_data = Configuration.bone_slip_data(bones, lattice, transv_fcs, cutoff='med')
        if len(slip_data) == 0:
            warnings.warn('slip analysis cutoff switch to loose')
            slip_data = Configuration.bone_slip_data(bones, lattice, transv_fcs, cutoff='loose')
        if len(slip_data) == 0:
            warnings.warn('cannot find slip even with loose criteria')
            sys.exit(1)
        dimers = []
        slip_vis_mols = organicmolsincell
        for entry in slip_data:
            i = int(entry[0])
            j = int(entry[1])
            k = int(entry[2])
            transv = lattice.get_cartesian_coords(transv_fcs[k])
            mol_ref = organicmolsincell[i]
            mol_oth = IMolecule.from_sites([Site(s.species_string, s.coords + transv)
                                            for s in organicmolsincell[j].mol.sites])

            slip_keys = ['oslip', 'pslip', 'qslip', 'pslip_ratio', 'qslip_ratio', 'h_angle', 'p_angle', 'ref_comx',
                         'ref_comy', 'ref_comz', 'oth_comx', 'oth_comy', 'oth_comz']
            slip_vals = entry[6:]
            slip = dict(zip(slip_keys, slip_vals))
            slip['jmol'] = '# draw arrow {' + ','.join([str(round(ii, 2)) for ii in [slip['ref_comx'],
                                                                                     slip['ref_comy'],
                                                                                     slip['ref_comz']]]) + '} {' \
                           + ','.join([str(round(ii, 2)) for ii in [slip['oth_comx'],
                                                                    slip['oth_comy'],
                                                                    slip['oth_comz']]]) + '} \r\n'
            dimer = Dimer(slip, mol_ref, OrganicMolecule(mol_oth), label=str(i) + '-' + str(j) + '-' + str(k))
            slip_vis_mols = slip_vis_mols + [OrganicMolecule(mol_oth)]
            dimers.append(dimer)
        slip_vis_sites = [item for sublist in [mol.mol.sites for mol in slip_vis_mols] for item in sublist]
        vis = IMolecule.from_sites(slip_vis_sites)
        return dimers, vis

    def write_slip(self):
        outstring = '# slip \r\n'
        dimers = self.dimers
        vis_sites = self.slipvis_structure.sites
        for dimer in dimers:
            slip = dimer.slip

            pslip = slip['pslip']
            qslip = slip['qslip']
            oslip = slip['oslip']
            pslip_ratio = slip['pslip_ratio']
            qslip_ratio = slip['qslip_ratio']
            h_angle = slip['h_angle']
            p_angle = slip['p_angle']
            jmol_command = slip['jmol']
            entry_items = ['%.2f' % round(i, 2) for i
                           in [pslip, qslip, oslip, pslip_ratio, qslip_ratio, h_angle, p_angle]]
            entry_items.append(jmol_command)
            entry_string = '; '.join(entry_items) + '\r\n'
            outstring = outstring + entry_string
        out = open('slip_data.txt', 'w')
        out.write(outstring)
        out.close()
        IMolecule.from_sites(vis_sites).to(fmt='xyz', filename='slip_vis.xyz')
        return outstring
