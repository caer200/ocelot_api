import math
from copy import deepcopy
from collections import OrderedDict
import itertools
import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen.symmetry.groups import SpaceGroup, SYMM_DATA
import re
import warnings
from itertools import groupby
from pymatgen.io.cif import CifFile, str2float, _get_cod_data, sub_spgrp
from pymatgen.core.structure import PeriodicSite, Structure, Lattice, lattice_points_in_supercell
from ocelot.routines.pbc import PBCparser

"""
DisParser: parse cif file into a list of configurations with no disorder

this is used to parse disordered cif file to configurations of a supercell

we assume the cif file contains following disorder-related fields
    _atom_site_occupancy 
    _atom_site_disorder_group 

if above fields do not exist, it is still possible the cif file contains disorder where
    1. sites of alternative configuration are labeled by special character (? or ') as appendix, e.g. ALOVOO
    2. sites of alternative configuration are labeled by capital letters as appendix, e.g. ASIXEH
    3. disordered sites are simply ignored in the cif file s.t. the molecule is not complete, e.g. ABEGET01
    4. disordered sites are simply ignored in the cif file, but the molecule looks legit, e.g. AGAJUO
    5. disordered sites are not labeled, e.g. ANOPEA

I cannot find a way to process 2., 3., 4., 5. without checking deposition records or chemistry analysis, so this parser
will ignore them

I will deal with 1. by artificially set (average) _atom_site_occupancy and _atom_site_disorder_group based on special character 
appendix, this should be done intentionally

configuration will be written with keys:
['_cell_length_a'],
['_cell_length_b'],
['_cell_length_c'],
['_cell_angle_alpha'],
['_cell_angle_beta'],
['_cell_angle_gamma'],
['_space_group_symop_operation_xyz'] or other symm labels
['_atom_site_label'],
['_atom_site_type_symbol'],
['_atom_site_fract_x'],
['_atom_site_fract_y'],
['_atom_site_fract_z'],
"""

possible_symm_labels = [
    "_symmetry_equiv_pos_as_xyz",
    "_symmetry_equiv_pos_as_xyz_",
    "_space_group_symop_operation_xyz",
    "_space_group_symop_operation_xyz_",
    "_symmetry_space_group_name_H-M",
    "_symmetry_space_group_name_H_M",
    "_symmetry_space_group_name_H-M_",
    "_symmetry_space_group_name_H_M_",
    "_space_group_name_Hall",
    "_space_group_name_Hall_",
    "_space_group_name_H-M_alt",
    "_space_group_name_H-M_alt_",
    "_symmetry_space_group_name_hall",
    "_symmetry_space_group_name_hall_",
    "_symmetry_space_group_name_h-m",
    "_symmetry_space_group_name_h-m_",
    "_space_group_IT_number",
    "_space_group_IT_number_",
    "_symmetry_Int_Tables_number",
    "_symmetry_Int_Tables_number_"]

latt_labels = [
    '_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma',
]

chemistry_labels = [
    '_cell_formula_units_Z', '_chemical_formula_moiety', '_chemical_formula_sum'
]

coord_labels = [
    '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y',
    '_atom_site_fract_z'
]

disorder_labels = [
    '_atom_site_occupancy', '_atom_site_disorder_group',
]

space_groups = {sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()}  # type: ignore
space_groups.update({sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()})  # type: ignore


def get_pmg_dict(cifstring):
    """
    use pmg dict to parse cifstring, only deal with one structure per file

    :param cifstring:
    :return:
    """
    cifdata = CifFile.from_string(cifstring).data
    idnetifiers = list(cifdata.keys())
    if len(idnetifiers) > 1:
        warnings.warn('W: find more than 1 structures in this cif file!')
    elif len(idnetifiers) == 0:
        warnings.warn('W: no structure found by pymatgen parser!')
    identifier = idnetifiers[0]
    pymatgen_dict = list(cifdata.items())[0][1].data
    return identifier, pymatgen_dict


class DisParser:  # chaos parser sounds cooler?

    def __init__(self, cifstring):
        self.cifstring = cifstring
        self.identifier, self.data = get_pmg_dict(self.cifstring)

    @classmethod
    def from_ciffile(cls, fn):
        with open(fn, 'r') as f:
            s = f.read()
        return cls(s)

    def classify(self):
        """
        one cif file belongs to one of the following 3 categories
        1 has _atom_site_occupancy AND _atom_site_disorder_group
        2 some of the labels cannot match r"^\w+$"
        3 else

        :return:
        """
        if '_atom_site_occupancy' in self.data.keys() and '_atom_site_disorder_group' in self.data.keys():
            return '1'
        else:
            for label in self.data['_atom_site_label']:
                if not re.search(r"^\w+$", label):
                    return '2'
            return '3'

    @staticmethod
    def class2_to_class1(cifdata):
        coord_data = map(list, zip(cifdata['_atom_site_fract_x'],
                                   cifdata['_atom_site_fract_y'],
                                   cifdata['_atom_site_fract_z'],
                                   cifdata['_atom_site_type_symbol']
                                   ))
        coord_data = list(coord_data)
        coord_data = OrderedDict(zip(cifdata['_atom_site_label'], coord_data))
        group_by_word = [list(v) for l, v in
                         groupby(sorted(coord_data.keys(), key=lambda x: re.search(r"\w+", x).group(0)),
                                 key=lambda x: re.search(r"\w+", x).group(0))]
        for i in range(len(group_by_word)):
            labels = sorted(group_by_word[i])  # this works if all disorder charaters are added as appendix
            if len(labels) == 1:
                coord_data[labels[0]].append('1')
                coord_data[labels[0]].append('.')  # convention for no disorder
            else:
                occu = 1.0 / len(labels)
                idisg = 1  # convention start from 1
                for lab in labels:
                    coord_data[lab].append(str(occu))
                    coord_data[lab].append(str(idisg))
                    idisg += 1
        newdata = OrderedDict()
        for ciflabel in possible_symm_labels + latt_labels + chemistry_labels:
            if ciflabel in cifdata.keys():
                newdata[ciflabel] = cifdata[ciflabel]
        labs = list(coord_data.keys())
        xs = []
        ys = []
        zs = []
        symbols = []
        occus = []
        idisgs = []
        for lab in labs:
            x, y, z, symb, occu, idisg = coord_data[lab]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            symbols.append(symb)
            occus.append(occu)
            idisgs.append(idisg)
        newdata['_atom_site_label'] = labs
        newdata['_atom_site_type_symbol'] = symbols
        newdata['_atom_site_fract_x'] = xs
        newdata['_atom_site_fract_y'] = ys
        newdata['_atom_site_fract_z'] = zs
        newdata['_atom_site_occupancy'] = occus
        newdata['_atom_site_disorder_group'] = idisgs
        return newdata

    @staticmethod
    def class3_to_class1(cifdata):
        coord_data = map(list, zip(cifdata['_atom_site_fract_x'],
                                   cifdata['_atom_site_fract_y'],
                                   cifdata['_atom_site_fract_z'],
                                   cifdata['_atom_site_type_symbol']
                                   ))
        coord_data = list(coord_data)
        coord_data = OrderedDict(zip(cifdata['_atom_site_label'], coord_data))
        newdata = OrderedDict()
        for ciflabel in possible_symm_labels + latt_labels + chemistry_labels:
            if ciflabel in cifdata.keys():
                newdata[ciflabel] = cifdata[ciflabel]
        labs = list(coord_data.keys())
        xs = []
        ys = []
        zs = []
        symbols = []
        occus = []
        idisgs = []
        for lab in labs:
            x, y, z, symb = coord_data[lab]
            xs.append(x)
            ys.append(y)
            zs.append(z)
            symbols.append(symb)
            occus.append('1')
            idisgs.append('.')
        newdata['_atom_site_label'] = labs
        newdata['_atom_site_type_symbol'] = symbols
        newdata['_atom_site_fract_x'] = xs
        newdata['_atom_site_fract_y'] = ys
        newdata['_atom_site_fract_z'] = zs
        newdata['_atom_site_occupancy'] = occus
        newdata['_atom_site_disorder_group'] = idisgs
        return newdata

    @staticmethod
    def class1_to_configs(cifdata, write_files=False, scaling_mat=(1, 1, 1)):
        psites, disunit_pairs, n_symmops = get_psites(cifdata)
        pstructure = Structure.from_sites(psites)
        mols, unwrap_str_sorted, unwrap_pblock_list = PBCparser.unwrap(pstructure)
        sc, n_unitcell = build_supercell_full_disorder(unwrap_str_sorted, scaling_mat)
        conf_ins = gen_instructions(disunit_pairs, n_symmops, n_unitcell)
        iconf = 0
        confs = []
        for confin in conf_ins:
            conf, conf_occu = dissc_to_config(sc, disunit_pairs, confin)
            confs.append([conf, conf_occu])
            if write_files:
                conf.to('cif', 'conf_{}.cif'.format(iconf))
            iconf += 1
        if write_files:
            pstructure.to('cif', 'confgen_ps.cif')
            unwrap_str_sorted.to('cif', 'confgen_unwrap.cif')
        return sorted(confs, key=lambda x: x[1], reverse=True)

    def to_configs(self, write_files=False, scaling_mat=(1, 1, 1)):
        classification = self.classify()
        if classification == '1':
            cifdata = self.data
        elif classification == '2':
            cifdata = self.class2_to_class1(self.data)
        elif classification == '3':
            cifdata = self.class3_to_class1(self.data)
        else:
            return None
        return self.class1_to_configs(cifdata, write_files, scaling_mat)


def apply_symmop(psites, ops):
    new_psites = []
    for ps in psites:
        iasym = 0
        for op in ops:
            new_coord = op.operate(ps.frac_coords)
            new_coord = np.array([i - math.floor(i) for i in new_coord])
            new_properties = deepcopy(ps.properties)
            new_properties['iasym'] = iasym
            new_ps = PeriodicSite(ps.species_string, new_coord, ps.lattice, properties=deepcopy(new_properties))
            new_psites.append(new_ps)
            iasym += 1
    return new_psites


def get_psites(pmg_dict):
    """
    get pmg_psites from cif data

    :param pmg_dict:
    :return: fin_psites, a list of psites with 'occu', 'disg', 'label' in properties
            disunit_pairs, a list of disorderunit pairs in one aym unit, e.g. [[u1, u1'], [u2, u2']] where u1 xor u1'
            len(symmops), number of aym units in a cell
    """
    labels = pmg_dict['_atom_site_label']
    symbols = pmg_dict['_atom_site_type_symbol']
    xs = [str2float(j) for j in pmg_dict['_atom_site_fract_x']]
    ys = [str2float(j) for j in pmg_dict['_atom_site_fract_y']]
    zs = [str2float(j) for j in pmg_dict['_atom_site_fract_z']]
    occus = [str2float(j) for j in pmg_dict['_atom_site_occupancy']]
    disgs = pmg_dict['_atom_site_disorder_group']
    psites = []
    length_strings = ("a", "b", "c")
    angle_strings = ("alpha", "beta", "gamma")
    lengths = [str2float(pmg_dict["_cell_length_" + i]) for i in length_strings]
    angles = [str2float(pmg_dict["_cell_angle_" + i]) for i in angle_strings]
    lattice = Lattice.from_parameters(*lengths, *angles)
    for i in range(len(labels)):
        ps = PeriodicSite(symbols[i], [xs[i], ys[i], zs[i]], lattice,
                          properties={'occu': occus[i], 'disg': disgs[i], 'label': labels[i]})
        psites.append(ps)
    disunit_pairs = DisorderUnit.form_pairs_from_asym(psites)
    symmops = get_symmop(pmg_dict)
    fin_psites = apply_symmop(psites, symmops)
    return fin_psites, disunit_pairs, len(symmops)


def get_symmop(data):
    symops = []
    for symmetry_label in ["_symmetry_equiv_pos_as_xyz",
                           "_symmetry_equiv_pos_as_xyz_",
                           "_space_group_symop_operation_xyz",
                           "_space_group_symop_operation_xyz_"]:
        if data.get(symmetry_label):
            xyz = data.get(symmetry_label)
            if isinstance(xyz, str):
                msg = "A 1-line symmetry op P1 CIF is detected!"
                warnings.warn(msg)
                xyz = [xyz]
            try:
                symops = [SymmOp.from_xyz_string(s)
                          for s in xyz]
                break
            except ValueError:
                continue
    if not symops:
        # Try to parse symbol
        for symmetry_label in ["_symmetry_space_group_name_H-M",
                               "_symmetry_space_group_name_H_M",
                               "_symmetry_space_group_name_H-M_",
                               "_symmetry_space_group_name_H_M_",
                               "_space_group_name_Hall",
                               "_space_group_name_Hall_",
                               "_space_group_name_H-M_alt",
                               "_space_group_name_H-M_alt_",
                               "_symmetry_space_group_name_hall",
                               "_symmetry_space_group_name_hall_",
                               "_symmetry_space_group_name_h-m",
                               "_symmetry_space_group_name_h-m_"]:
            sg = data.get(symmetry_label)

            if sg:
                sg = sub_spgrp(sg)
                try:
                    spg = space_groups.get(sg)
                    if spg:
                        symops = SpaceGroup(spg).symmetry_ops
                        msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
                              "Spacegroup from %s used." % symmetry_label
                        warnings.warn(msg)
                        break
                except ValueError:
                    # Ignore any errors
                    pass

                try:
                    for d in _get_cod_data():
                        if sg == re.sub(r"\s+", "",
                                        d["hermann_mauguin"]):
                            xyz = d["symops"]
                            symops = [SymmOp.from_xyz_string(s)
                                      for s in xyz]
                            msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
                                  "Spacegroup from %s used." % symmetry_label
                            warnings.warn(msg)
                            break
                except Exception:
                    continue

                if symops:
                    break
    if not symops:
        # Try to parse International number
        for symmetry_label in ["_space_group_IT_number",
                               "_space_group_IT_number_",
                               "_symmetry_Int_Tables_number",
                               "_symmetry_Int_Tables_number_"]:
            if data.get(symmetry_label):
                try:
                    i = int(str2float(data.get(symmetry_label)))
                    symops = SpaceGroup.from_int_number(i).symmetry_ops
                    break
                except ValueError:
                    continue

    if not symops:
        msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
              "Defaulting to P1."
        warnings.warn(msg)
        symops = [SymmOp.from_xyz_string(s) for s in ['x', 'y', 'z']]

    return symops


class DisorderUnit:
    def __repr__(self):
        return "\n".join(self.labels)

    def __init__(self, labels, occu, disorder_appendix="'|\?"):
        """
        #TODO what if one portion is occupied by >2 units?
        a collection of psites representing one possiblity for a portion of an asym unit
        a pair of xor DisorderUnit objs on an asym unit is the basis of disorder, assuming maximal entropy
        xor pair is defined by the label, e.g. [A1, B2, C3] is paired with [A1?, B2?, C3?]

        """
        self.disorder_appendix = disorder_appendix
        self.occu = occu
        self.labels = labels
        self.label_numbers = [re.search(r"\d+", lab).group(0) for lab in self.labels]
        self.number_set = set(self.label_numbers)

    @staticmethod
    def form_pairs_from_asym(psites):
        """
        get a list of xor disunit pairs

        :param psites:
        :return:
        """
        units = DisorderUnit.labels_from_aym_sites(psites)
        pairs = []
        assigned = []
        for i in range(len(units)):
            if i not in assigned:
                u1 = units[i]
                for j in range(i + 1, len(units)):
                    u2 = units[j]
                    if abs(u1.occu + u2.occu - 1) < 1e-5 and u1.number_set == u2.number_set:
                        pairs.append([u1, u2])
                        assigned += [i, j]
        return pairs

    @staticmethod
    def labels_from_aym_sites(psites):
        """
        all disorder units in a flat list within in one asym unit

        :param psites:
        :return:
        """
        disordered_sites = [ps for ps in psites if abs(ps.properties['occu'] - 1) > 1e-3]
        # first group by occu
        group_by_occu = [list(v) for l, v in groupby(sorted(disordered_sites, key=lambda x: x.properties['occu']),
                                                     lambda x: x.properties['occu'])]
        # within one group, group again by connection
        units = []
        for i in range(len(group_by_occu)):
            group = group_by_occu[i]
            mols, unwrap_str_sorted, unwrap_pblock_list = PBCparser.unwrap(Structure.from_sites(group))
            for pblock in unwrap_pblock_list:
                units.append(DisorderUnit([s.properties['label'] for s in pblock], pblock[0].properties['occu']))
        return units


def build_supercell_full_disorder(pstructure, scaling_matrix):
    """
    get a supercell with all disordered sites inside, should be used to generate a certain config based on instruction

    :param pstructure:
    :param scaling_matrix:
    :return:
    """
    scale_matrix = np.array(scaling_matrix, np.int16)
    if scale_matrix.shape != (3, 3):
        scale_matrix = np.array(scale_matrix * np.eye(3), np.int16)
    new_lattice = Lattice(np.dot(scale_matrix, pstructure._lattice.matrix))
    f_lat = lattice_points_in_supercell(scale_matrix)
    c_lat = new_lattice.get_cartesian_coords(f_lat)
    new_sites = []
    for site in pstructure.sites:
        icell = 0
        for v in c_lat:
            site_properties = deepcopy(site.properties)
            site_properties['icell'] = icell
            s = PeriodicSite(site.species, site.coords + v,
                             new_lattice, properties=site_properties,
                             coords_are_cartesian=True, to_unit_cell=False)
            new_sites.append(s)
            icell += 1
    new_charge = pstructure._charge * np.linalg.det(scale_matrix) if pstructure._charge else None
    return Structure.from_sites(new_sites, charge=new_charge), len(c_lat)


def dissc_to_config(sc, disunit_pairs, instructions):
    """
    take instructions to generate a certain config from super cell

    :param sc:
    :param disunit_pairs:
    :param instructions:
    :return:
    """
    pool = []
    config_occu = 1
    for icell in range(len(instructions)):
        for iasym in range(len(instructions[icell])):
            for idis in range(len(instructions[icell][iasym])):
                disu = disunit_pairs[idis][instructions[icell][iasym][idis]]
                config_occu *= disu.occu
                for label in disu.labels:
                    pool.append((label, iasym, icell))
    conf_sites = []
    for ps in sc.sites:
        if abs(ps.properties['occu'] - 1) < 1e-3:
            conf_sites.append(ps)
        else:
            plabel = ps.properties['label']
            piasym = ps.properties['iasym']
            picell = ps.properties['icell']
            if (plabel, piasym, picell) in pool:
                conf_sites.append(ps)
    return Structure.from_sites(conf_sites), config_occu


def gen_instructions(disunitpairs, n_asym, n_cell):
    """
    get all possible instructions, notice this is expentially scaled

    # of configs = 2^len(disunitpairs)^n_asym^n_cell where 2 comes from pairwise disordered occupancies

    :param disunitpairs:
    :param n_asym:
    :param n_cell:
    :return:
    """
    # icell, iasym, ipair, idis

    pair_size = []
    for ipair in range(len(disunitpairs)):
        pair_size.append(len(disunitpairs[ipair]))
    dis_ins_combinations = list(itertools.product(*[list(range(n)) for n in pair_size]))
    results = list(itertools.product(dis_ins_combinations, repeat=n_asym))
    results = list(itertools.product(results, repeat=n_cell))
    return results



