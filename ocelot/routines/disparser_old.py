import itertools
import math
import re
import warnings
from collections import OrderedDict
from copy import deepcopy
from itertools import groupby

import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure import Lattice
from pymatgen.core.structure import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure import lattice_points_in_supercell
from pymatgen.io.cif import CifFile
from pymatgen.io.cif import _get_cod_data
from pymatgen.io.cif import str2float
from pymatgen.io.cif import sub_spgrp
from pymatgen.symmetry.groups import SYMM_DATA
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.util.coord import pbc_shortest_vectors

from ocelot.routines.pbc import PBCparser


class DisorderParserError(Exception):
    pass


"""
DisParser: parse cif file into a list of configurations with no disorder

this is used to parse disordered cif file to configurations of a supercell

plz read 
[CSD API convention](https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/crystal.html#disorder) 
first 

we assume the cif file contains following disorder-related fields
    _atom_site_occupancy 
    _atom_site_disorder_group 
if both of the two keys exist, then we can extract configuration from the cif file with explicit occupancies

if at least one of above fields does not exist, it is still possible the cif file contains disorder, this depends on
the the values of dict[_atom_site_label], btw, '_atom_site_label' may not exist when using jmol to write cif

the values of dict[_atom_site_label] form a list of *unique* labels, a label can be considered as a concatenation of 
    
    E (an element symbol) + I (an index) + T (the third tag)
    
while CSD claims the alternative configurations should be labelled with a suffix question mark (?),
this is certainly not true for e.g. ASIXEH

so here a summary of representing disorder in cif file

    1. EI is unique for each site, sites of the alternative configuration are labeled by a non-word T
        in this case you cannot get >1 alternative configurations, e.g. ALOVOO
        
    2. EI is not unique for each site, site with a non-word T is an alternative config for the site having the same EI, 
        e.g. x17059
    
    3. EI is not unique for each site, there is no site with a non-word T. Sites sharing the same EI representing 
        disorder, e.g. ASIXEH
    
    4. disordered sites are simply ignored in the cif file s.t. the molecule is not complete, e.g. ABEGET01
    
    5. disordered sites are simply ignored in the cif file, but the molecule looks legit, e.g. AGAJUO
    
    6. disordered sites are present in the cif file but not labeled at all, e.g. ANOPEA

in case 1. it is possible to get the occupancies via CSD API, not always tho.

I cannot find a way to process 3., 4., 5., 6. without checking deposition records or chemistry analysis, so this parser
will ignore them. That is, there can be only one alternative configuration.

I will deal with 1. by artificially set (average) _atom_site_occupancy and _atom_site_disorder_group based on special character 
suffix, this should be done intentionally

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

space_groups = {sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()}
space_groups.update({sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()})


def get_pmg_dict(cifstring: str):
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

    # jmol writes '_atom_site_type_symbol', but not '_atom_site_label'
    if '_atom_site_label' not in pymatgen_dict.keys():
        warnings.warn('W: _atom_site_label not found in parsed dict')
        atom_site_label = []
        symbols = pymatgen_dict['_atom_site_type_symbol']
        for i in range(len(symbols)):
            s = symbols[i]
            atom_site_label.append('{}{}'.format(s, i))
        pymatgen_dict['_atom_site_label'] = atom_site_label
    return identifier, pymatgen_dict


class CifFileError(Exception):
    pass


class AtomLabel:
    def __init__(self, label: str):
        """
        a class for atom label in cif file, overkill I guess...

        :param label:
        """
        self.label = label

        tmplabel = self.label
        self.element = re.findall(r"^[a-zA-Z]+", tmplabel)[0]

        tmplabel = tmplabel.lstrip(self.element)
        self.index = re.findall(r"^\d+", tmplabel)[0]

        tmplabel = tmplabel.lstrip(str(self.index))
        self.tag = tmplabel

        self.index = int(self.index)

        self.ei = "{}{}".format(self.element, self.index)

    def __str__(self):
        return self.label

    def __repr__(self):
        return self.label

    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return self.label == other.label


class CifLegitError(Exception):
    pass


def is_nonword(s: str):
    return re.search(r"^\W$", s)


def get_labels_by_ei(al: AtomLabel, als):
    """
    get a list of atomlabel whose ei == al.ei

    :param al:
    :param als:
    :return:
    """
    sameei = []
    for alj in als:
        if al.ei == alj.ei and al != alj:
            sameei.append(alj)
    return sameei


class DisParser:  # chaos parser sounds cooler?

    def __init__(self, cifstring: str):
        self.cifstring = cifstring
        self.identifier, self.data = get_pmg_dict(self.cifstring)
        self.data['_atom_site_fract_x'] = [braket2float(x) for x in self.data['_atom_site_fract_x']]
        self.data['_atom_site_fract_y'] = [braket2float(x) for x in self.data['_atom_site_fract_y']]
        self.data['_atom_site_fract_z'] = [braket2float(x) for x in self.data['_atom_site_fract_z']]

        try:
            labels = self.data['_atom_site_label']
            self.labels = [AtomLabel(lab) for lab in labels]
            self.eis = [al.ei for al in self.labels]
            self.tags = [al.tag for al in self.labels]
            if len(self.labels) != len(set(self.labels)):
                raise CifFileError('duplicate labels found in the cifstring!')
        except KeyError:
            self.labels = None
            self.eis = None
            self.tags = None
            raise CifFileError('no _atom_site_type_symbol field in the cifstring!')
        self.latparams = [self.data[k] for k in latt_labels]
        self.lattice = Lattice.from_parameters(*[braket2float(p) for p in self.latparams], True)

    @property
    def labels_with_nonword_suffix(self):
        return [l for l in self.labels if is_nonword(l.tag)]

    @property
    def labels_with_word_or_no_suffix(self):
        return [l for l in self.labels if not is_nonword(l.tag)]

    @classmethod
    def from_ciffile(cls, fn):
        with open(fn, 'r') as f:
            s = f.read()
        return cls(s)

    def classify(self):
        """
        one cif file belongs to one of the following categories:

            dis-0: has _atom_site_occupancy AND _atom_site_disorder_group

            nodis-0: no dup in self.eis, set(self.tags) is {""}

            dis-1:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- EI

            dis-2:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- E'I' e.g. ALOVOO.cif

            nodis-1: dup in self.eis, set(self.tags) is ["", <word>, ...], this could be a dis as in ASIXEH

            weird: else
        """
        if '_atom_site_occupancy' in self.data.keys() and '_atom_site_disorder_group' in self.data.keys():
            return 'dis-0'
        else:
            tag_set = set(self.tags)
            if tag_set == {""}:
                return 'nodis-0'
            else:
                len_nonword_tagset = len([tag for tag in tag_set if is_nonword(tag)])
                if len_nonword_tagset == 1:
                    if all(len(get_labels_by_ei(al, self.labels)) == 0 for al in self.labels_with_nonword_suffix):
                        return 'dis-2'
                    elif all(len(get_labels_by_ei(al, self.labels)) == 1 for al in self.labels_with_nonword_suffix):
                        return 'dis-1'
                    else:
                        return 'weird'
                elif len_nonword_tagset > 1:
                    raise DisorderParserError('more than one possible nonword suffix')
                else:
                    return 'nodis-1'

    def get_coord_data(self):
        """
        coord_data[<atom_label>] = [x, y, z, symbol]
        """
        coord_data = map(list, zip(self.data['_atom_site_fract_x'],
                                   self.data['_atom_site_fract_y'],
                                   self.data['_atom_site_fract_z'],
                                   self.data['_atom_site_type_symbol']
                                   ))
        coord_data = list(coord_data)
        coord_data = OrderedDict(zip(self.labels, coord_data))
        return coord_data

    @staticmethod
    def coorddata2newdata(coord_data, cifdata):
        """
        coord_data[atomlable] = x, y, z, symbol, occu, group
        """
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
        newdata['_atom_site_label'] = [l.label for l in labs]
        newdata['_atom_site_type_symbol'] = symbols
        newdata['_atom_site_fract_x'] = xs
        newdata['_atom_site_fract_y'] = ys
        newdata['_atom_site_fract_z'] = zs
        newdata['_atom_site_occupancy'] = occus
        newdata['_atom_site_disorder_group'] = idisgs
        return newdata

    def dis2_to_dis0(self, cutoff=1.5):
        """
        c20 <--> c16?

        dis-1: no dup in self.eis, set(self.tags) is ["", <non-word>, <word>, ...]
        dis-0: has _atom_site_occupancy AND _atom_site_disorder_group
        e.g. ALOVOO.cif

        this is rather unreliable, user should be warned
        for each label with <non-word> tag, find nearest label (within the cutoff) as a counterpart
        exception includes AGUHUG.cif

        """
        warnings.warn('W: trying to convert dis-2 to dis-0 with cutoff {}, this is unreliable!'.format(cutoff))
        ali2alj = OrderedDict()
        coord_data = self.get_coord_data()
        for ali in self.labels_with_nonword_suffix:
            ali_nbs_dictance = []
            xi, yi, zi, symboli = coord_data[ali]
            fci = [xi, yi, zi]
            for alj in self.labels_with_word_or_no_suffix:
                xj, yj, zj, symbolj = coord_data[alj]
                fcj = [xj, yj, zj]
                if symboli != symbolj:
                    continue
                v, d2 = pbc_shortest_vectors(self.lattice, fci, fcj, return_d2=True)
                dij = math.sqrt(d2[0, 0])
                ali_nbs_dictance.append([alj, dij])
            if len(ali_nbs_dictance) == 0:
                raise DisorderParserError(
                    'dis1 to dis0 failed as cannot find any same-symbol neighbors of label {}'.format(ali.label))
            ali_nbs_dictance = sorted(ali_nbs_dictance, key=lambda x: x[1])
            nnb, nnbdis = ali_nbs_dictance[0]
            if nnbdis > cutoff:
                raise DisorderParserError(
                    'dis1 to dis0 failed as cannot find any same-symbol neighbors of label {} within {}'.format(
                        ali.label, cutoff))
            else:
                ali2alj[ali] = nnb

        if len(ali2alj.keys()) == len(self.labels_with_nonword_suffix):
            for ali in self.labels:
                if ali in ali2alj.keys():
                    coord_data[ali].append(0.5)
                    coord_data[ali].append('2')
                    coord_data[ali2alj[ali]].append(0.5)
                    coord_data[ali2alj[ali]].append('1')
                elif ali not in ali2alj.values():
                    coord_data[ali].append(1)
                    coord_data[ali].append('.')
            return self.coorddata2newdata(coord_data, self.data)

        else:
            raise DisorderParserError('not all disordered label have a counterpart')

    # def dis2_to_dis0(self):
    #     """
    #     c20 <--> c20?
    #
    #     dis-0: has _atom_site_occupancy AND _atom_site_disorder_group
    #
    #     dis-2: dup in self.eis, set(self.tags) is ["", <non-word>]
    #     """
    #     coord_data = self.get_coord_data()
    #     group_by_ei = [list(v) for l, v in
    #                      groupby(sorted(coord_data.keys(), key=lambda x: x.ei,),
    #                              key=lambda x: x.ei)]
    #     for i in range(len(group_by_ei)):
    #         labels = sorted(group_by_ei[i])  # this works if all disorder charaters are added as suffix
    #         if len(labels) == 1:
    #             coord_data[labels[0]].append('1')
    #             coord_data[labels[0]].append('.')  # convention for no disorder
    #         else:
    #             occu = 1.0 / len(labels)
    #             idisg = 1  # convention start from 1
    #             for lab in labels:
    #                 coord_data[lab].append(str(occu))
    #                 coord_data[lab].append(str(idisg))
    #                 idisg += 1

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
            labels = sorted(group_by_word[i])  # this works if all disorder charaters are added as suffix
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
    def dis0_to_configs(cifdata, write_files=False, scaling_mat=(1, 1, 1)):
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
        if classification == 'dis-0':
            cifdata = self.data
        elif classification == 'dis-2':
            cifdata = self.dis2_to_dis0()
        elif classification == '3':
            cifdata = self.class3_to_class1(self.data)
        else:
            raise NotImplementedError('to_configs is not implemented for classification {}'.format(classification))
        return self.dis0_to_configs(cifdata, write_files, scaling_mat)


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


def braket2float(s):
    try:
        return float(s)
    except ValueError:
        if isinstance(s, str):
            return str2float(s)
        raise TypeError('cannot parse {} into float'.format(s))


def get_psites(pmg_dict):
    """
    get pmg_psites from cif data of dis-0

    :param pmg_dict:
    :return: fin_psites, a list of psites with 'occu', 'disg', 'label' in properties
            disunit_pairs, a list of disorderunit pairs in one aym unit, e.g. [[u1, u1'], [u2, u2']] where u1 xor u1'
            len(symmops), number of aym units in a cell
    """
    labels = pmg_dict['_atom_site_label']
    symbols = pmg_dict['_atom_site_type_symbol']
    xs = [braket2float(j) for j in pmg_dict['_atom_site_fract_x']]
    ys = [braket2float(j) for j in pmg_dict['_atom_site_fract_y']]
    zs = [braket2float(j) for j in pmg_dict['_atom_site_fract_z']]
    occus = [braket2float(j) for j in pmg_dict['_atom_site_occupancy']]
    disgs = pmg_dict['_atom_site_disorder_group']
    psites = []
    length_strings = ("a", "b", "c")
    angle_strings = ("alpha", "beta", "gamma")
    lengths = [braket2float(pmg_dict["_cell_length_" + i]) for i in length_strings]
    angles = [braket2float(pmg_dict["_cell_angle_" + i]) for i in angle_strings]
    lattice = Lattice.from_parameters(*lengths, *angles, vesta=True)
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
                    i = int(braket2float(data.get(symmetry_label)))
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
        return "\n".join(str(self.labels))

    def __init__(self, psites):
        """
        #TODO what if one portion is occupied by >2 units?
        a collection of psites representing one possiblity for a portion of an asym unit

        a pair of xor DisorderUnit objs on an asym unit is the basis of disorder, assuming maximal entropy

        xor pair is defined by the label, e.g. [A1, B2, C3] is paired with [A1?, B2?, C3?]

        """
        self.psites = psites
        self.occu = [s.properties['occu'] for s in psites]
        self.labels= [s.properties['label'] for s in psites]

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
    get all possible instructions, notice this is exponentially scaled

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
