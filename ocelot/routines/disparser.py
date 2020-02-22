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
        still you cannot get >1 alternative configurations e.g. x17059
    
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

        if len(self.tag) > 1:
            raise AtomLabelError('tag for {} is {}, this is unlikely'.format(label, self.tag))

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

    @property
    def is_tag_nonword(self):
        return re.search(r"^\W$", self.tag)

    def get_labels_with_same_ei(self, als):
        """
        get a list of atomlabel whose ei == al.ei

        :param als:
        :return:
        """
        sameei = []
        for alj in als:
            if self.ei == alj.ei and self != alj:
                sameei.append(alj)
        return sameei


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


def apply_symmop(psites, ops):
    """
    symmop and xyz in cif file:

    lets say xyz -- op1 --> x'y'z' and xyz -- op2 --> x!y!z! and
    it is possible to have x'y'z' is_close x!y!z!

    this means one should take only x'y'z' or x!y!z!, aka op1 is equivalent to op2 due to the symmetry implicated by
    xyz/the asymmectric unit, e.g. ALOVOO.cif -- Z=2, asymmectric unit given by cif is one molecule, but there're 4 ops

    so we need first check if the cif file behaves like this

    """
    op_xyzs = []
    for op in ops:
        n_xyzs = []
        for ps in psites:
            new_coord = op.operate(ps.frac_coords)
            # new_coord = np.array([i - math.floor(i) for i in new_coord])
            n_xyzs.append(new_coord)
        op_xyzs.append(n_xyzs)

    latt = psites[0].lattice

    def pbc_dist(fc1, fc2, lattice):
        v, d2 = pbc_shortest_vectors(lattice, fc1, fc2, return_d2=True)
        return math.sqrt(d2[0, 0])

    def pbc_distmat(fcl1, fcl2):
        distmat = np.zeros((len(fcl1), len(fcl1)))
        for i in range(len(fcl1)):
            for j in range(i, len(fcl1)):
                distmat[i][j] = pbc_dist(fcl1[i], fcl2[j], latt)
                distmat[j][i] = distmat[i][j]
        return distmat

    def two_xyzs_close(xyzs1, xyzs2, tol=1e-5):
        dmat = pbc_distmat(xyzs1, xyzs2)
        almost_zeros = dmat[(dmat < tol)]
        if len(almost_zeros) > 0:
            return True
        return False

    op_identities = np.zeros((len(ops), len(ops)), dtype=bool)
    for i, j in itertools.combinations(range(len(ops)), 2):
        ixyzs = op_xyzs[i]
        jxyzs = op_xyzs[j]
        if two_xyzs_close(ixyzs, jxyzs):
            op_identities[i][j] = True
            op_identities[j][i] = True

    groups = [[0]]
    for i in range(len(ops)):
        for ig in range(len(groups)):
            if all(op_identities[i][j] for j in groups[ig]):
                groups[ig].append(i)
        if i not in [item for sublist in groups for item in sublist]:
            groups.append([i])
    unique_ops = [ops[g[0]] for g in groups]

    new_psites = []
    for ps in psites:
        iasym = 0
        for op in unique_ops:
            new_coord = op.operate(ps.frac_coords)
            new_coord = np.array([i - math.floor(i) for i in new_coord])
            new_properties = deepcopy(ps.properties)
            new_properties['iasym'] = iasym
            new_ps = PeriodicSite(ps.species_string, new_coord, ps.lattice, properties=deepcopy(new_properties))
            new_psites.append(new_ps)
            iasym += 1

    return new_psites, unique_ops


def braket2float(s):
    try:
        return float(s)
    except ValueError:
        if isinstance(s, str):
            return str2float(s)
        raise TypeError('cannot parse {} into float'.format(s))


def get_psite_by_atomlable(psites, al):
    for s in psites:
        if s.properties['label'] == str(al):
            return s
    raise AtomLabelError('cannot find psite with atomlable {}'.format(str(al)))


def get_psite_label(s):
    return AtomLabel(s.properties['label'])


def get_nearest_label(data: dict, ali: AtomLabel, neighbor_labels: [AtomLabel], lattice: Lattice, cutoff=1.5):
    """
    given a list of potential neighbouring AtomLable, get the one that is closest and has the same element

    :param data:
    :param ali:
    :param neighbor_labels:
    :param lattice:
    :param cutoff:
    :return:
    """
    ali_nbs_dictance = []
    xi, yi, zi, symboli, occu, disgrp = data[ali]
    fci = [xi, yi, zi]
    for alj in neighbor_labels:
        if alj == ali:
            continue
        xj, yj, zj, symbolj = data[alj][:4]
        fcj = [xj, yj, zj]
        if symboli != symbolj:
            continue
        v, d2 = pbc_shortest_vectors(lattice, fci, fcj, return_d2=True)
        dij = math.sqrt(d2[0, 0])
        ali_nbs_dictance.append([alj, dij])
    if len(ali_nbs_dictance) == 0:
        raise DisorderParserError(
            'cannot find any same-symbol neighbors of label {}'.format(ali.label))
    ali_nbs_dictance = sorted(ali_nbs_dictance, key=lambda x: x[1])
    nnb, nnbdis = ali_nbs_dictance[0]
    if nnbdis > cutoff:
        raise DisorderParserError(
            'cannot find any same-symbol neighbors of label {} within {}'.format(
                ali.label, cutoff))
    return nnb


class CifFileError(Exception): pass


class DisorderParserError(Exception): pass


class AtomLabelError(Exception): pass


class DisParser:  # chaos parser sounds cooler?

    labels: [AtomLabel]

    def __init__(self, cifstring: str):
        """
        this can only handle one alternative configuration for the asymmetric unit

        one asymmetric unit == inv_conf + disg1 + disg2

        disg1 = disunit_a + disunit_b + ...

        disg2 = disunit_a' + disunit_b' + ...

        DisorderPair_a = disunit_a + disunit_a'

        if the cif file contains previouly fitted occu and disg,
        we call it dis-0 and we use occu/disg info to get inv_conf, disg1, disg2

        if there is no previously fitted info, we deal with the following situations:

            dis-1:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- EI,

            note x17059.cif is dis-1 but it has hydrogens like H12A -- H12D, this can only be captured
            by previously fitted disg and occu,

            dis-2:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- E'I' e.g. ALOVOO.cif

            nodis-0: no dup in self.eis, set(self.tags) is {""}

            nodis-1: dup in self.eis, set(self.tags) is ["", <word>, ...], this could be a dis as in ASIXEH

            weird: else

        for dis-1, dis-2, we fill the occu, disg fields in cifdata, so they can be coonverted to dis-0

        Attributes:
            data[atomlabel] = [x, y, z, symbol, occu, disgrp]
        """
        self.cifstring = cifstring
        self.identifier, self.cifdata = get_pmg_dict(self.cifstring)

        if '_atom_site_occupancy' in self.cifdata.keys() and '_atom_site_disorder_group' in self.cifdata.keys():
            self.was_fitted = True
        else:
            self.was_fitted = False

        if self.was_fitted:
            self.cifdata['_atom_site_occupancy'] = [braket2float(x) for x in self.cifdata['_atom_site_occupancy']]
            self.cifdata['_atom_site_disorder_group'] = [braket2float(x) for x in
                                                         self.cifdata['_atom_site_disorder_group']]

        self.cifdata['_atom_site_fract_x'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_x']]
        self.cifdata['_atom_site_fract_y'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_y']]
        self.cifdata['_atom_site_fract_z'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_z']]
        try:
            labels = self.cifdata['_atom_site_label']
        except KeyError:
            raise CifFileError('no _atom_site_label field in the cifstring!')
        self.labels = [AtomLabel(lab) for lab in labels]
        self.eis = [al.ei for al in self.labels]
        self.tags = [al.tag for al in self.labels]
        if len(self.labels) != len(set(self.labels)):
            raise CifFileError('duplicate labels found in the cifstring!')
        self.latparams = [self.cifdata[k] for k in latt_labels]
        self.lattice = Lattice.from_parameters(*[braket2float(p) for p in self.latparams], True)

        data = map(list, zip(
            self.cifdata['_atom_site_fract_x'],
            self.cifdata['_atom_site_fract_y'],
            self.cifdata['_atom_site_fract_z'],
            self.cifdata['_atom_site_type_symbol'],
        ))
        data = list(data)
        data = OrderedDict(zip(self.labels, data))
        self.data = data

        for i in range(len(self.labels)):
            al = self.labels[i]
            if self.was_fitted:
                self.data[al].append(self.cifdata['_atom_site_occupancy'][i])
                self.data[al].append(self.cifdata['_atom_site_disorder_group'][i])
            else:
                self.data[al].append(None)
                self.data[al].append(None)

    def dis0data_to_disunit_pairs(self):
        psites = self.get_psites_from_data()
        disu_pairs, inv_conf = DisUnit.get_disunit_pairs_from_asym(psites)
        return disu_pairs, inv_conf

    def get_psites_from_data(self):
        """
        get psites from self.data, each psite will be assigned properties with fields
        occu
        disg
        label

        :return:
        """
        ps = []
        for k in self.data.keys():
            x, y, z, symbol, occu, disgrp = self.data[k]
            if occu is None or disgrp is None:
                raise DisorderParserError('getting psites with None occu/disgrp')
            ps.append(PeriodicSite(symbol, [x, y, z], properties={'occu': occu, 'disg': disgrp, 'label': k.label},
                                   lattice=self.lattice))
        return ps

    @property
    def labels_with_nonword_suffix(self):
        return [l for l in self.labels if l.is_tag_nonword]

    @property
    def labels_with_word_or_no_suffix(self):
        return [l for l in self.labels if not l.is_tag_nonword]

    @classmethod
    def from_ciffile(cls, fn):
        with open(fn, 'r') as f:
            s = f.read()
        return cls(s)

    def classify(self):
        """
        one cif file belongs to one of the following categories:

            nodis-0: no dup in self.eis, set(self.tags) is {""}

            dis-1:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- EI,

            note x17059.cif has hydrogens like H12A -- H12D, this can only be captured
            by previously fitted disg and occu, in general there
            should be a check on whether disgs identified by the parser are identical to previously fitted

            dis-2:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- E'I' e.g. ALOVOO.cif

            nodis-1: dup in self.eis, set(self.tags) is ["", <word>, ...], this could be a dis as in ASIXEH

            weird: else
        """
        tag_set = set(self.tags)
        if self.was_fitted:
            return "dis-0"

        elif tag_set == {""}:
            return "nodis-0"
        else:
            len_nonword_tagset = len(set([l.tag for l in self.labels if l.is_tag_nonword]))
            if len_nonword_tagset == 1:
                if all(len(al.get_labels_with_same_ei(self.labels)) == 0 for al in self.labels_with_nonword_suffix):
                    return 'dis-2'
                if all(len(al.get_labels_with_same_ei(self.labels)) == 1 for al in self.labels_with_nonword_suffix):
                    return 'dis-1'
                else:
                    return 'weird'
            elif len_nonword_tagset > 1:
                raise DisorderParserError('more than one possible nonword suffix')
            else:
                return 'nodis-1'

    @staticmethod
    def data2cifdata(data, cifdata):
        """
        data[atomlable] = x, y, z, symbol, occu, group
        """
        newdata = OrderedDict()
        for ciflabel in possible_symm_labels + latt_labels + chemistry_labels:
            if ciflabel in cifdata.keys():
                newdata[ciflabel] = cifdata[ciflabel]
        labs = list(data.keys())
        xs = []
        ys = []
        zs = []
        symbols = []
        occus = []
        idisgs = []
        for lab in labs:
            x, y, z, symb, occu, idisg = data[lab]
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

    def parse(self):
        classification = self.classify()
        print('{} thinks this cif file belongs to class {}'.format(self.__class__.__name__, classification))
        if classification in ['nodis-1', 'nodis-0']:
            self.nodis_to_dis0()
        elif classification == 'weird':
            raise CifFileError('classified as weird, check your cif input please')

        elif classification == 'dis-0':
            pass

        elif classification == 'dis-1':
            self.dis1_to_dis0()

        elif classification == 'dis-2':
            self.dis2_to_dis0()
        else:
            raise DisorderParserError('unknown classification!')
        disunit_pairs, inv_conf = self.dis0data_to_disunit_pairs()
        return disunit_pairs, inv_conf

    def nodis_to_dis0(self):
        for al in self.labels:
            self.data[al][4] = 1
            self.data[al][5] = '.'

    def alijdict_to_disgs_and_update_data(self, ali2alj: dict):
        inv_conf = []
        disg_a = []  # alj, ends with word/empty suffix
        disg_b = []  # ali, ends with non-word suffix
        for al in self.labels:
            if al in ali2alj.keys():
                disg_b.append(al)
                self.data[al][4] = 0.49
                self.data[al][5] = '2'
            elif al in ali2alj.values():
                disg_a.append(al)
                self.data[al][4] = 0.51
                self.data[al][5] = '1'
            else:
                inv_conf.append(al)
                if not self.was_fitted:
                    self.data[al][4] = 1
                    self.data[al][5] = '.'
        al_b2a = ali2alj
        al_a2b = {v: k for k, v in ali2alj.items()}

        return al_a2b, al_b2a, disg_a, disg_b, inv_conf

    def dis1_to_dis0(self):
        """
        C20 -- C20?

        dis-1:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- EI,
        """
        ali2alj = OrderedDict()
        for ali in self.labels_with_nonword_suffix:
            potential_matches = self.labels_with_word_or_no_suffix
            possible_alj = ali.get_labels_with_same_ei(potential_matches)
            if len(possible_alj) != 1:
                raise DisorderParserError('possible match for {} is not 1'.format(ali))
            alj = possible_alj[0]
            ali2alj[ali] = alj

        return self.alijdict_to_disgs_and_update_data(ali2alj)

    def dis2_to_dis0(self, cutoff=1.5):
        """
        c20 <--> c16?

        dis-2:  set(self.tags) is ["", <non-word>, <word>, ...], EI<non-word> -- E'I' e.g. ALOVOO.cif

        will write default occu and disg to self.data

        this is rather unreliable, user should be warned
        for each label with <non-word> tag, find nearest label (within the cutoff) as a counterpart
        exception includes AGUHUG.cif (H -- OH disorder)

        """
        warnings.warn('W: trying to parse disorder in dis-2 with cutoff {}, this is unreliable!'.format(cutoff))
        ali2alj = OrderedDict()
        for ali in self.labels_with_nonword_suffix:
            potential_matches = self.labels_with_word_or_no_suffix
            alj = get_nearest_label(self.data, ali, [l for l in potential_matches if
                                                     l not in ali2alj.keys() and l not in ali2alj.values()],
                                    self.lattice, cutoff)
            ali2alj[ali] = alj
        return self.alijdict_to_disgs_and_update_data(ali2alj)

    def to_configs(self, write_files=False, scaling_mat=(1, 1, 1), assign_siteid=True):
        disunit_pairs, inv_conf = self.parse()
        cc = ConfigConstructor(disunit_pairs, inv_conf)

        psites = self.get_psites_from_data()

        raw_symmops = get_symmop(self.cifdata)
        psites, symmops = apply_symmop(psites, raw_symmops)

        if assign_siteid:
            print('siteid is assigned in dp.to_configs')
            for isite in range(len(psites)):
                psites[isite].properties['siteid'] = isite

        pstructure = Structure.from_sites(psites, to_unit_cell=True)

        # sc, n_unitcell = ConfigConstructor.build_supercell_full_disorder(pstructure, scaling_mat)
        mols, unwrap_str_sorted, unwrap_pblock_list = PBCparser.unwrap(pstructure)
        sc, n_unitcell = ConfigConstructor.build_supercell_full_disorder(unwrap_str_sorted, scaling_mat)
        conf_ins = cc.gen_instructions(disunit_pairs, len(symmops), n_unitcell)
        iconf = 0
        confs = []
        for confin in conf_ins:
            conf, conf_occu = ConfigConstructor.dissc_to_config(sc, disunit_pairs, confin)
            confs.append([conf, conf_occu, mols])  # this molecule still contains disordered sites!
            if write_files:
                conf.to('cif', 'conf_{}.cif'.format(iconf))  # pymatgen somehow does not write disg field in the cif
            iconf += 1
        if write_files:
            pstructure.to('cif', 'confgen_ps.cif')
            # unwrap_str_sorted.to('cif', 'confgen_unwrap.cif')
        return sorted(confs, key=lambda x: x[1], reverse=True)


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


class DisUnitError(Exception): pass


class DisUnit:
    def __repr__(self):
        return " ".join([str(l) for l in self.labels])

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return set(self.labels) == set(other.labels)

    def __hash__(self):
        return hash(set(self.labels))

    def __init__(self, psites: [PeriodicSite]):
        """
        a portion of an asymmetric unit representing one possibility

        a pair of DisUnit is the basis of disorder, assuming maximal entropy
        """
        self.sites = psites
        self.labels = [AtomLabel(s.properties['label']) for s in self.sites]
        occus = [s.properties['occu'] for s in self.sites]
        disgs = [s.properties['disg'] for s in self.sites]
        self.occu = occus[0]
        self.disg = disgs[0]
        self.symbols = [l.element for l in self.labels]
        self.composition = tuple(sorted(self.symbols))
        if len(set(occus)) != 1:
            raise DisUnitError('occu is not uniform for {}'.format(str(self)))
        if len(set(disgs)) != 1:
            raise DisUnitError('disg is not uniform for {}'.format(str(self)))

    @property
    def geoc(self):
        c = np.zeros(3)
        for s in self.sites:
            c += s.coords
        return c / len(self.sites)

    @staticmethod
    def get_disunit_pairs_from_asym(psites):
        """
        get a list of xor disunit pairs

        :param psites:
        :return: [[ua, ua'], [ub, ub'], ...]
        """
        units, inv_conf = DisUnit.partition_asymmetric_unit(psites)
        pairs = []
        assigned = []
        for i in range(len(units)):
            if i not in assigned:
                u1 = units[i]
                potential_u2 = []
                for j in range(i + 1, len(units)):
                    if j not in assigned:
                        u2 = units[j]
                        if abs(u1.occu + u2.occu - 1) < 1e-5 and u1.composition == u2.composition:
                            potential_u2.append([u2, np.linalg.norm(u2.geoc - u1.geoc), j])
                potential_u2.sort(key=lambda x: x[1])
                u2_real, u2_dist, u2_j = potential_u2[0]
                if u2_dist > 1.5:
                    warnings.warn('at least one disunit is paired with another that is >1.5 A far away, unlikely')
                    # raise DisUnitError('at least one disunit is paired with another that is >1.5 A far away, highly unlikely')
                pairs.append([u1, u2_real])
                assigned += [i, u2_j]
        return pairs, inv_conf

    @staticmethod
    def partition_asymmetric_unit(psites):
        """
        all disorder units in a flat list within in one asym unit

        :param psites:
        :return:
        """
        inv_conf = []
        disordered_sites = []
        for ps in psites:
            if abs(ps.properties['occu'] - 1) > 1e-3:
                disordered_sites.append(ps)
            else:
                inv_conf.append(ps)
        # first group by occu
        group_by_occu = [list(v) for l, v in groupby(sorted(disordered_sites, key=lambda x: x.properties['occu']),
                                                     lambda x: x.properties['occu'])]
        # within one group, group again by connection
        units = []
        for i in range(len(group_by_occu)):
            group = group_by_occu[i]
            mols, unwrap_str_sorted, unwrap_pblock_list = PBCparser.unwrap(Structure.from_sites(group))
            """
            the problem here is it thinks a(a')-b-c(c')-d(d') has two pairs of disunit, as there is no disorder at b
            if a pblock is not far away from another, they should be one pblock
            """
            connected_blocks = DisUnit.get_connected_pblock(unwrap_pblock_list)
            for pblock in connected_blocks:
                units.append(DisUnit(pblock))
        return units, inv_conf

    @staticmethod
    def get_connected_pblock(pblocks, cutoff=4.0):
        from scipy.spatial.distance import cdist
        import networkx as nx

        def get_coords(pb: [PeriodicSite]):
            coords = np.zeros((len(pb), 3))
            for i in range(len(pb)):
                coords[i] = pb[i].coords
            return coords

        def get_shortest_distance_between_blocks(pb1, pb2):
            pb1coords = get_coords(pb1)
            pb2coords = get_coords(pb2)
            distmat = cdist(pb1coords, pb2coords)
            minid = np.unravel_index(np.argmin(distmat, axis=None), distmat.shape)
            return distmat[minid]

        def get_block_graph(distmat, cutoff):
            g = nx.Graph()
            for i in range(len(distmat)):
                g.add_node(i)
                for j in range(i + 1, len(distmat)):
                    g.add_node(j)
                    if distmat[i][j] < cutoff:
                        g.add_edge(i, j)
            return g

        distmat_pblocks = np.zeros((len(pblocks), len(pblocks)))
        for i in range(len(pblocks)):
            for j in range(i + 1, len(pblocks)):
                distmat_pblocks[i][j] = get_shortest_distance_between_blocks(pblocks[i], pblocks[j])
                distmat_pblocks[j][i] = distmat_pblocks[i][j]
        block_graph = get_block_graph(distmat_pblocks, cutoff)
        connected_block_ids = nx.connected_components(block_graph)  # [[1,3], [2, 4, 5], ...]
        merged_blocks = []
        for ids in connected_block_ids:
            merged_block = []
            for i in ids:
                merged_block += pblocks[i]
            merged_blocks.append(merged_block)
        return merged_blocks


class ConfigConstructor:

    def __init__(self, inv_conf: [PeriodicSite], disunit_pairs):
        self.inv_conf = inv_conf
        self.disunit_pairs = disunit_pairs

    @staticmethod
    def gen_instructions(disunit_pairs, n_asym, n_cell):
        """
        get all possible instructions, notice this is exponentially scaled

        # of configs = 2^len(self.disunitpairs)^n_asym^n_cell where 2 comes from pairwise disordered occupancies

        return a list of instruction to build a config

        instruction is a nested tuple, instruction[icell][iasym][ipair] gives idis,
        that is, instruction[1][2][3] == 4 means in the 2nd unitcell (icell=1),
        the 3rd asymm unit (iasym=2),
        the 4th disorder portion (ipair=3),
        the 5th disorder unit (disunit_pairs[3][4], idis=4) is going to present

        :param disunit_pairs:
        :param n_asym:
        :param n_cell:
        """
        pair_size = []
        for ipair in range(len(disunit_pairs)):
            pair_size.append(len(disunit_pairs[ipair]))  # should always be 2
        dis_ins_combinations = itertools.product(*[list(range(n)) for n in pair_size])
        results = itertools.product(dis_ins_combinations, repeat=n_asym)
        results = list(itertools.product(results, repeat=n_cell))
        return results

    @staticmethod
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

    @staticmethod
    def dissc_to_config(sc, disunit_pairs, instruction):
        """
        take instructions to generate a certain config from super cell

        :param Structure sc: supercell structure
        :param disunit_pairs:
        :param instruction:
        :return:
        """
        pool = []
        config_occu = 1
        for icell in range(len(instruction)):
            for iasym in range(len(instruction[icell])):
                for ipair in range(len(instruction[icell][iasym])):
                    idis = instruction[icell][iasym][ipair]
                    disu = disunit_pairs[ipair][idis]
                    config_occu *= disu.occu
                    for label in disu.labels:
                        pool.append((label.label, iasym, icell))
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
