from collections import OrderedDict
from itertools import groupby

from pymatgen.core.sites import Site
from pymatgen.core.structure import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.structure import lattice_points_in_supercell

from ocelot.routines.disparser_functions import *
from ocelot.routines.pbc import PBCparser
from ocelot.schema.conformer import ConformerInitError
from ocelot.schema.conformer import MolConformer

"""
DisParser: prepare_data cif file into a list of disorder-free configurations

## notions and terms

- cifstring: 
    the raw string of the cif file
    
- configuration: 
    a disorder-free structure describing the crystal, 
    it has an occupancy,
    it could be a supercell
    
- unit_configuration: 
    a configuration that is not a supercell
    
- asymmetric_unit: 
    the atomic sites specified in cifstring, 
    please notice that the real, crystallographic asymmetric unit may be a subset of asymmetric_unit defined here

- disordered_sites:
    every atomic site with disorder

- counterpart:
    if atomic site X cannot co-exist with Y, Y is a counterpart of X and vice versa

- disorder_group:
    the largest subset of disordered_sites, 
    the atomic sites in disorder_unit can co-exist,
    *usually* there are two disorder_groups in a disordered cif file
    
- disorder_group_label: 
    an atomic site may be assigned a string to indicate the disorder_group it belongs to 
    defaults are ".", "1", "2".
    
- disorder_unit:
    a subset of a disorder_group
    each atomic site in the disorder_unit must share the same occupancy
    the sites in disorder_unit must be chemically connected/related (subject to a cutoff).

- disorder_portion:
    a list of disorder_units, each of them belongs to a different disorder_group
    in a certain configuration only one of the disorder_units will present to represent the disorder_portion

## disorder represenstation

- situation a: 

    the cif file may contain the following disorder-related fields

    _atom_site_occupancy 
    _atom_site_disorder_group 
    
    if both of the two keys exist, then we can extract configuration from the cif file with explicit occupancies
    if only **_atom_site_disorder_group** exists, **_atom_site_occupancy** will be set to 0.5 for disordered sites
    
- situation b:

    if not, it is still possible to prepare_data the cif file contains disorder, 
    this depends on the the values of **_atom_site_label**, 
    btw, '_atom_site_label' may not exist when using jmol to write cif

    the values of **_atom_site_label** form a list of *unique* labels, a label can be considered as a concatenation of 
    
    **E (an element symbol) + I (an index) + T (the third tag)**
    
    while CSD claims the alternative configurations should be labelled 
    with a suffix question mark (?), this is certainly not true for e.g. ASIXEH
    
#### classifications of situation b. cif files (that I've seen)

    1. EI is unique for each site, 
        only two disorder_groups
        sites of one disorder_group (the minor configuration) share the same T
        e.g. ALOVOO
        
    2. EI is **not** unique for each site, 
        only two disorder_groups
        sites of one disorder_group (the minor configuration) share the same T
        a disordered site has only one counterpart, and they share the same EI
        e.g. x17059
    
    3. EI is **not** unique for each site, 
        there is no site with a non-word T 
        Sites sharing the same EI representing disorder, 
        e.g. ASIXEH
    
    4. disordered sites are simply ignored in the cif file 
        s.t. the molecule is not complete, 
        e.g. ABEGET01
    
    5. disordered sites are simply ignored in the cif file, 
        but the molecule looks legit, 
        e.g. AGAJUO
    
    6. disordered sites are present in the cif file 
        but not labeled at all, e.g. ANOPEA

    notes: in case 1. it is possible to get the occupancies via CSD API, not always tho.

I cannot find a way to process 3., 4., 5., 6. without checking deposition records 
or chemistry analysis, so this parser will ignore them. 
That is, there can be only one alternative configuration.

## convert b1/b2 to a.
The idea is to find disorder_groups in b1 and b2, 
then artificially set (average) **_atom_site_occupancy** and **_atom_site_disorder_group**, 
this should be done intentionally.


## parsing a.

1. pairing/subgrouping disordered sites in an asymmetric unit:

    i. for each disordered atomic site, 
        find its counterpart if the other site is disordered and share the same EI
        
    ii. if i. failed at leaset once, find the nearest disordered site that share the same element

2. find disorder units:

    within one disorder_group, get connected components as disorder_units, 
    then find their pairs to form disorder_portion

3. config generation:
    
    - config_asymm = asymm_inv + asymm_disorder_portion1[0] + asymm_disorder_portion2[1] + ...
    - config_cell = config_asymm1 + config_asymm2 + ... 
    - config = config_cell1 + config_cell2 + ...
    
    the instruction for generating a configuration is a dictionary:
    
        conf_instruction[icell][iasymm][idisorder_portion] = j

    configuration as a cif file will be written with keys:
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

Ref:
[CSD API convention](https://downloads.ccdc.cam.ac.uk/documentation/API/descriptive_docs/crystal.html#disorder)  

"""


class DisorderParserError(Exception): pass




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
        try:
            self.index = re.findall(r"^\d+", tmplabel)[0]
        except IndexError:
            self.index = "0"

        tmplabel = tmplabel.lstrip(str(self.index))
        self.tag = tmplabel

        # if len(self.tag) > 1:
        #     raise ValueError('tag for {} is {}, this is unlikely'.format(label, self.tag))

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

    @staticmethod
    def get_labels_with_tag(tag, als):
        sametag = []
        for alj in als:
            if tag == alj.tag:
                sametag.append(alj)
        return sametag

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

    @staticmethod
    def get_psite_by_atomlable(psites: [PeriodicSite], al):
        for s in psites:
            if s.properties['label'] == str(al):
                return s
        raise ValueError('cannot find psite with atomlable {}'.format(str(al)))

    @classmethod
    def from_psite(cls, s: PeriodicSite):
        return cls(s.properties['label'])


class AsymmUnit:

    def __init__(self, psites: [PeriodicSite]):
        """
        an asymmetric unit with paired disorder units

        :param psites:
        # :param supress_sidechain_disorder: can be a function takes a psite and returns bool,
        #     if True the psite will always be considered as a non-disordered site
        """
        self.sites = psites
        self.labels = [AtomLabel(s.properties['label']) for s in self.sites]
        self.symbols = [l.element for l in self.labels]
        self.composition = tuple(sorted(self.symbols))

        # make sure occu and label and disg are present in prop
        for s in self.sites:
            for k in ['occu', 'label', 'disg']:
                if k not in s.properties.keys():
                    raise KeyError('{} is not in the properties of site {}!'.format(k, s))

        self.inv = []
        self.disordered_sites = []
        for ps in self.sites:
            if ps.properties['disg'] == '.' or abs(ps.properties['occu'] - 1) < 1e-5:
                self.inv.append(ps)
            else:
                self.disordered_sites.append(ps)

        # group by disg
        disordered_groups = [list(v) for l, v in
                             groupby(sorted(self.disordered_sites, key=lambda x: x.properties['disg']),
                                     lambda x: x.properties['disg'])]

        self.disordered_groups = {}
        # within one group, group again by connection
        for i in range(len(disordered_groups)):
            units = []
            group = disordered_groups[i]
            disg_label = group[0].properties['disg']

            block_list = find_connected_psites(group)
            """
            the problem here is it thinks a(a')-b-c(c')-d(d') has two pairs of disunit, as there is no disorder at b
            if a pblock is not far away from another, they should be one pblock
            """
            connected_blocks = get_connected_pblock(block_list, cutoff=4.0)
            for pblock in connected_blocks:
                units.append(DisorderUnit(pblock, disg=disg_label))
            self.disordered_groups[disg_label] = units
        # pairing
        self.disordered_portions = []
        init_disg_label = sorted(list(self.disordered_groups.keys()))[0]
        for u1 in self.disordered_groups[init_disg_label]:
            portion = [u1]
            for k in self.disordered_groups.keys():
                if k != init_disg_label:
                    u2 = DisorderUnit.find_counterpart(u1, self.disordered_groups[k])
                    portion.append(u2)
            self.disordered_portions.append(portion)

    def infodict(self, is_sidechain=None):
        au_labels = self.labels
        au_inv_labels = [AtomLabel(ps.properties['label']) for ps in self.inv]
        disordered_portions = []
        if is_sidechain is None:
            disordered_portions = [[u for u in portion] for portion in self.disordered_portions]
        else:
            for portion in self.disordered_portions:
                p = []
                for disorderunit in portion:
                    if all(is_sidechain(l) for l in disorderunit.labels):
                        au_inv_labels += disorderunit.labels
                        break
                    else:
                        p.append(disorderunit)
                if len(p) != 0:
                    disordered_portions.append(p)
        return au_labels, au_inv_labels, disordered_portions


class DisorderUnit:
    def __repr__(self):
        return "DisorederUnit -- disg: {}\n".format(self.disg) + " ".join([str(l) for l in self.labels])

    def __str__(self):
        return self.__repr__()

    def __eq__(self, other):
        return set(self.labels) == set(other.labels)

    def __len__(self):
        return len(self.sites)

    def __hash__(self):
        return hash(set(self.labels))

    def is_likely_counterpart(self, other):
        if not isinstance(other, DisorderUnit):
            return False
        if self == other:
            return False
        if self.composition != other.composition:
            return False
        if set(self.eis) == set(other.eis):
            return True
        if abs(self.occu + other.occu - 1) < 1e-5:
            return True

    @staticmethod
    def find_counterpart(u1, u2s):
        """

        :type u1: DisorderUnit
        :type u2s: [DisorderUnit]
        """
        k = lambda x: np.linalg.norm(x.geoc - u1.geoc)
        potential_u2 = sorted(u2s, key=k)
        potential_u2 = [u2 for u2 in potential_u2 if u1.is_likely_counterpart(u2)]
        try:
            u2 = potential_u2[0]
            if np.linalg.norm(u2.geoc - u1.geoc) > 1.5:
                warnings.warn('at least one disunit is paired with another that is >1.5 A far away, unlikely')
            return u2
        except IndexError:
            raise DisorderParserError('cannot find counterpart for DisorderUnit: {}'.format(u1))

    def __init__(self, psites: [PeriodicSite], disg: str):
        """
        a portion of an asymmetric unit representing one possibility

        a pair of DisUnit is the basis of disorder, assuming maximal entropy
        """
        self.disg = disg
        self.sites = psites
        self.labels = [AtomLabel(s.properties['label']) for s in self.sites]
        self.eis = [al.ei for al in self.labels]
        self.symbols = [l.element for l in self.labels]
        self.composition = tuple(sorted(self.symbols))
        occus = [s.properties['occu'] for s in self.sites]
        disgs = [s.properties['disg'] for s in self.sites]
        if len(set(occus)) != 1:
            raise DisorderParserError('occu is not uniform for {}'.format(str(self)))
        if len(set(disgs)) != 1:
            raise DisorderParserError('disg is not uniform for {}'.format(str(self)))
        self.occu = occus[0]
        self.disg = disgs[0]

    @property
    def geoc(self):
        """
        in cart
        """
        c = np.zeros(3)
        for s in self.sites:
            c += s.coords
        return c / len(self.sites)


class DisParser:  # chaos parser sounds cooler?

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
            data[atomlabel] = [x, y, z, symbol, occu, disgrp] this will be used to write config cif file
        """
        # prepare_data into cifdata
        self.classification = None
        self.cifstring = cifstring
        self.identifier, self.cifdata = get_pmg_dict(self.cifstring)
        self.cifdata['_atom_site_fract_x'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_x']]
        self.cifdata['_atom_site_fract_y'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_y']]
        self.cifdata['_atom_site_fract_z'] = [braket2float(x) for x in self.cifdata['_atom_site_fract_z']]

        for i in range(len(self.cifdata['_atom_site_type_symbol'])):
            if self.cifdata['_atom_site_type_symbol'][i] == 'D':
                warnings.warn('D is considered as H in the ciffile _atom_site_type_symbol!')
                self.cifdata['_atom_site_type_symbol'][i] = 'H'

        # check cif file
        try:
            labels = self.cifdata['_atom_site_label']
        except KeyError:
            raise CifFileError('no _atom_site_label field in the cifstring!')
        if len(labels) != len(set(labels)):
            warnings.warn('duplicate labels found in the cifstring, reassign labels!')
            for i in range(len(self.cifdata['_atom_site_label'])):
                label = AtomLabel(self.cifdata['_atom_site_label'][i])
                self.cifdata['_atom_site_label'][i] = label.element + str(i) + label.tag
            # raise CifFileError('duplicate labels found in the cifstring!')

        # "global" info
        self.labels = [AtomLabel(lab) for lab in labels]
        self.eis = [al.ei for al in self.labels]
        self.tags = [al.tag for al in self.labels]
        self.latparams = [self.cifdata[k] for k in latt_labels]
        self.lattice = Lattice.from_parameters(*[braket2float(p) for p in self.latparams], True)

        if '_atom_site_disorder_group' in self.cifdata.keys():
            self.was_fitted = True

            # deal with e.g. k06071
            for i in range(len(self.cifdata['_atom_site_disorder_group'])):
                if self.cifdata['_atom_site_disorder_group'][i] != ".":
                    dv = self.cifdata['_atom_site_disorder_group'][i]
                    dv = abs(int(dv))
                    self.cifdata['_atom_site_disorder_group'][i] = dv

            disgs = self.cifdata['_atom_site_disorder_group']

            disg_vals = list(set([disg for disg in disgs if disg != "."]))
            disg_vals.sort()
            for i in range(len(self.cifdata['_atom_site_disorder_group'])):
                for j in range(len(disg_vals)):
                    dv = disg_vals[j]
                    if self.cifdata['_atom_site_disorder_group'][i] == dv:
                        self.cifdata['_atom_site_disorder_group'][i] = str(j + 1)
                        break
            if '_atom_site_occupancy' not in self.cifdata.keys():
                occus = []
                for disg in self.cifdata['_atom_site_disorder_group']:
                    if disg == '.':
                        occus.append(1)
                    else:
                        if int(disg) & 1:
                            occus.append(0.51)
                        else:
                            occus.append(0.49)
                self.cifdata['_atom_site_occupancy'] = occus
        else:
            self.was_fitted = False

        if self.was_fitted:
            self.cifdata['_atom_site_occupancy'] = [braket2float(x) for x in self.cifdata['_atom_site_occupancy']]
            # self.cifdata['_atom_site_disorder_group'] = [braket2float(x) for x in
            #                                              self.cifdata['_atom_site_disorder_group']]

        # prepare self.data to be used in parsing
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

    @classmethod
    def from_ciffile(cls, fn):
        with open(fn, 'r') as f:
            s = f.read()
        return cls(s)

    def classify(self):
        """
        one cif file belongs to one of the following categories:

            dis-alpha: was fitted with, at least, disorder group

            nodis-0: no dup in eis, one type of unique tag

            dis-beta: dup in eis, EIx -- EIy, x or y can be empty, e.g. c17013.cif

            dis-gamma: no dup in eis, EIx -- EJy, x or y can be empty, e.g. ALOVOO.cif

            weird: else

            note x17059.cif has hydrogens like H12A -- H12D, this can only be captured
            by previously fitted disg and occu (dis-a), in general there
            should be a check on whether disgs identified by the parser are identical to previously fitted
        """
        if self.was_fitted:
            return "dis-alpha"

        # we look at labels of non-hydrogen sites
        labels = [l for l in self.labels if l.element != 'H']
        tags = [l.tag for l in labels]
        eis = [l.ei for l in labels]
        unique_tags = list(set(tags))
        unique_eis = list(set(eis))

        if len(unique_tags) == 1 and len(unique_eis) == len(eis):
            for al in self.labels:
                self.data[al][4] = 1
                self.data[al][5] = '.'
            return "nodis-0"

        u_tags = list(set(self.tags))
        u_nonalphabet_nonempty_tags = [t for t in u_tags if not re.search(r"[a-zA-Z]", t) and t != ""]
        if len(u_nonalphabet_nonempty_tags) == 1:
            minor_tag = u_nonalphabet_nonempty_tags[0]
            minor_als = AtomLabel.get_labels_with_tag(minor_tag, self.labels)
            disg2 = []
            disg1 = []
            classification = ""
            # try pairing based on ei, EIx -- EIy
            for minor_al in minor_als:
                potential_aljs = [l for l in self.labels if l.ei == minor_al.ei and l != minor_al]
                if len(potential_aljs) == 0:
                    classification = "dis-gamma"
                    break
                potential_aljs.sort(key=lambda x: x.tag)
                alj = potential_aljs[0]
                disg1.append(alj)
                disg2.append(minor_al)
            if classification == "":
                classification = "dis-beta"
            else:
                # EIx -- EJy
                disg2 = []
                disg1 = []
                cutoff = 2.5
                warnings.warn(
                    'W: trying to prepare_data disorder in dis-gamma with cutoff {}, this is unreliable!'.format(
                        cutoff))
                major_als = [l for l in self.labels if l not in minor_als]
                assigned_major = []
                for minor_al in minor_als:
                    alj = self.get_nearest_label(minor_al, [l for l in major_als if l not in assigned_major],
                                                 self.lattice, cutoff)
                    assigned_major.append(alj)
                    disg1.append(alj)
                    disg2.append(minor_al)
                classification = "dis-gamma"
            inv = [l for l in self.labels if l not in disg1 + disg2]
            for al in inv:
                self.data[al][4] = 1
                self.data[al][5] = '.'
            for al in disg1:
                self.data[al][4] = 0.51
                self.data[al][5] = '1'
            for al in disg2:
                self.data[al][4] = 0.49
                self.data[al][5] = '2'
            return classification

        if len(unique_tags) > 0 and len(eis) != len(unique_eis):
            # try pairing based on ei, EIx -- EIy
            for ei in self.eis:
                als = [l for l in self.labels if l.ei == ei]
                if len(als) == 1:
                    al = als[0]
                    self.data[al][4] = 1
                    self.data[al][5] = '.'
                else:
                    als.sort(key=lambda x: x.tag)
                    ali = als[0]
                    alj = als[1]
                    self.data[ali][4] = 0.51
                    self.data[ali][5] = '1'
                    self.data[alj][4] = 0.49
                    self.data[alj][5] = '2'
            return "dis-beta"

        raise NotImplementedError("this structure is classified as weird, better to check ciffile")

    def get_nearest_label(self, ali: AtomLabel, neighbor_labels: [AtomLabel], lattice: Lattice, cutoff=2.5):
        """
        given a list of potential neighbouring AtomLable, get the one that is closest and has the same element
        """
        data = self.data
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
            raise ValueError(
                'cannot find any same-symbol neighbors of label {}'.format(ali.label))
        ali_nbs_dictance = sorted(ali_nbs_dictance, key=lambda x: x[1])
        nnb, nnbdis = ali_nbs_dictance[0]
        if nnbdis > cutoff:
            raise ValueError(
                'cannot find any same-symbol neighbors of label {} within {}'.format(
                    ali.label, cutoff))
        return nnb

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

    def prepare_data(self):
        self.classification = self.classify()
        print('{} thinks this cif file belongs to class {}'.format(self.__class__.__name__, self.classification))
        # from pprint import pprint
        # pprint(self.data)

    @staticmethod
    def get_vanilla_configs(psites: [PeriodicSite]):
        """
        must have 'disg' in properties
        """
        disg_vs = []
        for i in range(len(psites)):
            disg = psites[i].properties['disg']
            disg_vs.append(disg)
        disg_vs = list(set(disg_vs))
        disg_vs.sort()
        if len(disg_vs) == 1:
            return [Structure.from_sites(psites)]
        disgs = []
        inv = [s for s in psites if s.properties['disg'] == disg_vs[0]]
        for disg in disg_vs[1:]:
            disg_sites = []
            for s in psites:
                if s.properties['disg'] == disg:
                    disg_sites.append(s)
            disgs.append(Structure.from_sites(disg_sites + inv))
        return disgs

    @staticmethod
    def get_site_location_by_key(psites: [PeriodicSite], key='label'):
        """
        loc[key] = 'bone'/'sidechain'
        must have 'imol' 'disg' 'siteid' and <key> assigned
        """
        res = {}
        k = lambda x: x.properties['imol']
        psites.sort(key=k)
        for imol, group in groupby(psites, key=k):
            obc_sites = []
            for ps in group:
                obc_sites.append(Site(ps.species_string, ps.coords, properties=ps.properties))
            try:
                molconformer = MolConformer.from_sites(obc_sites, siteids=[s.properties['siteid'] for s in obc_sites])
            except ConformerInitError:
                raise DisorderParserError('cannot init legit conformer to identify site location')
            for sid in molconformer.siteids:
                site = molconformer.get_site_byid(sid)
                if molconformer.backbone is None:
                    res[site.properties[key]] = 'sidechain'
                else:
                    if sid in molconformer.backbone.siteids:
                        res[site.properties[key]] = 'bone'
                    else:
                        res[site.properties[key]] = 'sidechain'
        return res

    def to_configs(self, write_files=False, scaling_mat=(1, 1, 1), assign_siteids=True, vanilla=True,
                   supressdisorder=True):
        """
        return
            pstructure, pmg structure is the unit cell structure with all disordered sites
            unwrap_str, pmg unwrap structure, with all disordered sites
            mols, a list of pmg mol with disordered sties
            confs, [[conf1, occu1], ...], conf1 is a clean structure
        """
        self.prepare_data()  # edit self.data
        asym_psites = self.get_psites_from_data()  # assign fields: occu, disg, label
        raw_symmops = get_symmop(self.cifdata)

        psites, symmops = apply_symmop(asym_psites, raw_symmops)  # assign field: iasym
        pstructure = Structure.from_sites(psites, to_unit_cell=True)
        if write_files:
            pstructure.to('cif', 'dp_wrapcell.cif')

        # assign siteids here if unit cell
        if np.prod(scaling_mat) == 1.0 and assign_siteids:
            print('siteid is assigned to unitcell in dp.to_configs')
            for isite in range(len(pstructure)):
                pstructure[isite].properties['siteid'] = isite

        # sc, n_unitcell = ConfigConstructor.build_supercell_full_disorder(pstructure, scaling_mat)
        mols, unwrap_str, unwrap_pblock_list = PBCparser.unwrap_and_squeeze(pstructure)  # assign field: imol

        if write_files:
            unwrap_str.to('cif', 'dp_unwrapunit.cif')

        sc, n_unitcell = ConfigConstructor.build_supercell_full_disorder(unwrap_str, scaling_mat)  # assign field: icell

        if write_files:
            sc.to('cif', 'dp_supercell.cif')

        if np.prod(scaling_mat) != 1.0 and assign_siteids:
            warnings.warn(
                'siteid is assigned to supercell in dp.to_configs, notice in this case, sites in mols, unwrap_str, unwrap_pblock_list mols do not have siteids!'
            )
            for isite in range(len(sc)):
                sc[isite].properties['siteid'] = isite

        if vanilla:
            iconf = 0
            conf_structures = self.get_vanilla_configs(sc.sites)
            # # it is possible to have the situation where unwrapping disordered sites lead to super-large molecule
            # # e.g. x15029, this can be dealt with unwrpa again, I do not think this is common tho.
            # # if use this you may want to reassign imol field
            # for c in conf_structures:
            #     _, struct, _ = PBCparser.unwrap(c)
            #     unwrap_again.append(struct)
            # confs = [[c, 0.5] for c in conf_structures]
            def get_average_occu(structure:Structure):
                try:
                    return np.mean([s.properties['occu'] for s in structure.sites])
                except KeyError:
                    return 0.5
            confs = [[c, get_average_occu(c)] for c in conf_structures]
            if write_files:
                for conf in conf_structures:
                    conf.to('cif',
                            'vconf_{}.cif'.format(iconf))  # pymatgen somehow does not write disg field in the cif
                    iconf += 1
            if len(set(c.composition for c in conf_structures)) > 1:
                warnings.warn('different compositions in confs! better use vconf_0!')
                # raise DisorderParserError('different compositions in confs!')
            return pstructure, unwrap_str, mols, sorted(confs, key=lambda x: x[1], reverse=True)

        if supressdisorder:
            conf_structures = self.get_vanilla_configs(sc.sites)
            locbysitelabel = {}
            for conf in conf_structures:
                try:
                    loc_conf = self.get_site_location_by_key(conf.sites, 'label')
                except DisorderParserError:
                    warnings.warn('cannot get site location, supressdisorder is disabled')
                    loc_conf = {}
                for k in loc_conf:
                    locbysitelabel[k] = loc_conf[k]

            def is_sidechain(label: AtomLabel):
                if locbysitelabel[str(label)] == 'sidechain':
                    return True
                return False
        else:
            is_sidechain = None

        iconf = 0
        confs = []
        conf_structures = []
        au = AsymmUnit(asym_psites)
        au_labels, au_inv_labels, disordered_portions = au.infodict(is_sidechain)
        inv_labels = [l.label for l in au_inv_labels]
        conf_ins = ConfigConstructor.gen_instructions(disordered_portions, len(symmops), n_unitcell)
        for confin in conf_ins:
            conf, conf_occu = ConfigConstructor.dissc_to_config(sc, inv_labels, disordered_portions, confin)
            confs.append([conf, conf_occu])
            conf_structures.append(conf)
            if write_files:
                conf.to('cif', 'conf_{}.cif'.format(iconf))  # pymatgen somehow does not write disg field in the cif
            iconf += 1
        if len(set(c.composition for c in conf_structures)) > 1:
            warnings.warn('different compositions in confs! better use vconf_0!')
            # raise DisorderParserError('different compositions in confs!')
        return pstructure, unwrap_str, mols, sorted(confs, key=lambda x: x[1], reverse=True)


class ConfigConstructor:

    @staticmethod
    def gen_instructions(asym_portions, n_asym, n_cell):
        """
        get all possible instructions, notice this is exponentially scaled

        # of configs = 2^n_port*n_asym*n_cell where 2 comes from pairwise disordered occupancies
        """
        pair_size = []
        for ipair in range(len(asym_portions)):
            pair_size.append(len(asym_portions[ipair]))  # should always be 2
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
    def dissc_to_config(sc: Structure, inv_labels: [str], disportions, instruction):
        """
        take instructions to generate a certain config from super cell
        """
        pool = []
        config_occu = 1
        for icell in range(len(instruction)):
            for iasym in range(len(instruction[icell])):
                for iport in range(len(instruction[icell][iasym])):
                    idis = instruction[icell][iasym][iport]
                    disu = disportions[iport][idis]
                    config_occu *= disu.occu
                    for label in disu.labels:
                        pool.append((label.label, iasym, icell))
        conf_sites = []
        for ps in sc.sites:
            prop = ps.properties
            if prop['label'] in inv_labels:
                conf_sites.append(ps)
            else:
                plabel = prop['label']
                piasym = prop['iasym']
                picell = prop['icell']
                if (plabel, piasym, picell) in pool:
                    conf_sites.append(ps)
        return Structure.from_sites(conf_sites), config_occu
