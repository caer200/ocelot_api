import math
from ocelot.schema.element import Element
from ocelot.schema.omol import OMol
from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.core.structure import Site, PeriodicSite, IMolecule, IStructure, Molecule, Structure
from pymatgen.io.cif import CifFile

"""
CIFparser: parse cif file into a list of configurations with no disorder
PBCparser: get unwrapped structure and mols
"""


class CIFparser:

    def __init__(self, cifstring, identifier, pymatgen_dict, clean_dicts, clean_cifstrings, nconfig, is_disorder):

        self.cifstring = cifstring
        self.identifier = identifier
        self.pymatgen_dict = pymatgen_dict
        self.clean_dicts = clean_dicts
        self.clean_cifstrings = clean_cifstrings
        self.nconfig = nconfig
        self.is_disorder = is_disorder

    @classmethod
    def from_cifstring(cls, cifstring):
        cifdata = CifFile.from_string(cifstring).data
        identifier = '_'.join(list(cifdata.keys()))
        pymatgen_dict = list(cifdata.items())[0][1].data
        is_disorder = CIFparser.isdisorder(pymatgen_dict)
        clean_dicts = CIFparser.split_pymatgen_dict(pymatgen_dict, is_disorder, identifier)
        clean_cifstrings = [CIFparser.get_config_dict_string(d) for d in clean_dicts]
        nconfig = len(clean_dicts)
        return cls(cifstring, identifier, pymatgen_dict, clean_dicts, clean_cifstrings, nconfig, is_disorder)

    @staticmethod
    def get_config_dict_string(d):
        s = d['identifier'] + '\n'
        s += 'loop_\n _symmetry_equiv_pos_as_xyz\n' + '\n'.join(["'" + xyz + "'" for xyz in d['symop']]) + '\n\n'
        s += '_cell_length_a\t' + d['a'] + '\n'
        s += '_cell_length_b\t' + d['b'] + '\n'
        s += '_cell_length_c\t' + d['c'] + '\n'
        s += '_cell_angle_alpha\t' + d['A'] + '\n'
        s += '_cell_angle_beta\t' + d['B'] + '\n'
        s += '_cell_angle_gamma\t' + d['C'] + '\n'
        s += 'loop_\n _atom_site_label\n _atom_site_type_symbol\n _atom_site_fract_x\n _atom_site_fract_y\n _atom_site_fract_z\n _atom_site_occupancy\n _atom_site_disorder_group\n'

        for site in d['msites']:
            line = '\t'.join(site)
            s += line + '\n'
        return s

    @staticmethod
    def write_config_dict(d, cifn):
        """
        write config_dict into minimum readable cif file
        :return:
        """
        s = CIFparser.get_config_dict_string(d)
        with open(cifn, 'w') as f:
            f.write(s)

    def get_clean_cifs_stringlist(self):
        res = []
        for i in range(len(self.clean_dicts)):
            cleandict = self.clean_dicts[i]
            res.append(self.get_config_dict_string(cleandict))
        return res

    def write_clean_cifs(self, prefix='cleanconfig'):
        fns = []
        for i in range(self.nconfig):
            fn = '{}-{}.cif'.format(prefix, i)
            self.write_config_dict(self.clean_dicts[i], fn)
            fns.append(fn)
        return fns

    @staticmethod
    def split_pymatgen_dict(pymatgen_dict, isdisorder, identifier):
        if not isdisorder:
            config_d = dict(
                identifier='data_{}_config-{}'.format(identifier, 0),
                a=pymatgen_dict['_cell_length_a'],
                b=pymatgen_dict['_cell_length_b'],
                c=pymatgen_dict['_cell_length_c'],
                A=pymatgen_dict['_cell_angle_alpha'],
                B=pymatgen_dict['_cell_angle_beta'],
                C=pymatgen_dict['_cell_angle_gamma'],
            )
            if '_space_group_symop_operation_xyz' in pymatgen_dict.keys():
                config_d['symop'] = pymatgen_dict['_space_group_symop_operation_xyz']
            elif '_symmetry_equiv_pos_as_xyz' in pymatgen_dict.keys():
                config_d['symop'] = pymatgen_dict['_symmetry_equiv_pos_as_xyz']
            sites = []
            for j in range(len(pymatgen_dict['_atom_site_label'])):
                site = [
                    pymatgen_dict['_atom_site_label'][j],
                    pymatgen_dict['_atom_site_type_symbol'][j],
                    pymatgen_dict['_atom_site_fract_x'][j],
                    pymatgen_dict['_atom_site_fract_y'][j],
                    pymatgen_dict['_atom_site_fract_z'][j],
                    '1.0',
                    '.'
                ]
                sites.append(site)
            config_d['msites'] = sites
            return [config_d]
        else:
            config_ds = []
            dg_labels = sorted(
                [line for line in list(set([s.strip('-') for s in pymatgen_dict['_atom_site_disorder_group']]))
                 if line != u'.'])
            for l in dg_labels:
                config_d = dict(
                    identifier='data_{}_config-{}'.format(identifier, l),
                    a=pymatgen_dict['_cell_length_a'],
                    b=pymatgen_dict['_cell_length_b'],
                    c=pymatgen_dict['_cell_length_c'],
                    A=pymatgen_dict['_cell_angle_alpha'],
                    B=pymatgen_dict['_cell_angle_beta'],
                    C=pymatgen_dict['_cell_angle_gamma'],
                )
                if '_space_group_symop_operation_xyz' in pymatgen_dict.keys():
                    config_d['symop'] = pymatgen_dict['_space_group_symop_operation_xyz']
                elif '_symmetry_equiv_pos_as_xyz' in pymatgen_dict.keys():
                    config_d['symop'] = pymatgen_dict['_symmetry_equiv_pos_as_xyz']
                sites = []
                for j in range(len(pymatgen_dict['_atom_site_label'])):
                    if pymatgen_dict['_atom_site_disorder_group'][j].strip('-') == l \
                            or pymatgen_dict['_atom_site_disorder_group'][j] == u'.':
                        site = [
                            pymatgen_dict['_atom_site_label'][j],
                            pymatgen_dict['_atom_site_type_symbol'][j],
                            pymatgen_dict['_atom_site_fract_x'][j],
                            pymatgen_dict['_atom_site_fract_y'][j],
                            pymatgen_dict['_atom_site_fract_z'][j],
                            '1.0',
                            pymatgen_dict['_atom_site_disorder_group'][j]
                        ]
                        sites.append(site)
                config_d['msites'] = sites
                config_ds.append(config_d)
            return config_ds

    @staticmethod
    def isdisorder(pymatgen_dict):
        if '_atom_site_disorder_group' in pymatgen_dict.keys():
            dg_labels = sorted(
                [line for line in list(set([s.strip('-') for s in pymatgen_dict['_atom_site_disorder_group']])) if
                 line != u'.'])
            if len(dg_labels) > 0:
                return True
        return False


class PBCparser:
    #
    # def __init__(self, pstructure):
    #     self.nsites = len(pstructure)
    #     self.structure = pstructure

    @staticmethod
    def get_dist_and_trans(lattice, fc1, fc2):
        """
        get the shortest distance and corresponding translation vector between two frac coords

        :param lattice: pmg lattic obj
        :param fc1:
        :param fc2:
        :return:
        """
        v, d2 = pbc_shortest_vectors(lattice, fc1, fc2, return_d2=True)
        fc = lattice.get_fractional_coords(v[0][0]) + fc1 - fc2
        return math.sqrt(d2[0, 0]), fc

    @staticmethod
    def unwrap(pstructure):
        """
        unwrap the structure, extract isolated mols

        :param pstructure: periodic structure obj from pymatgen
        :return:
        """
        psites = pstructure.sites
        pindices = range(len(psites))
        visited = []
        block_list = []
        unwrap = []
        unwrap_block_list = []
        unwrap_pblock_list = []
        while len(visited) != len(psites):
            # initialization
            unvisited = [idx for idx in pindices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            unwrap.append(psites[ini_idx])
            unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords)]
            unwrap_pblock = [psites[ini_idx]]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    distance, fctrans = PBCparser.get_dist_and_trans(pstructure.lattice,
                                                                     psites[block[pointer]]._frac_coords,
                                                                     psites[i]._frac_coords, )

                    cutoff = Element.covalent_radii[psites[block[pointer]].species_string] + Element.covalent_radii[
                        psites[i].species_string]
                    cutoff *= 1.3
                    if distance < cutoff:
                        block.append(i)
                        psites[i] = PeriodicSite(psites[i].species_string, psites[i]._frac_coords + fctrans,
                                                 pstructure.lattice)
                        unwrap_block.append(Site(psites[i].species_string, psites[i].coords))
                        unwrap.append(psites[i])
                        unwrap_pblock.append(psites[i])
                visited.append(block[pointer])
                pointer += 1
            unwrap_block_list.append(unwrap_block)
            unwrap_pblock_list.append(unwrap_pblock)
            block_list.append(block)
        mols = [IMolecule.from_sites(i) for i in unwrap_block_list]
        unwrap = sorted(unwrap, key=lambda x: x.species_string)
        unwrap_str_sorted = IStructure.from_sites(unwrap)
        return mols, unwrap_str_sorted, unwrap_pblock_list

    @staticmethod
    def squeeze(pstructure):
        """
        after unwrapping, the mols can be far away from each other, this tries to translate them s.t. they stay together

        :param pstructure:
        :return:
        """
        mols, unwrap_structure, psiteblocks = PBCparser.unwrap(pstructure)

        omols = []
        for mm in mols:
            omol = OMol.from_pymatgen_mol(mm)
            if not omol.is_solvent:
                omols.append(omol)

        if len(omols) > 1:
            refpoint = omols[0].backbone.geoc
            refpoint = unwrap_structure.lattice.get_fractional_coords(refpoint)
            for i in range(1, len(omols)):
                varmol = omols[i]
                varpoint = varmol.backbone.geoc
                varpoint = unwrap_structure.lattice.get_fractional_coords(varpoint)
                distance, fctrans = PBCparser.get_dist_and_trans(unwrap_structure.lattice, refpoint, varpoint)
                for j in range(len(psiteblocks[i])):
                    psiteblocks[i][j]._frac_coords += fctrans
            psites = []
            mols = []
            for pblock in psiteblocks:
                block = []
                for ps in pblock:
                    psites.append(ps)
                    block.append(Site(ps.species_string, ps.coords))
                mol = Molecule.from_sites(block)
                mols.append(mol)
            unwrap_structure = Structure.from_sites(sorted(psites, key=lambda x: x.species_string))

            omols = []
            for m in mols:
                omol = OMol.from_pymatgen_mol(m)
                if not omol.is_solvent:
                    omols.append(omol)
            return mols, omols, unwrap_structure

        return mols, omols, unwrap_structure
