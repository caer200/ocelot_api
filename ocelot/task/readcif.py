import warnings
from collections import OrderedDict

from pymatgen.core.composition import CompositionError
from pymatgen.core.structure import Composition

from ocelot.routines.disparser import DisParser
from ocelot.schema.configuration import Config
from ocelot.schema.conformer import MolConformer

"""
ReadCif implements a set of checkers/functions as the first step of reading cif file
1. is one type of molecule?
2. is the molecule legit? (can it be parsed to rdmol)
3. where is the disorder --> bone or side group or no
4. get all configurations (during which molconformers for each config will be obtained)
5. is there any hydrogen missing?
"""


class ReadCif:

    def __init__(self, cifstring, source, identifier=None):
        self.cifstring = cifstring
        self.source = source
        self.dp = DisParser(self.cifstring)
        if identifier is None:
            self.identifier = self.dp.identifier
        else:
            self.identifier = identifier
        self.lattice = self.dp.lattice
        self.was_fitted = self.dp.was_fitted
        self.disorder_class = self.dp.classification
        self.results = OrderedDict()

    def read(self):
        dis_pstructure, dis_unwrap_str, dis_mols, config_infos = self.dp.to_configs(write_files=False, vanilla=True)  # if True writes conf_x.cif, configs is a list of pmg Structure
        self.disorder_class = self.dp.classification
        self.results['disordered_pstructure'] = dis_unwrap_str
        self.results['disordered_pmgmols'] = dis_mols

        config_structures = []
        config_occupancies = []
        for item in config_infos:
            config_structures.append(item[0])
            config_occupancies.append(item[1])

        self.results['config_sturcutures'] = config_structures
        self.results['config_occupancies'] = config_occupancies


        configs = []
        missh = []
        for i in range(len(config_structures)):
            structure = config_structures[i]
            conf = Config.from_labeled_clean_pstructure(structure, occu=config_occupancies[i])
            config_missingh = False
            for conformer in conf.molconformers:
                if conformer.is_missing_hydrogen():
                    config_missingh = True
                    break
            if config_missingh:
                conf.pstructure.to('cif', '{}_mhconf_{}.cif'.format(self.identifier, i))
                warnings.warn('missing hydrogens in {}_mhconf_{}.cif'.format(self.identifier, i))
            configs.append(conf)
            missh.append(config_missingh)
        self.results['configurations'] = configs
        self.results['missingh'] = missh

        # these are checked against to configs[0]
        check_config = configs[0]
        try:
            self.results['n_unique_molecule'] = len(check_config.molgraph_set())
            self.results['n_molconformers'] = len(check_config.molconformers)
            self.results['all_molconformers_legit'] = check_config.molconformers_all_legit()
            self.results['disorder_location'] = self.where_is_disorder(check_config)
        except:
            warnings.warn('there are problems in readcif.results, some fileds will be missing!')
        
        try:
            comp = Composition(self.dp.cifdata['_chemical_formula_sum'])
            self.results['cif_sum_composition'] = comp
            if not all(self.results['cif_sum_composition'] == mc.composition for mc in check_config.molconformers):
                self.results['sum_composition_match'] = False
                print('cif_sum_composition: {}'.format(self.results['cif_sum_composition']))
                for mc in check_config.molconformers:
                    print('mc composition: {}'.format(mc.composition))
                warnings.warn('moiety sum composition does not match that specified in cif file!')
            else:
                self.results['sum_composition_match'] = True
        except (KeyError, CompositionError) as e:
            self.results['cif_sum_composition'] = None
            self.results['sum_composition_match'] = None

        try:
            comp_str = self.dp.cifdata['_chemical_formula_moiety']
            comps = [Composition(s) for s in comp_str.split(',')]
            comps = sorted(comps, key=lambda x:len(x), reverse=True)
            if len(comps) > 1:
                warnings.warn('more than 1 moiety from cif file! only the largest one is checked!')
            self.results['cif_moiety_composition'] = comps[0]
            if not all(self.results['cif_moiety_composition'] == mc.composition for mc in check_config.molconformers):
                self.results['moiety_composition_match'] = False
                print('cif_moiety_composition: {}'.format(self.results['cif_moiety_composition']))
                for mc in check_config.molconformers:
                    print('mc composition: {}'.format(mc.composition))
                warnings.warn('moiety composition does not match that specified in cif file!')
            else:
                self.results['moiety_composition_match'] = True
        except (KeyError, CompositionError) as e:
            self.results['cif_moiety_composition'] = None
            self.results['moiety_composition_match'] = None

    # def as_dict(self):
    #     d = OrderedDict()
    #     d['cifstring'] = self.cifstring
    #     d['clean_pstructures'] = [s.as_dict() for s in self.config_structures]
    #     d['occus'] = self.occus
    #     d['disordered_pmgmols'] = [m.as_dict() for m in self.disordered_pmgmols]
    #     d['disordered_pstructure'] = self.disordered_pstructure.as_dict()
    #     d['disparser'] = self.dp.as_dict()
    #     d['configs'] = [c.as_dict() for c in self.configs]
    #     d['properties'] = self.properties
    #     return d
    #
    # @classmethod
    # def from_dict(cls, d):
    #     cifstring = d['cifstring']
    #     return cls(cifstring)


    @classmethod
    def from_ciffile(cls, ciffile, source, identifier=None):
        with open(ciffile, 'r') as f:
            s = f.read()
        return cls(s, source, identifier)

    @staticmethod
    def where_is_disorder(c: Config):
        """
        data[imol] = disorder info in conformer_properties
        """
        disorderinfo = {}
        mc: MolConformer
        for imol in range(len(c.molconformers)):
            mc = c.molconformers[imol]
            try:
                disordered_siteid = [s for s in mc if abs(s.properties['occu'] - 1) > 1e-3]
            except KeyError:
                warnings.warn('not all sites have occu field, cannot decide disorder location!')
                disorderinfo[imol] = 'not sure'
                continue
            if len(disordered_siteid) == 0:
                disorderinfo[imol] = 'no disorder'
            else:
                if mc.backbone is None:
                    disorderinfo[imol] = 'sc disorder'
                elif set(mc.backbone.siteids).intersection(set(disordered_siteid)):
                    disorderinfo[imol] = 'bone disorder'
                else:
                    disorderinfo[imol] = 'sc disorder'
        disorderinfo = dict((str(k), v) for k, v in disorderinfo.items())
        return disorderinfo
