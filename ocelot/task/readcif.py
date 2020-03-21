from itertools import groupby
from collections import OrderedDict

from ocelot.schema.configuration import Config
from ocelot.routines.disparser import DisParser
from ocelot.schema.conformer import MolConformer, ConformerInitError
from ocelot.routines.pbc import Site

"""
ReadCif implements a set of checkers/functions as the first step of reading cif file
1. is one type of molecule?
2. is the molecule legit? (can it be parsed to rdmol)
3. where is the disorder --> bone or side group or no
4. get all configurations (during which molconformers for each config will be obtained)
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
        for i in range(len(config_structures)):
            structure = config_structures[i]
            configs.append(Config.from_labeled_clean_pstructure(structure, occu=config_occupancies[i]))
        self.results['configuraions'] = configs

        # these are checked against to configs[0]
        self.results['n_unique_molecule'] = len(configs[0].molgraph_set())
        self.results['n_molconformers'] = len(configs[0].molconformers)
        self.results['all_molconformers_legit'] = configs[0].molconformers_all_legit()
        self.results['disorder_location'] = self.where_is_disorder(config_structures[0])

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
    def where_is_disorder(config_structure):
        """
        data[imol] = disorder info in conformer_properties
        """
        c = config_structure
        dic = {}
        k = lambda x: x.properties['imol']
        psites = sorted(c.sites, key=k)
        disorderinfo = {}
        for imol, group in groupby(psites, key=k):
            dic[imol] = list(group)  # e.g. dic[0] is psites with the imol=0
            disordered_siteid = []
            obc_sites = []
            for ps in dic[imol]:
                if abs(ps.properties['occu'] - 1.0) > 1e-3:
                    disordered_siteid.append(ps.properties['siteid'])
                obc_sites.append(Site(ps.species_string, ps.coords, properties=ps.properties))
            try:
                molconformer = MolConformer.from_sites(obc_sites, siteids=[s.properties['siteid'] for s in obc_sites])
            except ConformerInitError:
                disorderinfo[imol] = 'not sure'
                continue
            if len(disordered_siteid) > 0:
                if set(molconformer.backbone.siteids).intersection(set(disordered_siteid)):
                    mol_disorder_info = 'bone disorder'
                else:
                    mol_disorder_info = 'sc disorder'
            else:
                mol_disorder_info = 'no disorder'
            disorderinfo[imol] = mol_disorder_info
        disorderinfo = dict((str(k), v) for k, v in disorderinfo.items())
        return disorderinfo
