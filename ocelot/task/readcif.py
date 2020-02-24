from itertools import groupby
from collections import OrderedDict

from ocelot.schema.configuration import Config
from ocelot.routines.disparser import DisParser
from ocelot.schema.conformer import MolConformer
from ocelot.routines.pbc import Site

"""
ReadCif implements a set of checkers/functions as the first step of reading cif file
1. is one type of molecule?
2. is the molecule legit? (can it be parsed to rdmol)
3. where is the disorder --> bone or side group or no
4. get all configurations (during which molconformers for each config will be obtained)
"""


class ReadCif:

    def __init__(self, cifstring):
        self.cifstring = cifstring
        self.dp = DisParser(self.cifstring)
        dis_pstructure, dis_unwrap_str, dis_mols, config_infos = self.dp.to_configs(write_files=False)  # if True writes conf_x.cif, configs is a list of pmg Structure
        self.disordered_pstructure = dis_unwrap_str
        self.disordered_pmgmols = dis_mols
        self.config_structures = []
        self.occus = []
        for item in config_infos:
            self.config_structures.append(item[0])
            self.occus.append(item[1])

        self.configs = []
        for i in range(len(self.config_structures)):
            structure = self.config_structures[i]
            self.configs.append(Config(structure, occu=self.occus[i], assign_siteids=False))

        self.properties = OrderedDict()
        self.properties['is_one_type_mol'] = all(len(c.molgraph_set()) == len(c.molconformers) for c in self.configs)
        self.properties['is_all_mol_legit'] = all(c.molconformers_all_legit() for c in self.configs)
        self.properties['where_is_disorder'] = self.where_is_disorder()

    def as_dict(self):
        d = OrderedDict()
        d['cifstring'] = self.cifstring
        d['clean_pstructures'] = [s.as_dict() for s in self.config_structures]
        d['occus'] = self.occus
        d['disordered_pmgmols'] = [m.as_dict() for m in self.disordered_pmgmols]
        d['disordered_pstructure'] = self.disordered_pstructure.as_dict()
        d['disparser'] = self.dp.as_dict()
        d['configs'] = [c.as_dict() for c in self.configs]
        d['properties'] = self.properties
        return d

    @classmethod
    def from_dict(cls, d):
        cifstring = d['cifstring']
        return cls(cifstring)


    @classmethod
    def from_ciffile(cls, ciffile):
        with open(ciffile, 'r') as f:
            s = f.read()
        return cls(s)

    def where_is_disorder(self):
        """
        data[imol] = MolConformer with disorder info in conformer_properties
        """
        c = self.config_structures[0]
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
            molconformer = MolConformer.from_sites(obc_sites, siteids=[s.properties['siteid'] for s in obc_sites])
            if len(disordered_siteid) > 0:
                if set(molconformer.backbone.siteids).intersection(set(disordered_siteid)):
                    mol_disorder_info = 'bone disorder'
                else:
                    mol_disorder_info = 'sc disorder'
            else:
                mol_disorder_info = 'no disorder'
            disorderinfo[imol] = mol_disorder_info
        return disorderinfo
