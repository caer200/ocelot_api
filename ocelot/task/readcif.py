from ocelot.routines.disparser import DisParser
from ocelot.routines.pbc import Site, MolConformer
from itertools import groupby

class ReadCif:

    def __init__(self, cifstring):
        self.cifstring = cifstring
        dp = DisParser(self.cifstring)
        config_infos = dp.to_configs(write_files=False)  # if True writes conf_x.cif, configs is a list of pmg Structure
        self.config_structures = [item[0] for item in config_infos]
        self.occus = [item[1] for item in config_infos]

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
        insepect_data = {}
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
            molconformer.conformer_properties['disorder'] = mol_disorder_info
            insepect_data[imol] = molconformer
        return insepect_data
