import hashlib  # 3.8
from abc import ABCMeta
from abc import abstractmethod

import numpy as np
from monty.json import MSONable
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure


def sha1hash(obj):
    m = hashlib.sha1()
    m.update(bytes(repr(obj), 'utf-8'))
    return m.hexdigest()


class DBschema(MSONable, metaclass=ABCMeta):

    @abstractmethod
    def getid(self) -> str:
        pass

    @abstractmethod
    def __hash__(self) -> str:
        pass

    @property
    def unfilled(self):
        d = self.as_dict()
        uf = []
        for k in d.keys():
            if d[k] is None:
                uf.append(k)
        return uf


from collections import Counter


class DBCrystal(DBschema):
    _id: str

    def __init__(
            self,

            lattice: Lattice = None,
            smiles: str = None,
            hsmiles: str = None,
            chromophore: str = None,

            elements: (str) = None,  # e.g. ('c', 'h', ...), must be sorted
            identifier: str = None,
            source: str = None,
            cif_string: str = None,
            disorder_location: str = None,  # 'bone' or 'sidechain' or 'not sure' or 'no disorder
            disorder_class: str = None,
            configurations: [str] = None,
            major_configuration: [str] = None,

            results: dict = None
    ):
        self.elements = elements
        self.identifier = identifier
        self.source = source
        self.cif_string = cif_string
        self.disorder_location = disorder_location
        self.disorder_class = disorder_class
        self.configurations = configurations
        self.major_configuration = major_configuration
        self.lattice = lattice
        self.smiles = smiles
        self.hsmiles = hsmiles
        self.chromophore = chromophore
        self.results = results

        self.id = self.getid()  # source + identifier  # e.g. 'csd_ASFNC'

    def getid(self):
        return "{}_{}".format(self.source, self.identifier)

    def __hash__(self):
        f = ','.join(['{:.4}'] * 9)
        c = Counter(self.elements)
        els = sorted(c.most_common(), key=lambda x: x[0])
        els = str(els)
        latmat_format = f.format(*np.round(self.lattice.matrix.flatten(), decimals=4))
        return sha1hash('@'.join([els, latmat_format]))


class DBConfiguration(DBschema):

    def __init__(
            self,
            calculation_stage: str = None,  # denotes how many calculations have been done for this config
            unwrap_structure: Structure = None,  # this is the unwrap_clean_pstructure
            molconformers: [str] = None,  # use rmsd to get rid of identical conformers, see BasicConformer.compare()
            z: int = None,
            occupancy: float = None,
            backbone_structure: Structure = None,
            crystal: str = None,
            config_index: int = None,
            is_major: bool = None,
            packing_data: dict = None,
            hop_data_geodict: dict = None,
            hop_data_zindo: dict = None,
            hop_data_dft: dict = None,
            fb_opt_structure: Structure = None,
            ar_opt_structure: Structure = None,
            fb_opt_ebs: dict = None,
            ar_opt_ebs: dict = None,
            energetics: dict = None,  # keys are fb_opt, ar_opt, mol_in_box, cohesive
    ):
        self.calculation_stage = calculation_stage
        self.unwrap_structure = unwrap_structure
        self.molconformers = molconformers
        self.z = z
        self.occupancy = occupancy
        self.backbone_structure = backbone_structure
        self.crystal = crystal
        self.config_index = config_index
        self.is_major = is_major
        self.packing_data = packing_data
        self.hop_data_geodict = hop_data_geodict
        self.hop_data_zindo = hop_data_zindo
        self.hop_data_dft = hop_data_dft
        self.fb_opt_structure = fb_opt_structure
        self.ar_opt_structure = ar_opt_structure
        self.fb_opt_ebs = fb_opt_ebs
        self.ar_opt_ebs = ar_opt_ebs
        self.energetics = energetics
        self.id = self.getid()

    def getid(self):
        return '_'.join([self.crystal, self.config_index])

    def __hash__(self):
        return self.id


class DBMolConformer(DBschema):
    def __init__(
            self,
            mc_index: int = None,
            molconformer: dict = None,
            Configuration: str = None,
            GeoBoneConf: dict = None,
            GeoBoneConf_descriptors: dict = None,
            GeoScsConf: dict = None,
            GeoScsConf_descriptors: dict = None,
            ChromBoneConf: dict = None,
            ChromBoneConf_descriptors: dict = None,
            ChromScsConf: dict = None,
            ChromScsConf_descriptors: dict = None,
    ):
        self.mc_index = mc_index
        self.molconformer = molconformer
        self.Configuration = Configuration
        self.GeoBoneConf = GeoBoneConf
        self.GeoBoneConf_descriptors = GeoBoneConf_descriptors
        self.GeoScsConf = GeoScsConf
        self.GeoScsConf_descriptors = GeoScsConf_descriptors
        self.ChromBoneConf = ChromBoneConf
        self.ChromBoneConf_descriptors = ChromBoneConf_descriptors
        self.ChromScsConf = ChromScsConf
        self.ChromScsConf_descriptors = ChromScsConf_descriptors
        self.id = self.getid()

    def getid(self):
        return '_'.join([self.Configuration, self.mc_index])

    def __hash__(self):
        return self.id


class DBChromophore(MSONable):
    def __init__(
            self,
            hsmiles: str = None,
            smiles: str = None,
            ChromophoreConfs: [str] = None
    ):
        self.hsmiles = hsmiles
        self.smiles = smiles
        self.ChromophoreConfs = ChromophoreConfs
        self.id = self.getid()

    def getid(self):
        return self.hsmiles

    def __hash__(self):
        return self.id


class DBChromophoreConformer(DBschema):
    def __init__(
            self,
            geo: dict = None,
            index: int = None,
            anion_geo: int = None,
            cation_geo: int = None,
            AIP: float = None,
            AEA: float = None,
            reorg: float = None,
            tddft: dict = None,
            vs0s1: float = None,
            vs0t1: float = None,
            as0t1: float = None,
            Chromophore: str = None
    ):
        self.geo = geo
        self.index = index
        self.anion_geo = anion_geo
        self.cation_geo = cation_geo
        self.AIP = AIP
        self.AEA = AEA
        self.reorg = reorg
        self.tddft = tddft
        self.vs0s1 = vs0s1
        self.vs0t1 = vs0t1
        self.as0t1 = as0t1
        self.Chromophore = Chromophore
        self.id = self.getid()

    def getid(self):
        return '_'.join([self.Chromophore, self.index])

    def __hash__(self):
        return self.id

# dc = DBCrystal()
# dc.lattice = Lattice.from_parameters(1, 2, 3, 40, 50, 60)
# from pprint import pprint
# pprint(dc.as_dict())