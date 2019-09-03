from api.schema.omol import OMol
import numpy as np
import itertools
from copy import deepcopy
from api.routines.pbc import PBCparser
from pymatgen.core.structure import Structure, PeriodicSite
from api.schema.dimer import Dimer


class Config:

    def __init__(self, pstructure):
        self.pstructure = pstructure
        self.mols, self.omols, self.unwrap_structure = PBCparser.squeeze(self.pstructure)
        self.z = len(self.omols)  # no solvent!
        for i in range(self.z):
            self.omols[i].comment = {'index in the cell': i}
        self.dimers, self.transv_fcs = self.get_dimers_array(1)

    def as_dict(self):
        d = {"@module": self.__class__.__module__,
             "@class": self.__class__.__name__}
        d['pymatgen_structure'] = self.pstructure.as_dict()
        d['mols'] = [m.as_dict() for m in self.mols]
        d['omols'] = [m.as_dict() for m in self.omols]
        d['z'] = self.z

        dimers_dictarray = np.empty((self.z, self.z, len(self.transv_fcs)), dtype=dict)
        for i in range(self.z):
            for j in range(self.z):
                for k in range(len(self.transv_fcs)):
                    dimers_dictarray[i][j][k] = self.dimers[i][j][k].as_dict()
        d['dimers_dict_array'] = dimers_dictarray
        return d

    def get_dimers_array(self, maxfold=2):
        transv_1d = range(-maxfold, maxfold + 1)
        transv_fcs = list(v for v in itertools.product(transv_1d, transv_1d, transv_1d))
        dimers = np.empty((self.z, self.z, len(transv_fcs)), dtype=Dimer)

        for i in range(self.z):
            ref_omol = self.omols[i]
            for j in range(self.z):
                var_omol = self.omols[j]
                for k in range(len(transv_fcs)):
                    transv = self.unwrap_structure.lattice.get_cartesian_coords(transv_fcs[k])
                    msites = deepcopy(var_omol.sites)
                    for h in range(len(msites)):
                        msites[h].coords += transv
                    dimer_ijk = Dimer(deepcopy(ref_omol), OMol(msites), self, label="{}-{}_{}".format(i, j, k))
                    dimers[i][j][k] = dimer_ijk
        return dimers, transv_fcs

    def get_bone_config(self):

        terminated_backbones = [m.backbone.terminate() for m in self.omols]
        backbone_sites = []
        for b in terminated_backbones:
            backbone_sites += b.sites

        boneonly_psites = [PeriodicSite(s.element.name, s.coords, self.unwrap_structure.lattice, to_unit_cell=True, coords_are_cartesian=True)
                           for s in backbone_sites]
        boneonly_pstructure = Structure.from_sites(boneonly_psites)
        return Config(boneonly_pstructure), boneonly_pstructure, terminated_backbones



    @classmethod
    def from_cifstring(cls, string):
        s = Structure.from_str(string, fmt='cif')
        return cls(s)

    @classmethod
    def from_file(cls, filename):
        s = Structure.from_file(filename)
        return cls(s)
