from pymatgen.core.sites import Site
from pymatgen.io.xyz import XYZ
from pymatgen.core.structure import Molecule


class DimerCollection:
    """
    they share the same ref_mol
    """

    def __init__(self, dimers):
        self.dimers = dimers

    def get_xyz_string(self):
        sites = []
        for d in self.dimers:
            sites += [s.to_pymatgen_site() for s in d.omol_var.msites]
        sites += [Site('La', ss.coords) for ss in self.dimers[0].omol_ref.msites]

        mol = Molecule.from_sites(sites)
        xyz = XYZ(mol)
        return str(xyz)

    def to_xyz(self, fn):
        sites = []
        for d in self.dimers:
            sites += [s.to_pymatgen_site() for s in d.omol_var.msites]
        sites += [Site('La', ss.coords) for ss in self.dimers[0].omol_ref.msites]
        Molecule.from_sites(sites).to('xyz', fn)
