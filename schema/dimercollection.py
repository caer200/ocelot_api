from pymatgen.core.sites import Site
from pymatgen.io.xyz import XYZ
from pymatgen.core.structure import Molecule


class DimerCollection:
    def __init__(self, dimers):
        """
        just a list of dimers, they should share the same ref_mol
        """
        self.dimers = dimers

    def get_xyz_string(self):
        """
        put a La atom at the center of omol_ref, used for visualization

        :return: xyz string to be written
        """
        sites = []
        for d in self.dimers:
            sites += [s.to_pymatgen_site() for s in d.omol_var.msites]
        sites += [Site('La', ss.coords) for ss in self.dimers[0].omol_ref.msites]

        mol = Molecule.from_sites(sites)
        xyz = XYZ(mol)
        return str(xyz)

    def to_xyz(self, fn):
        """
        put a La atom at the center of omol_ref, used for visualization

        :param fn: xyz file name, with extension
        """
        sites = []
        for d in self.dimers:
            sites += [s.to_pymatgen_site() for s in d.omol_var.msites]
        sites += [Site('La', ss.coords) for ss in self.dimers[0].omol_ref.msites]
        Molecule.from_sites(sites).to('xyz', fn)
