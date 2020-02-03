from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors3D as descriptors3d

"""
commonly used rdkit functions
"""


class RdFunc:

    @staticmethod
    def get_volume(m):
        """
        AllChem.ComputeMolVolume in rdkit

        :return: volume in \AA^3
        """
        return AllChem.ComputeMolVolume(m, confId=-1, gridSpacing=0.2, boxMargin=2.0)

    @staticmethod
    def get_3ddescriptors(mol):
        """
        using methods from rdkit

        https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors3D.html

        Asphericity: 0 corresponds to spherical top molecules and 1 to linear molecules, For prolate (cigar-sh aped) molecules, ~ 0.25, whereas for oblate (disk-shaped) molecules, ~ 1.

        Eccentricity: 0 corresponds to spherical top molecules and 1 to linear molecules.

        RadiusOfGyration: a size descriptor based on the distribution of atomic masses in a molecule, a measure of molecular compactness for long-chain molecules and, specifically, small values are obtained when most o f the atoms are close to the center of mass

        SpherocityIndex: spherosity index varies from zero for flat molecules, such as benzene, to unity for totally spherical molecules

        :key: Asphericity, Eccentricity, InertialShapeFactor, NPR1, NPR2, PMI1, PMI2, PMI3, RadiusOfGyration, SpherocityIndex

        :return: dict
        """
        d = dict(
            Asphericity=descriptors3d.Asphericity(mol),
            Eccentricity=descriptors3d.Eccentricity(mol),
            InertialShapeFactor=descriptors3d.InertialShapeFactor(mol),
            NPR1=descriptors3d.NPR1(mol),
            NPR2=descriptors3d.NPR2(mol),
            PMI1=descriptors3d.PMI1(mol),
            PMI2=descriptors3d.PMI2(mol),
            PMI3=descriptors3d.PMI3(mol),
            RadiusOfGyration=descriptors3d.RadiusOfGyration(mol),
            SpherocityIndex=descriptors3d.SpherocityIndex(mol),
        )
        return d

    @staticmethod
    def get_bond_order(m, atomidx1, atomidx2):
        bt = m.GetBondBetweenAtoms(atomidx1, atomidx2).GetBondType()
        if bt == Chem.BondType.SINGLE:
            return 1
        elif bt == Chem.BondType.DOUBLE:
            return 2
        elif bt == Chem.BondType.TRIPLE:
            return 3
        elif bt == Chem.BondType.AROMATIC:
            return 1.5
        else:
            return None

    @staticmethod
    def from_string(format, string):
        """
        construct a rdmol

        :param format: 'smiles', 'hsmiles', 'smarts', 'inchi'
        :param string:
        :return:
        """
        if format == 'smiles':
            rdmol = Chem.MolFromSmiles(string)
            rdmol = Chem.AddHs(rdmol)
        elif format == 'hsmiles':
            rdmol = Chem.MolFromSmiles(string)
            rdmol = Chem.AddHs(rdmol, explicitOnly=True)
        elif format == 'smarts':
            rdmol = Chem.MolFromSmarts(string)
        elif format == 'inchi':
            rdmol = Chem.MolFromInchi(string)
        else:
            rdmol = None
        return rdmol

    @staticmethod
    def substructurematch(m1, m2):
        """
        substructure search with smarts, True if m2 has sub graph of m1

        :param m1: pattern used to search
        :param m2:
        :return:
        """
        patt = m1.to_smarts()
        return m2.HasSubstructMatch(patt)

    @staticmethod
    def getmolweight(mol):
        return Descriptors.MolWt(mol)

    @staticmethod
    def to_smarts(rdmol):
        return Chem.MolToSmarts(rdmol)

    @staticmethod
    def n_nonsignlebond_electrons(m):
        nve = Descriptors.NumValenceElectrons(m)
        nre = Descriptors.NumRadicalElectrons(m)
        nbonds = len(m.GetBonds())
        return nve - nre - 2 * nbonds

    @staticmethod
    def unsaturation(m):
        return RdFunc.n_nonsignlebond_electrons(m) / len(m.GetAtoms())

    @staticmethod
    def fingerprint(mol):
        return FingerprintMols.FingerprintMol(mol)

    @staticmethod
    def fp_similarity(m1, m2, metric='Tanimoto'):
        """
        use RDK fingerprint similarity based on different metrics

        TODO add args to customize RDKfp, see https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RDKFingerprint

        see Landrum2012 for more details

        :param str metric: "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet",
        :return:
        """
        fps = [FingerprintMols.FingerprintMol(m1), FingerprintMols.FingerprintMol(m2)]
        for func in DataStructs.similarityFunctions:
            if func[0] == metric:
                metric_function = func[1]
                return DataStructs.FingerprintSimilarity(fps[0], fps[1], metric=metric_function)
        return None
