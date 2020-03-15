import numpy as np
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Descriptors3D as D3d
from rdkit.Chem.Fingerprints import FingerprintMols

"""
commonly used rdkit functions
"""

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
from rdkit.Chem import Mol
from rdkit.Chem import ResonanceMolSupplier


class RdFunc:

    @staticmethod
    def get_volume(m):
        """
        AllChem.ComputeMolVolume in rdkit

        :return: volume in A^3
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
            Asphericity=D3d.Asphericity(mol),
            Eccentricity=D3d.Eccentricity(mol),
            InertialShapeFactor=D3d.InertialShapeFactor(mol),
            NPR1=D3d.NPR1(mol),
            NPR2=D3d.NPR2(mol),
            PMI1=D3d.PMI1(mol),
            PMI2=D3d.PMI2(mol),
            PMI3=D3d.PMI3(mol),
            RadiusOfGyration=D3d.RadiusOfGyration(mol),
            SpherocityIndex=D3d.SpherocityIndex(mol),
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
    def from_string(fmt, string):
        """
        construct a rdmol

        :param fmt: 'smiles', 'hsmiles', 'smarts', 'inchi'
        :param string:
        :return:
        """
        if fmt == 'smiles':
            rdmol = Chem.MolFromSmiles(string)
            rdmol = Chem.AddHs(rdmol)
        elif fmt == 'hsmiles':
            rdmol = Chem.MolFromSmiles(string)
            rdmol = Chem.AddHs(rdmol, explicitOnly=True)
        elif fmt == 'smarts':
            rdmol = Chem.MolFromSmarts(string)
        elif fmt == 'inchi':
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

        :param m1:
        :param m2:
        :param str metric: "Tanimoto", "Dice", "Cosine", "Sokal", "Russel", "RogotGoldberg", "AllBit", "Kulczynski", "McConnaughey", "Asymmetric", "BraunBlanquet",
        :return:
        """
        fps = [FingerprintMols.FingerprintMol(m1), FingerprintMols.FingerprintMol(m2)]
        for func in DataStructs.similarityFunctions:
            if func[0] == metric:
                metric_function = func[1]
                return DataStructs.FingerprintSimilarity(fps[0], fps[1], metric=metric_function)
        return None

    @staticmethod
    def draw_smiles(smiles, fn='m.ps'):
        rdmol = RdFunc.from_string('smiles', smiles)
        Draw.MolToFile(rdmol, fn)

    @staticmethod
    def draw_mol(rdmol, fn='m.ps'):
        Draw.MolToFile(rdmol, fn)

    @staticmethod
    def get_sub_rdmol(m: Mol, atomids: [int]):
        atoms_in_old_mol: [Atom] = [a for a in m.GetAtoms() if a.GetIdx() in atomids]
        atom_numbers = [a.GetAtomicNum() for a in atoms_in_old_mol]

        old_id_2_new_id = {}
        newid = 0
        for oldatom in atoms_in_old_mol:
            old_id = oldatom.GetIdx()
            old_id_2_new_id[old_id] = newid
            newid += 1

        mol = Chem.MolFromSmarts("[#" + str(atom_numbers[0]) + "]")
        rwmol = Chem.RWMol(mol)
        for s in atom_numbers[1:]:
            rwmol.AddAtom(Chem.Atom(s))

        # print('new mol atom')
        # for a in rwmol.GetAtoms():
        #     print(a.GetIdx(), a.GetSymbol())
        # print('--')

        for aini, ainj in combinations(atomids, 2):
            b = m.GetBondBetweenAtoms(aini, ainj)
            if isinstance(b, Bond):
                # iatom = m.GetAtomWithIdx(aini).GetSymbol()
                # jatom = m.GetAtomWithIdx(ainj).GetSymbol()
                # print('found bond {} {} - {} {}, {}'.format(iatom, aini, jatom, ainj, b.GetBondType()))
                bt = b.GetBondType()
                newi = old_id_2_new_id[aini]
                newj = old_id_2_new_id[ainj]
                rwmol.AddBond(newi, newj, bt)
                # newatomi = rwmol.GetAtomWithIdx(newi).GetSymbol()
                # newatomj = rwmol.GetAtomWithIdx(newj).GetSymbol()
                # print('added {} {} - {} {}'.format(newatomi, newi, newatomj, newj))
        mol = rwmol.GetMol()
        return mol, old_id_2_new_id

    @staticmethod
    def get_conjugate_group(m: Mol):
        # supp = ResonanceMolSupplier(m, Chem.KEKULE_ALL)
        supp = ResonanceMolSupplier(m, Chem.ALLOW_CHARGE_SEPARATION)
        cg_dict = {}
        a: Atom
        for a in m.GetAtoms():
            aid = a.GetIdx()
            cgid = supp.GetAtomConjGrpIdx(aid)
            if cgid < 1e5:
                cg_dict[aid] = cgid
        cgids = set(cg_dict.values())
        cgs = []
        for cgid in cgids:
            cg = [i for i in cg_dict.keys() if cg_dict[i] == cgid]
            cgmol, old_id_2_new_id = RdFunc.get_sub_rdmol(m, cg)
            cgs.append([cgmol, old_id_2_new_id])
        return sorted(cgs, key=lambda x: x[0].GetNumAtoms(), reverse=True)

    @staticmethod
    def get_conjugate_group_with_halogen(m: Mol):
        natoms = len(m.GetAtoms())
        adjmat = np.zeros((natoms, natoms), dtype=bool)
        for i in range(natoms):
            for j in range(i + 1, natoms):
                if isinstance(m.GetBondBetweenAtoms(i, j), Bond):
                    adjmat[i][j] = True
                    adjmat[j][i] = True

        supp = ResonanceMolSupplier(m, )
        # supp = ResonanceMolSupplier(m, Chem.KEKULE_ALL)
        # supp = ResonanceMolSupplier(m, Chem.ALLOW_CHARGE_SEPARATION)
        cg_dict = {}
        a: Atom
        for a in m.GetAtoms():
            aid = a.GetIdx()
            cgid = supp.GetAtomConjGrpIdx(aid)
            if cgid < 1e5:
                cg_dict[aid] = cgid
        cgids = set(cg_dict.values())
        cgs = []
        for cgid in cgids:
            cg = [i for i in cg_dict.keys() if cg_dict[i] == cgid]
            atom: Atom
            for atom in m.GetAtoms():
                if atom.GetIdx() not in cg:
                    if any(adjmat[atom.GetIdx()][cg_aid] for cg_aid in cg) and atom.GetSymbol() in (
                    "I", "F", "Cl", "Br"):
                        cg.append(atom.GetIdx())
            cgmol, old_id_2_new_id = RdFunc.get_sub_rdmol(m, cg)
            cgs.append([cgmol, old_id_2_new_id])
        return sorted(cgs, key=lambda x: x[0].GetNumAtoms(), reverse=True)

    @staticmethod
    def vis_partition(m: Mol, fnprefix='m'):
        cgs = RdFunc.get_conjugate_group_with_halogen(m)
        # cgs = RdFunc.get_conjugate_group(m)
        Draw.MolToFile(m, "{}.ps".format(fnprefix), kekulize=False)
        kk = 0
        for cg in cgs:
            Chem.SanitizeMol(cg[0])
            Draw.MolToFile(cg[0], "{}-{}.ps".format(fnprefix, kk), kekulize=False, )
            kk += 1

    @staticmethod
    def get_aids(m: Mol):
        aids = []
        for a in m.GetAtoms():
            aid = a.GetIdx()
            aids.append(aid)
        return aids

    @staticmethod
    def mol2xyz_by_confid(molecule: Mol, prefix='rdmol', confid=0, comment_line=''):
        natoms = molecule.GetNumAtoms()
        filename = "{}_{}.xyz".format(prefix, confid)
        s = "{}\n{}\n".format(natoms, comment_line)
        for i in range(natoms):
            position = molecule.GetConformer(confid).GetAtomPosition(i)
            symbol = molecule.GetAtomWithIdx(i).GetSymbol()
            s += "{}\t{:.6} {:.6} {:.6}\n".format(symbol, position.x, position.y, position.z)
        with open(filename, 'w') as f:
            f.write(s)

    @staticmethod
    def conf2xyz(conf: Chem.Conformer, outputname, atom_list:list, comment_line=''):
        natoms = conf.GetNumAtoms()
        s = "{}\n{}\n".format(natoms, comment_line)
        for i in range(natoms):
            position = conf.GetAtomPosition(i)
            symbol = atom_list[i]
            s += "{}\t{:.6} {:.6} {:.6}\n".format(symbol, position.x, position.y, position.z)
        with open(outputname, 'w') as f:
            f.write(s)
