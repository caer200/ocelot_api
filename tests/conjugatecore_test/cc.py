"""
this is used to extract chromophore that would be used for molecular calculations
"""

TEST_SMILES = [
    "c1ccc(c2ccccc2)cc1",
    "Oc1cccc(c2ccccc2)c1",
    "O=C1C=CC(=O)C(=C1)c1ccccc1",
    "O=C1C=C(c2ccccc2)C(=O)c2sccc12",
    "C1=CC(c2ccccc2)C=C1",
    "Clc1ccc2ccccc2c1",
    "C=C/C(=C\C)/CCc1cc2ccccc2cc1Cl",
    "C=Nc1cc(OCC)cc2cc(C#N)c(F)cc12"
]

# https://sourceforge.net/p/rdkit/mailman/message/34580312/

from itertools import combinations

from rdkit import Chem
from rdkit.Chem import Atom
from rdkit.Chem import Bond
from rdkit.Chem import Draw
from rdkit.Chem import Mol
from rdkit.Chem import ResonanceMolSupplier

from ocelot.schema.rdfunc import RdFunc


def draw_smiles(smiles, fn='m.ps'):
    rdmol = RdFunc.from_string('smiles', smiles)
    Draw.MolToFile(rdmol, fn)


# draw_smiles(TEST_SMILES[0])

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
    return mol


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
        cgmol = get_sub_rdmol(m, cg)
        cgs.append(cgmol)
    return sorted(cgs, key=lambda x: x.GetNumAtoms(), reverse=True)


import numpy as np


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
                if any(adjmat[atom.GetIdx()][cg_aid] for cg_aid in cg) and atom.GetSymbol() in ("I", "F", "Cl", "Br"):
                    cg.append(atom.GetIdx())
        cgmol = get_sub_rdmol(m, cg)
        cgs.append(cgmol)
    return sorted(cgs, key=lambda x: x.GetNumAtoms(), reverse=True)


jj = 0
for smiles in TEST_SMILES:
    m = RdFunc.from_string('smiles', smiles)
    # cgs = get_conjugate_group(m)
    cgs = get_conjugate_group_with_halogen(m)
    Draw.MolToFile(m, "{}.ps".format(jj), kekulize=False)
    kk = 0
    cg: Mol
    for cg in cgs:
        Chem.SanitizeMol(cg)
        Draw.MolToFile(cg, "{}-{}.ps".format(jj, kk), kekulize=False, )
        kk += 1
    jj += 1
