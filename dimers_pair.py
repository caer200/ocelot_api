from ocelot.task.zindo import ZindoJob, conver2zindo
from pymatgen.core.structure import Molecule
import os
from scipy.spatial.distance import cdist
import numpy as np

ZINDOLIB = '/home/vinayak/ZINDO'
ZINDOBIN = '/home/vinayak/ZINDO/zindo'
ZINDOCTBIN = '/home/vinayak/ZINDO/zindo-ct'


def get_hh_coupling(coupling_data):
    """

    :param coupling_data: the return value of Zindo.dimer_run
    :return: the HOMO-HOMO electronic coupling in meV
    """
    nmo_a = coupling_data[1]
    nmo_b = coupling_data[1]
    hh_coupling = coupling_data[0][int(nmo_a / 2 - 1)][int(nmo_b / 2 - 1)]['ti']
    return hh_coupling


def get_close_molecules(dimer_array, mols_in_cell):
    """

    :param dimer_array: Array of dimers and their translation vectors
    :param mols_in_cell: Number of molecules in the unit cell
    :return: A cluster of dimers less than 5.5 Anstrong apart. Each element
            of the list is of the form [Dimer, trans_vector]
    """
    for mol in range(mols_in_cell):
        dimer_close_mol = []
        for dimer in range(dimer_array[0].size):
            if dimer_array[0][mol][mol][dimer].is_close:
                if dimer_array[0][mol][mol][dimer].is_not_identical:
                    dimer_close_mol.append([dimer_array[0][mol][mol][dimer], dimer_array[1][dimer]])
    return dimer_close_mol


def get_vectors(dimer_close):
    """

    :param dimer_close: A list of list of dimer and its translation vector i.e. [Dimer, trans_v]
    :return: list of translation vectors from the dimer_close list
    """
    vectors = []
    for index in range(len(dimer_close)):
        vectors.append(dimer_close[index][1])
    return vectors


def get_unique_dimers(dimer_array, mol):
    """
    The unique dimer pairs can be extracted from all possible dimer pairs
    :param dimer_array: The dimer array object
    :param mol: number of molecules in the unit cell
    :return: list of list of dimer that are unique. The element of the list is of the form [Dimer, trans_v]
    """
    dimer_close_unique = get_close_molecules(dimer_array, mol)
    for index_1 in range(len(dimer_close_unique)):
        for index_2 in range(index_1 + 1, len(dimer_close_unique)):
            if abs(dimer_close_unique[index_1][0].oslip - dimer_close_unique[index_2][0].oslip) < 10e-5:
                if abs(dimer_close_unique[index_1][0].qslip - dimer_close_unique[index_2][0].qslip) < 10e-5:
                    if abs(dimer_close_unique[index_1][0].pslip - dimer_close_unique[index_2][0].pslip) < 10e-5:
                        del dimer_close_unique[index_2]
    return dimer_close_unique


def write_dimer_close(dimer_close):
    """
    Writes a xyz file for the cluster of dimers in a folder named 'structures'

    :param dimer_close: a list with each element of the form [Dimer, trans_v]
    :return: None
    """
    os.system('mkdir structures')
    os.chdir('structures')
    sites = []
    for i in range(len(dimer_close)):
        dimer_close[i][0].to_xyz(dimer_close[i][0].label + '.xyz')
        sites += [s.to_pymatgen_site() for s in dimer_close[i][0].omol_var.msites]
    sites += [s.to_pymatgen_site() for s in dimer_close[0][0].omol_ref.msites]
    Molecule.from_sites(sites).to('xyz', "dimer_close_collection.xyz")


def run_zindo(dimer_close, workdir):
    """
    Run zindo for all the dimers provided and calculate the HOMO-HOMO electronic coupling.
    xzy files are written for each of the dimer in the respective folder. xyz file for the dimer cluster is also written.

    :param dimer_close: A list with elements of form [Dimer, trans_v]
    :param workdir: working directory to run zindo
    :return: a list of coupling values for the dimers
    """
    hh_coupling = []
    sites = []
    os.system('mkdir ' + workdir + 'dimers')
    for dimer in range(len(dimer_close)):
        mol_A = conver2zindo(dimer_close[dimer][0].omol_ref.to_pymatgen_mol())
        mol_D = conver2zindo(dimer_close[dimer][0].omol_var.to_pymatgen_mol())
        os.system('mkdir ' + workdir + 'dimers/' + dimer_close[dimer][0].label)
        os.chdir(workdir + 'dimers/' + dimer_close[dimer][0].label)
        coupling_data = ZindoJob.dimer_run(dimer_close[dimer][0].label, workdir + 'dimers/' + dimer_close[dimer][0].label,
                                  ZINDOBIN,
                                  ZINDOCTBIN, ZINDOLIB, mol_A, mol_D)
        dimer_close[dimer][0].to_xyz(dimer_close[dimer][0].label + '.xyz')
        sites += [s.to_pymatgen_site() for s in dimer_close[dimer][0].omol_var.msites]
        os.chdir(workdir + 'dimers')
        hh_coupling.append(get_hh_coupling(coupling_data))
    sites += [s.to_pymatgen_site() for s in dimer_close[0][0].omol_ref.msites]
    Molecule.from_sites(sites).to('xyz', "dimercollection.xyz")
    return hh_coupling


def write_hopping_mols(dimer_array, mol_id):
    """
    Writes a xyz file with combined coordinates from the provided Dimer.label. Works for one molecule in unit cell
    :param dimer_array: Dimer array object
    :param mol_id: Dimer.label
    :return: None
    """
    sites = []
    for mol in mol_id:
        sites += [s.to_pymatgen_site for s in dimer_array[0][0][0][int(mol)].omol_var.msites]
    sites += [s.to_pymatgen_site for s in dimer_array[0][0][0][0].omol_ref.msites]
    Molecule.from_sites(sites).to('xyz','path.xyz')
    return

def is_molecule_close(ref_coord, var_coord, cutoff=5.5):
    """
    Check whether molecules are closer the cutoff

    :param ref_coord: Dimer.omol_ref.coord
    :param var_coord:  Dimer.omol_var.coord
    :param cutoff: in Anstrong
    :return: Bool
    """
    distmat = cdist(ref_coord, var_coord)
    return np.min(distmat) < cutoff



def sym_mol(dimer_array, vector):
    """

    :param dimer_array: Dimer array object
    :param vector: trans_v
    :return: inverted trans_v. For [1,1,0], [-1,-1,0] is returned
    """
    for i in range(len(dimer_array[1])):
        if np.all(dimer_array[1][i] == -1 * vector):
            return i


def hopping(dimer_array, hop_vector):
    """
    The dimer.omol_ref is considered as the initial position. The dimer.omol_var that are close are considered. The hop
    is in the direction of the vector. The molecule hopped to is the considered as the initial position for next hop.

    :param dimer_array: Dimer array object
    :param hop_vector: trans_v for hopping
    :return: COM of the molecules that are connected by hops and their labels
    """
    com_hop = []
    mol_id = []
    ref_v = [0, 0, 0]
    ref_coord = dimer_array[0][0][0][0].omol_ref.backbone.coordmat
    ref_id = dimer_array[0][0][0][0].label.split('-')[-1]
    ref_com = dimer_array[0][0][0][0].omol_ref.geoc
    hop = 0

    while hop < 10:
        if ref_id not in mol_id:
            mol_id.append(ref_id)
            com_hop.append(ref_com)
            if not np.all(ref_v == [0, 0, 0]):
                mol_id.append(dimer_array[0][0][0][sym_mol(dimer_array, ref_v)].label.split('-')[-1])
                com_hop.append(dimer_array[0][0][0][sym_mol(dimer_array, ref_v)].omol_var.geoc)
            var_list = []
            for i in range(len(dimer_array[1])):
                if dimer_array[0][0][0][i].is_not_identical:
                    var_coord = dimer_array[0][0][0][i].omol_var.backbone.coordmat
                    var_id = dimer_array[0][0][0][i].label.split('-')[-1]
                if var_id not in mol_id:
                    if is_molecule_close(ref_coord, var_coord):
                        var_com = dimer_array[0][0][0][i].omol_var.geoc
                        var_v = dimer_array[1][i]
                        var_list.append([var_v, var_id, var_coord, var_com])

            for molecule in var_list:
                if np.all(molecule[0] - ref_v == hop_vector):
                    ref_v = molecule[0]
                    ref_id = molecule[1]
                    ref_coord = molecule[2]
                    ref_com = molecule[3]

        hop += 1

    return com_hop, mol_id

