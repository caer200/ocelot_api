import itertools
import math
import re
import warnings
from copy import deepcopy

import networkx as nx
import numpy as np
from pymatgen.core.operations import SymmOp
from pymatgen.core.periodic_table import Element
from pymatgen.core.structure import PeriodicSite
from pymatgen.io.cif import CifFile
from pymatgen.io.cif import str2float
from pymatgen.io.cif import sub_spgrp
from pymatgen.symmetry.groups import SYMM_DATA
from pymatgen.symmetry.groups import SpaceGroup
from pymatgen.util.coord import pbc_shortest_vectors
from scipy.spatial.distance import cdist

from ocelot.routines.pbc import PBCparser


_COD_DATA = None


def _get_cod_data():
    global _COD_DATA
    if _COD_DATA is None:
        import pymatgen
        with open(os.path.join(pymatgen.symmetry.__path__[0],
                               "symm_ops.json")) \
                as f:
            import json
            _COD_DATA = json.load(f)

    return _COD_DATA

class CifFileError(Exception): pass

possible_symm_labels = [
    "_symmetry_equiv_pos_as_xyz",
    "_symmetry_equiv_pos_as_xyz_",
    "_space_group_symop_operation_xyz",
    "_space_group_symop_operation_xyz_",
    "_symmetry_space_group_name_H-M",
    "_symmetry_space_group_name_H_M",
    "_symmetry_space_group_name_H-M_",
    "_symmetry_space_group_name_H_M_",
    "_space_group_name_Hall",
    "_space_group_name_Hall_",
    "_space_group_name_H-M_alt",
    "_space_group_name_H-M_alt_",
    "_symmetry_space_group_name_hall",
    "_symmetry_space_group_name_hall_",
    "_symmetry_space_group_name_h-m",
    "_symmetry_space_group_name_h-m_",
    "_space_group_IT_number",
    "_space_group_IT_number_",
    "_symmetry_Int_Tables_number",
    "_symmetry_Int_Tables_number_"]

latt_labels = [
    '_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma',
]

chemistry_labels = [
    '_cell_formula_units_Z', '_chemical_formula_moiety', '_chemical_formula_sum'
]

coord_labels = [
    '_atom_site_label', '_atom_site_type_symbol', '_atom_site_fract_x', '_atom_site_fract_y',
    '_atom_site_fract_z'
]

disorder_labels = [
    '_atom_site_occupancy', '_atom_site_disorder_group',
]

space_groups = {sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()}
space_groups.update({sub_spgrp(k): k for k in SYMM_DATA['space_group_encoding'].keys()})


def get_pmg_dict(cifstring: str):
    """
    use pmg dict to parse cifstring, only deal with one structure per file

    :param cifstring:
    :return:
    """
    cifdata = CifFile.from_string(cifstring).data
    idnetifiers = list(cifdata.keys())
    if len(idnetifiers) > 1:
        warnings.warn('W: find more than 1 structures in this cif file!')
    elif len(idnetifiers) == 0:
        warnings.warn('W: no structure found by pymatgen parser!')
    try:
        identifier = idnetifiers[0]
    except IndexError:
        raise CifFileError('no identifier found in the ciffile!')
    pymatgen_dict = list(cifdata.items())[0][1].data

    # jmol writes '_atom_site_type_symbol', but not '_atom_site_label'
    if '_atom_site_label' not in pymatgen_dict.keys():
        warnings.warn('W: _atom_site_label not found in parsed dict')
        atom_site_label = []
        symbols = pymatgen_dict['_atom_site_type_symbol']
        for i in range(len(symbols)):
            s = symbols[i]
            atom_site_label.append('{}{}'.format(s, i))
        pymatgen_dict['_atom_site_label'] = atom_site_label
    return identifier, pymatgen_dict


def apply_symmop(psites, ops):
    """
    symmop and xyz in cif file:

    lets say xyz -- op1 --> x'y'z' and xyz -- op2 --> x!y!z! and
    it is possible to have x'y'z' is_close x!y!z!

    this means one should take only x'y'z' or x!y!z!, aka op1 is equivalent to op2 due to the symmetry implicated by
    xyz/the asymmectric unit, e.g. ALOVOO.cif -- Z=2, asymmectric unit given by cif is one molecule, but there're 4 ops

    so we need first check if the cif file behaves like this

    """
    op_xyzs = []
    for op in ops:
        n_xyzs = []
        for ps in psites:
            new_coord = op.operate(ps.frac_coords)
            # new_coord = np.array([i - math.floor(i) for i in new_coord])
            n_xyzs.append(new_coord)
        op_xyzs.append(n_xyzs)

    latt = psites[0].lattice

    def pbc_dist(fc1, fc2, lattice):
        v, d2 = pbc_shortest_vectors(lattice, fc1, fc2, return_d2=True)
        return math.sqrt(d2[0, 0])

    def pbc_distmat(fcl1, fcl2):
        distmat = np.zeros((len(fcl1), len(fcl1)))
        for i in range(len(fcl1)):
            for j in range(i, len(fcl1)):
                distmat[i][j] = pbc_dist(fcl1[i], fcl2[j], latt)
                distmat[j][i] = distmat[i][j]
        return distmat

    def two_xyzs_close(xyzs1, xyzs2, tol=1e-5):
        dmat = pbc_distmat(xyzs1, xyzs2)
        almost_zeros = dmat[(dmat < tol)]
        if len(almost_zeros) > 0:
            return True
        return False

    op_identities = np.zeros((len(ops), len(ops)), dtype=bool)
    for i, j in itertools.combinations(range(len(ops)), 2):
        ixyzs = op_xyzs[i]
        jxyzs = op_xyzs[j]
        if two_xyzs_close(ixyzs, jxyzs):
            op_identities[i][j] = True
            op_identities[j][i] = True

    groups = [[0]]
    for i in range(len(ops)):
        for ig in range(len(groups)):
            if all(op_identities[i][j] for j in groups[ig]):
                groups[ig].append(i)
        if i not in [item for sublist in groups for item in sublist]:
            groups.append([i])
    unique_ops = [ops[g[0]] for g in groups]

    new_psites = []
    for ps in psites:
        iasym = 0
        for op in unique_ops:
            new_coord = op.operate(ps.frac_coords)
            new_coord = np.array([i - math.floor(i) for i in new_coord])
            new_properties = deepcopy(ps.properties)
            new_properties['iasym'] = iasym
            new_ps = PeriodicSite(ps.species_string, new_coord, ps.lattice, properties=deepcopy(new_properties))
            new_psites.append(new_ps)
            iasym += 1

    return new_psites, unique_ops


def braket2float(s):
    try:
        return float(s)
    except ValueError:
        if isinstance(s, str):
            return str2float(s)
        raise TypeError('cannot parse {} into float'.format(s))


def get_symmop(data):
    symops = []
    for symmetry_label in ["_symmetry_equiv_pos_as_xyz",
                           "_symmetry_equiv_pos_as_xyz_",
                           "_space_group_symop_operation_xyz",
                           "_space_group_symop_operation_xyz_"]:
        if data.get(symmetry_label):
            xyz = data.get(symmetry_label)
            if isinstance(xyz, str):
                msg = "A 1-line symmetry op P1 CIF is detected!"
                warnings.warn(msg)
                xyz = [xyz]
            try:
                symops = [SymmOp.from_xyz_string(s)
                          for s in xyz]
                break
            except ValueError:
                continue
    if not symops:
        # Try to parse symbol
        for symmetry_label in ["_symmetry_space_group_name_H-M",
                               "_symmetry_space_group_name_H_M",
                               "_symmetry_space_group_name_H-M_",
                               "_symmetry_space_group_name_H_M_",
                               "_space_group_name_Hall",
                               "_space_group_name_Hall_",
                               "_space_group_name_H-M_alt",
                               "_space_group_name_H-M_alt_",
                               "_symmetry_space_group_name_hall",
                               "_symmetry_space_group_name_hall_",
                               "_symmetry_space_group_name_h-m",
                               "_symmetry_space_group_name_h-m_"]:
            sg = data.get(symmetry_label)

            if sg:
                sg = sub_spgrp(sg)
                try:
                    spg = space_groups.get(sg)
                    if spg:
                        symops = SpaceGroup(spg).symmetry_ops
                        msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
                              "Spacegroup from %s used." % symmetry_label
                        warnings.warn(msg)
                        break
                except ValueError:
                    # Ignore any errors
                    pass

                try:
                    for d in _get_cod_data():
                        if sg == re.sub(r"\s+", "",
                                        d["hermann_mauguin"]):
                            xyz = d["symops"]
                            symops = [SymmOp.from_xyz_string(s)
                                      for s in xyz]
                            msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
                                  "Spacegroup from %s used." % symmetry_label
                            warnings.warn(msg)
                            break
                except Exception:
                    continue

                if symops:
                    break
    if not symops:
        # Try to parse International number
        for symmetry_label in ["_space_group_IT_number",
                               "_space_group_IT_number_",
                               "_symmetry_Int_Tables_number",
                               "_symmetry_Int_Tables_number_"]:
            if data.get(symmetry_label):
                try:
                    i = int(braket2float(data.get(symmetry_label)))
                    symops = SpaceGroup.from_int_number(i).symmetry_ops
                    break
                except ValueError:
                    continue

    if not symops:
        msg = "No _symmetry_equiv_pos_as_xyz type key found. " \
              "Defaulting to P1."
        warnings.warn(msg)
        symops = [SymmOp.from_xyz_string(s) for s in ['x', 'y', 'z']]

    return symops


def get_connected_pblock(pblocks, cutoff=4.0):
    def get_coords(pb: [PeriodicSite]):
        coords = np.zeros((len(pb), 3))
        for i in range(len(pb)):
            coords[i] = pb[i].coords
        return coords

    def get_shortest_distance_between_blocks(pb1, pb2):
        pb1coords = get_coords(pb1)
        pb2coords = get_coords(pb2)
        distmat = cdist(pb1coords, pb2coords)
        minid = np.unravel_index(np.argmin(distmat, axis=None), distmat.shape)
        return distmat[minid]

    def get_block_graph(distmat, cutoff):
        g = nx.Graph()
        for i in range(len(distmat)):
            g.add_node(i)
            for j in range(i + 1, len(distmat)):
                g.add_node(j)
                if distmat[i][j] < cutoff:
                    g.add_edge(i, j)
        return g

    distmat_pblocks = np.zeros((len(pblocks), len(pblocks)))
    for i in range(len(pblocks)):
        for j in range(i + 1, len(pblocks)):
            distmat_pblocks[i][j] = get_shortest_distance_between_blocks(pblocks[i], pblocks[j])
            distmat_pblocks[j][i] = distmat_pblocks[i][j]
    block_graph = get_block_graph(distmat_pblocks, cutoff)
    connected_block_ids = nx.connected_components(block_graph)  # [[1,3], [2, 4, 5], ...]
    merged_blocks = []
    for ids in connected_block_ids:
        merged_block = []
        for i in ids:
            merged_block += pblocks[i]
        merged_blocks.append(merged_block)
    return merged_blocks


def find_connected_psites(psites: [PeriodicSite]):
    """
    simplified version of percolation in unwrap
    """
    latt = psites[0].lattice
    pindices = range(len(psites))
    visited = []
    block_list = []
    while len(visited) != len(psites):
        # initialization
        unvisited = [idx for idx in pindices if idx not in visited]
        ini_idx = unvisited[0]
        block = [ini_idx]
        pointer = 0
        while pointer != len(block):
            outside = [idx for idx in pindices if idx not in block and idx not in visited]
            for i in outside:
                distance, fctrans = PBCparser.get_dist_and_trans(latt,
                                                                 psites[block[pointer]]._frac_coords,
                                                                 psites[i]._frac_coords, )

                cutoff = Element(psites[block[pointer]].species_string).atomic_radius + Element(
                    psites[i].species_string).atomic_radius
                cutoff *= 1.3
                if distance < cutoff:
                    block.append(i)
            visited.append(block[pointer])
            pointer += 1
        block_list.append(block)
    return [[psites[j] for j in b] for b in block_list]
