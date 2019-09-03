import math
import itertools
import warnings
import sys
import numpy as np
import mathop
from pymatgen.util.coord import pbc_shortest_vectors
from pymatgen.core.structure import Molecule, Site, PeriodicSite, IMolecule, IStructure
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius


def compare_mol(mol1, mol2):
    # TODO compare with openbabel
    pass


def dimer_massage(mol, xsliprange, ysliprange, zsliprange, prange, hrange):
    # TODO return a list of dimers as imol objs,
    pass


def slip_ribbon_mol(mol, xslip, yslip, zslip, p_angle, h_angle):
    # TODO return a copy with translation & rotation
    pass


def auto_terminate(subset_mol, wholemol):
    """
    terminate the subset with hydrogens
    :param subset_mol
    :param wholemol
    :return:
    """
    ch_bond_length = 1.09
    term_sites = []
    for site in subset_mol.sites:
        nbs = get_nbs_in_mol(site, subset_mol)
        if identify_cnhybrid(site, wholemol) == 'csp2':
            if len(nbs) == 2:
                v1 = site.coords - nbs[0].coords
                v2 = site.coords - nbs[1].coords
                term_site = Site('H', mathop.unify(v1 + v2) * ch_bond_length + site.coords)
                term_sites.append(term_site)
            if len(nbs) == 1:
                # TODO add two hydrogens with proper angle
                pass
        elif identify_cnhybrid(site, wholemol) == 'csp':
            if len(nbs) == 1:
                # TODO add 3 hydrogens or 1 hydrogen denpend on what is the site's nb
                pass
    return Molecule.from_sites(subset_mol.sites + term_sites)


def plane_project(p1, p2, p3, mol):
    """
    shift a mol with site coords on a plane
    :param p1:
    :param p2:
    :param p3:
    :param mol:
    :return:
    """
    if mathop.iscolinear(p1, p2, p3):
        warnings.warn('plane project failed as 3 pts on a line')
        sys.exit(1)
    proj_sites = []
    for s in mol.sites:
        new_coords = mathop.kart_trans(p1, p2, p3, s.coords)
        new_s = Site(s.species_string, new_coords)
        proj_sites.append(new_s)
    return Molecule.from_sites(proj_sites), (p1, p2, p3)


def plane_project_self(mol):
    """
    shift a mol with site coords on a plane
    :param mol:
    :return:
    """
    p1 = np.zeros(3)
    p2 = np.zeros(3)
    p3 = np.zeros(3)
    for i in itertools.combinations(range(3), 3):
        p1 = mol.sites[i[0]].coords
        p2 = mol.sites[i[1]].coords
        p3 = mol.sites[i[2]].coords
        if not mathop.iscolinear(p1, p2, p3):
            break
    proj_sites = []
    for s in mol.sites:
        new_coords = mathop.kart_trans(p1, p2, p3, s.coords)
        new_s = Site(s.species_string, new_coords)
        proj_sites.append(new_s)
    return Molecule.from_sites(proj_sites), (p1, p2, p3)


def get_dist_and_trans(fc1, fc2, lat):
    """
    get the shortest distance and corresponding translation vector between two frac coords
    :param fc1:
    :param fc2:
    :param lat:
    :return:
    """
    v, d2 = pbc_shortest_vectors(lat, fc1, fc2, return_d2=True)
    fc = lat.get_fractional_coords(v[0][0]) + fc1 - fc2
    return math.sqrt(d2[0, 0]), fc


def get_cutoff(a, b, co=1.3):
    """
    get the cutoff for identifying bonding
    :param a:
    :param b:
    :param co:
    :return:
    """
    return (CovalentRadius.radius[str(a)] + CovalentRadius.radius[str(b)]) * co


def identify_cnhybrid(site, mol):
    """
    we assume the molecule is neutral, now dealing with C N only
    :param site:
    :param mol:
    :return:
    """
    if site.species_string != 'C' and site.species_string != 'N':
        warnings.warn('trying to identify a non-carbon/nitrogen atom hybrid!')
        # sys.exit(1)
    elif site.species_string == 'C':
        nbs = get_nbs_in_mol(site, mol)
        if len(nbs) == 2:
            return 'csp'
        elif len(nbs) == 3:
            return 'csp2'
        elif len(nbs) == 4:
            return 'csp3'
        else:
            warnings.warn('a carbon site hybrid is weird!')
            sys.exit(1)
    elif site.species_string == 'N':
        nbs = get_nbs_in_mol(site, mol)
        if len(nbs) == 2:
            return 'nsp2'
        elif len(nbs) == 3:
            return 'nsp3'
        elif len(nbs) == 1:
            return 'nsp'
        else:
            warnings.warn('a nitrogen site hybrid is weird!')
            sys.exit(1)
    else:
        return 'not hybrid'


def get_nbs_in_mol(site, mol):
    """
    get a list of nbs of a site in a certain mol
    :param site:
    :param mol:
    :return: a list of nbs of a site in a certain mol
    """
    nbs = []
    for i in mol.sites:
        if not np.allclose(i.coords, site.coords):
            cutoff = get_cutoff(site.species_string, i.species_string)
            vdiff = i.coords - site.coords
            distance = math.sqrt(vdiff[0]**2 + vdiff[1]**2 + vdiff[2]**2)
            if distance < cutoff:
                nbs.append(i)
    if len(nbs) == 0:
        warnings.warn('isolated site detected!')
    return nbs


def mol_nrnbr_table(mol):
    """
    table[i] = list of nb index
    :param mol:
    :return:
    """
    table = []
    for i in range(len(mol.sites)):
        nbs = []
        s1 = mol.sites[i]
        for j in range(len(mol.sites)):
            if j != i:
                s2 = mol.sites[j]
                if 0.00001 < s1.distance(s2) < get_cutoff(s1.species_string, s2.species_string):
                    nbs.append(j)
        if len(nbs) == 0:
            warnings.warn('isolated site detected!')
        table.append(nbs)
    return table


def group_close_sites(sites):
    """
    percolation on msites
    :param sites
    :return:
    """
    indices = range(len(sites))
    visited = []
    side_chains = []
    block_list = []
    while len(visited) != len(sites):
        # initialization
        unvisited = [idx for idx in indices if idx not in visited]
        ini_idx = unvisited[0]
        block = [ini_idx]
        pointer = 0
        while pointer != len(block):
            outside = [idx for idx in indices if idx not in block and idx not in visited]
            for i in outside:
                distance = sites[block[pointer]].distance(sites[i])
                cutoff = get_cutoff(sites[block[pointer]].species_string, sites[i].species_string)
                if distance < cutoff:
                    block.append(i)
            visited.append(block[pointer])
            pointer += 1
        block_list.append(block)
        side_chain = Molecule.from_sites([sites[i] for i in block])
        side_chains.append(side_chain)
    return side_chains


def unwrapper(pstructure):
    """
    unwrap the structure, extract isolated mols
    :param pstructure: periodic structure obj from pymatgen
    :return:
    """
    psites = pstructure.sites
    pindices = range(len(psites))
    visited = []
    block_list = []
    unwrap = []
    unwrap_block_list = []
    while len(visited) != len(psites):
        # initialization
        unvisited = [idx for idx in pindices if idx not in visited]
        ini_idx = unvisited[0]
        block = [ini_idx]
        unwrap.append(psites[ini_idx])
        unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx]._coords)]
        # unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords)]
        pointer = 0
        while pointer != len(block):
            outside = [idx for idx in pindices if idx not in block and idx not in visited]
            for i in outside:
                distance, fctrans = get_dist_and_trans(psites[block[pointer]]._fcoords, psites[i]._fcoords,
                                                       pstructure.lattice)
                # distance, fctrans = get_dist_and_trans(psites[block[pointer]].fcoords, psites[i].fcoords,
                #                                        pstructure.lattice)
                cutoff = get_cutoff(psites[block[pointer]].species_string, psites[i].species_string)
                if distance < cutoff:
                    block.append(i)
                    psites[i] = PeriodicSite(psites[i].species_string, psites[i]._fcoords + fctrans,
                    pstructure.lattice)
                    unwrap_block.append(Site(psites[i].species_string, psites[i]._coords))
                    # psites[i] = PeriodicSite(psites[i].species_string, psites[i].fcoords + fctrans, pstructure.lattice)
                    # unwrap_block.append(Site(psites[i].species_string, psites[i].coords))
                    unwrap.append(psites[i])
            visited.append(block[pointer])
            pointer += 1
        unwrap_block_list.append(unwrap_block)
        block_list.append(block)
    mols = [IMolecule.from_sites(i) for i in unwrap_block_list]
    # for i in range(len(mols)):
    #     mols[i].to(fmt='xyz', filename='mol_' + str(i) + '.xyz')
    unwrap = sorted(unwrap, key=lambda x: x.species_string)
    unwrap_str = IStructure.from_sites(unwrap)
    unwrap_str.to(fmt='poscar', filename='unwrap.poscar')
    return mols, unwrap_str


'''
def percolation(msites, rule):
    """
    a generic percolation scheme, just used as a reference in case you need to recode
    :param msites: an array-like of site, better to use immutable
    :param rule: rule[i][j] gives boolean for percolation
    :return: a list of blocks
    """
    pindices = range(len(msites))
    visited = []
    block_list = []
    while len(visited) != len(msites):
        # initialization
        unvisited = [idx for idx in pindices if idx not in visited]
        ini_idx = unvisited[0]
        block = [ini_idx]
        pointer = 0
        while pointer != len(block):
            outside = [idx for idx in pindices if idx not in block and idx not in visited]
            for i in outside:
                if rule[block[pointer], i]:
                    block.append(i)
            visited.append(block[pointer])
            pointer += 1
        block_list.append(block)
    return block_list
'''
