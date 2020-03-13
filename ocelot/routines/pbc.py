import warnings
from copy import deepcopy

import numpy as np
from pymatgen.core.structure import Lattice
from pymatgen.core.structure import Molecule
from pymatgen.core.structure import PeriodicSite
from pymatgen.core.structure import Site
from pymatgen.core.structure import Structure
from pymatgen.util.coord import pbc_shortest_vectors

from ocelot.schema.conformer import Element

"""
PBCparser: get unwrapped structure and mols

method here should not rely on ocelot schema
"""


class PBCparser:

    @staticmethod
    def get_dist_and_trans(lattice: Lattice, fc1, fc2):
        """
        get the shortest distance and corresponding translation vector between two frac coords

        :param lattice: pmg lattic obj
        :param fc1:
        :param fc2:
        :return:
        """
        v, d2 = pbc_shortest_vectors(lattice, fc1, fc2, return_d2=True)
        if len(np.array(fc1).shape) == 1:
            fc = lattice.get_fractional_coords(v[0][0]) + fc1 - fc2
        else:
            fc = np.zeros((len(fc1), len(fc2), 3))
            for i in range(len(fc1)):
                for j in range(len(fc2)):
                    fc_vector = np.dot(v[i][j], lattice.inv_matrix)
                    fc[i][j] = fc_vector + fc1[i] - fc2[j]
        return np.sqrt(d2), fc

    @staticmethod
    def unwrap(pstructure: Structure):
        """
        unwrap the structure, extract isolated mols

        this will modify psites *in-place*, properties will be inherited and a new property
        'imol' will be written
        psite with imol=x is an element of both mols[x] and unwrap_pblock_list[x]

        this method is not supposed to modify siteid!

        :param pstructure: periodic structure obj from pymatgen
        :return: mols, unwrap_str_sorted, unwrap_pblock_list
        """

        psites = pstructure.sites
        cutoff_matrix = np.zeros((len(psites), len(psites)))
        for i in range(len(psites)):
            for j in range(i + 1, len(psites)):
                cutoff = Element(psites[i].species_string).atomic_radius + Element(
                    psites[j].species_string).atomic_radius
                cutoff *= 1.3
                cutoff_matrix[i][j] = cutoff
                cutoff_matrix[j][i] = cutoff

        if all('iasym' in s.properties for s in psites):
            sameiasym = True
            warnings.warn('if a percolating step involves hydrogens, only percolate to sites with the same iasym')
        else:
            sameiasym = False

        pindices = range(len(psites))
        visited = []
        block_list = []
        # unwrap = []
        unwrap_block_list = []
        unwrap_pblock_list = []
        while len(visited) != len(psites):
            # initialization
            unvisited = [idx for idx in pindices if idx not in visited]
            ini_idx = unvisited[0]
            block = [ini_idx]
            # unwrap.append(psites[ini_idx])
            unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords,
                                 properties=deepcopy(psites[ini_idx].properties))]
            unwrap_pblock = [psites[ini_idx]]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    si = psites[i]
                    sj = psites[block[pointer]]
                    if sameiasym and (si.species_string == 'H' or sj.species_string == 'H'):
                        if si.properties['iasym'] != sj.properties['iasym']:
                            continue
                    distance, fctrans = PBCparser.get_dist_and_trans(pstructure.lattice,
                                                                     psites[block[pointer]].frac_coords,
                                                                     psites[i].frac_coords, )
                    fctrans = np.rint(fctrans)

                    cutoff = cutoff_matrix[block[pointer]][i]
                    if distance < cutoff:
                        block.append(i)
                        psites[i]._frac_coords += fctrans
                        # psites[i] = PeriodicSite(psites[i].species_string, psites[i].frac_coords + fctrans,
                        #                          pstructure.lattice, properties=deepcopy(psites[i].properties))
                        unwrap_block.append(
                            # Site(psites[i].species_string, psites[i].coords, properties=deepcopy(psites[i].properties))
                            Site(psites[i].species_string, psites[i].coords, properties=psites[i].properties)
                        )
                        unwrap_pblock.append(psites[i])
                visited.append(block[pointer])
                pointer += 1
            unwrap_block_list.append(unwrap_block)
            unwrap_pblock_list.append(unwrap_pblock)
            block_list.append(block)

        unwrap = []
        for i in range(len(unwrap_block_list)):
            for j in range(len(unwrap_block_list[i])):
                unwrap_block_list[i][j].properties['imol'] = i
                unwrap_pblock_list[i][j].properties['imol'] = i
                unwrap.append(unwrap_pblock_list[i][j])

        mols = [Molecule.from_sites(sites) for sites in unwrap_block_list]

        # unwrap_structure = Structure.from_sites(sorted(unwrap, key=lambda x: x.species_string))
        unwrap_structure = Structure.from_sites(unwrap, to_unit_cell=False)
        return mols, unwrap_structure, unwrap_pblock_list

    @staticmethod
    def unwrap_and_squeeze(pstructure: Structure):
        """
        after unwrapping, the mols can be far away from each other, this tries to translate them s.t. they stay together

        :param pstructure:
        :return:
        """
        mols, unwrap_structure, psiteblocks = PBCparser.unwrap(pstructure)
        if len(mols) > 1:
            refpoint = mols[0].center_of_mass
            refpoint = unwrap_structure.lattice.get_fractional_coords(refpoint)
            for i in range(1, len(mols)):
                varmol = mols[i]
                varpoint = varmol.center_of_mass
                varpoint = unwrap_structure.lattice.get_fractional_coords(varpoint)
                distance, fctrans = PBCparser.get_dist_and_trans(unwrap_structure.lattice, refpoint, varpoint)
                for j in range(len(psiteblocks[i])):
                    psiteblocks[i][j]._frac_coords += fctrans
            squeeze_psites = []
            squeeze_mols = []
            pblock: [PeriodicSite]
            for pblock in psiteblocks:
                squeeze_block = []
                for ps in pblock:
                    squeeze_psites.append(ps)
                    squeeze_block.append(
                        Site(ps.species_string, ps.coords, properties=ps.properties))  # do we need deepcopy?
                mol = Molecule.from_sites(squeeze_block)
                squeeze_mols.append(mol)
            squeeze_unwrap_structure = Structure.from_sites(squeeze_psites)
            return squeeze_mols, squeeze_unwrap_structure, psiteblocks


        else:
            return mols, unwrap_structure, psiteblocks
