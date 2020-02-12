import math
from copy import deepcopy

from pymatgen.core.structure import Site, PeriodicSite, Molecule, Structure
from pymatgen.util.coord import pbc_shortest_vectors

from ocelot.schema.conformer import MolConformer, Element

"""
PBCparser: get unwrapped structure and mols
"""


class PBCparser:
    #
    # def __init__(self, pstructure):
    #     self.nsites = len(pstructure)
    #     self.structure = pstructure

    @staticmethod
    def get_dist_and_trans(lattice, fc1, fc2):
        """
        get the shortest distance and corresponding translation vector between two frac coords

        :param lattice: pmg lattic obj
        :param fc1:
        :param fc2:
        :return:
        """
        v, d2 = pbc_shortest_vectors(lattice, fc1, fc2, return_d2=True)
        fc = lattice.get_fractional_coords(v[0][0]) + fc1 - fc2
        return math.sqrt(d2[0, 0]), fc

    @staticmethod
    def unwrap(pstructure):
        """
        unwrap the structure, extract isolated mols

        :param pstructure: periodic structure obj from pymatgen
        :return:
        """
        psites = pstructure.sites
        for isite in range(len(psites)):
            psites[isite].properties['siteid'] = isite
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
            unwrap_block = [Site(psites[ini_idx].species_string, psites[ini_idx].coords, properties=deepcopy(psites[ini_idx].properties))]
            unwrap_pblock = [psites[ini_idx]]
            pointer = 0
            while pointer != len(block):
                outside = [idx for idx in pindices if idx not in block and idx not in visited]
                for i in outside:
                    distance, fctrans = PBCparser.get_dist_and_trans(pstructure.lattice,
                                                                     psites[block[pointer]]._frac_coords,
                                                                     psites[i]._frac_coords, )

                    cutoff = Element(psites[block[pointer]].species_string).atomic_radius + Element(
                        psites[i].species_string).atomic_radius
                    cutoff *= 1.3
                    if distance < cutoff:
                        block.append(i)
                        psites[i] = PeriodicSite(psites[i].species_string, psites[i]._frac_coords + fctrans,
                                                 pstructure.lattice, properties=deepcopy(psites[i].properties))
                        unwrap_block.append(
                            Site(psites[i].species_string, psites[i].coords, properties=deepcopy(psites[i].properties)))
                        # unwrap.append(psites[i])
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

        # this does not work, from_sites cannot pickup properties
        mols = [Molecule.from_sites(sites) for sites in unwrap_block_list]
        # mols = []
        # for group in unwrap_block_list:
        #     property_list = []
        #     for i in range(len(group)):
        #         property_list.append(deepcopy(group[i].properties))
        #         group[i].properties = {}
        #     mol = Molecule.from_sites(group)
        #     for i in range(len(group)):
        #         mol._sites[i].properties = property_list[i]
        #     mols.append(mol)

        unwrap = sorted(deepcopy(unwrap), key=lambda x: x.species_string)
        unwrap_str_sorted = Structure.from_sites(unwrap)
        return mols, unwrap_str_sorted, unwrap_pblock_list

    @staticmethod
    def squeeze(pstructure: Structure):
        """
        after unwrapping, the mols can be far away from each other, this tries to translate them s.t. they stay together

        :rtype:
        :param pstructure:
        :return:
        """
        mols, unwrap_structure, psiteblocks = PBCparser.unwrap(pstructure)

        moleconformers = []
        for mm in mols:
            mc = MolConformer.from_pmgmol(mm)
            if not mc.is_solvent:
                moleconformers.append(mc)

        if len(moleconformers) > 1:
            refpoint = moleconformers[0].backbone.geoc
            refpoint = unwrap_structure.lattice.get_fractional_coords(refpoint)
            for i in range(1, len(moleconformers)):
                varmol = moleconformers[i]
                varpoint = varmol.backbone.geoc
                varpoint = unwrap_structure.lattice.get_fractional_coords(varpoint)
                distance, fctrans = PBCparser.get_dist_and_trans(unwrap_structure.lattice, refpoint, varpoint)
                for j in range(len(psiteblocks[i])):
                    psiteblocks[i][j]._frac_coords += fctrans
            psites = []
            mols = []
            for pblock in psiteblocks:
                block = []
                for ps in pblock:
                    psites.append(ps)
                    block.append(Site(ps.species_string, ps.coords))
                mol = Molecule.from_sites(block)
                mols.append(mol)
            unwrap_structure = Structure.from_sites(sorted(psites, key=lambda x: x.species_string))

            moleconformers = []
            for m in mols:
                mc = MolConformer.from_pmgmol(m)
                if not mc.is_solvent:
                    moleconformers.append(mc)
            return mols, moleconformers, unwrap_structure

        return mols, moleconformers, unwrap_structure
