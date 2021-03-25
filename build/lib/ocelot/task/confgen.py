"""
given a molecular specification (smiles) get distinct conformers via rdkit

the molecule should already have hydrogen explicitly added, i.e. from hsmiles or addhs

workflow:
    rdkit conf gen [prune_in_gen] --> [ff minimize] --> prune_by_energy --> prune_by_rmsd

ref
http://rdkit.blogspot.com/2014/08/optimizing-diversity-picking-in-rdkit.html
https://squonk.it/docs/cells/RDKit%20Diverse%20Subset%20Picker/
https://github.com/skearnes/rdkit-utils/blob/master/rdkit_utils/conformers.py
https://new.pharmacelera.com/scripts/rdkit-conformation-generation-script/
https://www.rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
http://rdkit.org/docs_temp/Cookbook.html#parallel-conformation-generation
"""
__author__ = "qianxiang ai, parker sornberger"

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from ocelot.schema.rdfunc import RdFunc
import warnings


class ConfGenError(Exception): pass


class ConfGen:

    def __init__(self, m: Chem.Mol, mol_name='rdmol', nconf=None, mingen=50, maxgen=300, nthread=0, rngseed=-1,
                 prune_in_gen=None, minimize=True, ff='MMFF94s', prune_by_energy=100):
        self.mol_name = mol_name
        self.m = Chem.Mol(m)  # get a new mol
        self.m.RemoveAllConformers()
        self.mingen = mingen
        self.maxgen = maxgen
        self.rngseed = rngseed
        try:
            self.nconf = int(nconf)
        except TypeError:
            self.nconf = self.ideal_conf_num(self.m, self.mingen, self.maxgen)

        try:
            self.prune_in_gen = float(prune_in_gen)  # this is the threshold
        except TypeError:
            warnings.warn('prune_in_gen will not be carried out')
            self.prune_in_gen = -1.0

        self.minimize = minimize
        if ff not in ['UFF', 'MMFF94', 'MMFF94s']:
            raise ConfGenError('force field: {} not understood!'.format(ff))
        self.ff = ff

        self.nthread = nthread
        self.prune_by_energy = float(prune_by_energy)

    @staticmethod
    def is_hydrogen_explicit(m: Chem.Mol):
        if not any('h' == a.GetSymbol().lower() for a in m.GetAtoms()):
            raise ConfGenError('no explicit hydrogen in this rdkit mol!')

    @classmethod
    def from_smiles(cls, smiles: str, mingen=50, maxgen=300):
        if 'h' not in smiles.lower():
            raise ConfGenError('no h in the smiles!')
        m = Chem.MolFromSmiles(smiles)
        n = Chem.AddHs(m)
        return cls(n, mingen=mingen, maxgen=maxgen)

    @staticmethod
    def ideal_conf_num(some_molecule: Chem.Mol, min=50, max=300):
        """
        return the number of conformers that need to be generated to sample the conf space of a molecule
        """
        rotatable_bonds = int(AllChem.CalcNumRotatableBonds(some_molecule))
        if rotatable_bonds < 3:
            conf_number = min
        elif rotatable_bonds > 6:
            conf_number = max
        else:
            conf_number = (rotatable_bonds ** 3)
        return conf_number

    @staticmethod
    def genconfs_ETKDG(m: Chem.Mol, nconf, nthreads, rmsprune, rngseed):
        # confids = AllChem.EmbedMultipleConfs(m, numConfs=nconf, params=AllChem.ETKDG(), pruneRmsThresh=rmsprune,
        #                                  numThreads=nthreads, randomSeed=rngseed)
        params = AllChem.ETKDG()
        params.numThreads = nthreads
        params.randomSeed = rngseed
        params.pruneRmsThresh = rmsprune
        confids = AllChem.EmbedMultipleConfs(m, numConfs=nconf, params=params, )
        return m, list(confids)

    @staticmethod
    def cal_energy(m: Chem.Mol, cid: int = 0, ff='MMFF94s'):
        if 'MMFF' in ff:
            molecule_properties = AllChem.MMFFGetMoleculeProperties(m, mmffVariant=ff)
            force_field = AllChem.MMFFGetMoleculeForceField(m, molecule_properties, confId=cid)
        else:
            force_field = AllChem.UFFGetMoleculeForceField(m, confId=cid)
        return float(force_field.CalcEnergy())

    def genconfs_prune_by_energy(self):
        # energy in kcal/mol
        m, confids = self.genconfs_ETKDG(self.m, self.nconf, self.nthread, self.prune_in_gen, self.rngseed)
        if self.minimize:
            if 'MMFF' in self.ff:
                ff_results = AllChem.MMFFOptimizeMoleculeConfs(m, mmffVariant=self.ff, numThreads=self.nthread)
            elif 'UFF' in self.ff:
                ff_results = AllChem.UFFOptimizeMoleculeConfs(m, numThreads=self.nthread)
            else:
                raise ConfGenError('force field {} not understood!'.format(self.ff))
            energies = [r[1] for r in ff_results]
        else:
            energies = [self.cal_energy(m, cid, self.ff) for cid in confids]
        energies = np.array(energies, dtype=float)
        minenergy = min(energies)
        new = Chem.Mol(m)
        new.RemoveAllConformers()
        new_energies = []
        for i in confids:
            if energies[i] - minenergy < self.prune_by_energy:
                conf = m.GetConformer(confids[i])
                new.AddConformer(conf, assignId=True)
                new_energies.append(energies[i])
        return new, new_energies

    @staticmethod
    def prune_by_rmsd(m, energies, max_conformers=None, rmsd_threshold=0.5):
        """
        conformer in m should be assigned consecutive ids with delta=1

        :param energies:
        :param m:
        :param max_conformers:
        :param rmsd_threshold:
        :return:
        """
        rmsd_condensed = AllChem.GetConformerRMSMatrix(m)
        # rmsd_condensed = ConfGen.GetBestRMSMatrix(m)
        if max_conformers is None:
            max_conformers = len(m.GetConformers())
        # Returns the RMS matrix of the conformers of a molecule.
        # As a side-effect, the conformers will be aligned to the first conformer (i.e. the reference)
        # and will left in the aligned state.
        from scipy.spatial.distance import squareform
        rmsd_mat = squareform(rmsd_condensed)
        # print(rmsd_mat[0,2])
        # print(ConfGen.GetBestRMSConf(m, 0, 2))
        # print(ConfGen.GetBestRMSConf(m, 2, 0))
        atom_list = [a.GetSymbol() for a in m.GetAtoms()]
        # RdFunc.conf2xyz(m.GetConformer(2), '2.xyz', atom_list, '2')
        # RdFunc.conf2xyz(m.GetConformer(0), '0.xyz', atom_list, '0')
        keep = []
        discard = []
        confids = []
        for c in m.GetConformers():
            cid = c.GetId()
            confids.append(cid)
            # always keep lowest-energy conformer
            if len(keep) == 0:
                keep.append(cid)
                continue

            # discard conformers after max_conformers is reached
            if len(keep) >= max_conformers:
                discard.append(cid)
                continue

            # get RMSD to selected conformers
            this_rmsd = rmsd_mat[cid][np.asarray(keep, dtype=int)]

            # discard conformers within the RMSD threshold
            if np.all(this_rmsd >= rmsd_threshold):
                keep.append(cid)
            else:
                discard.append(cid)
        m_after_prune = Chem.Mol(m)
        m_after_prune.RemoveAllConformers()
        n_energies = []
        for i in keep:
            conf = m.GetConformer(i)
            n_energies.append(energies[i])
            m_after_prune.AddConformer(conf, assignId=True)
        return m_after_prune, n_energies

    @staticmethod
    def GetBestRMSConf(mol, i, j):
        from rdkit.Chem.rdMolAlign import GetBestRMS
        # Chem.AssignStereochemistry(mol,cleanIt=True,force=True,flagPossibleStereoCenters=True)
        mol_prob = Chem.Mol(mol)
        # from collections import defaultdict
        # equivs = defaultdict(set)
        # matches = m.GetSubstructMatches(m,uniquify=False)
        # for match in matches:
        #     for idx1,idx2 in enumerate(match): equivs[idx1].add(idx2)
        # print(equivs)


        # from copy import deepcopy
        # mol_prob = deepcopy(mol)
        # map = [[(a.GetIdx(), a.GetIdx()) for a in mol.GetAtoms()]]
        # [atom.GetProp('_CIPRank') for atom in m.GetAtoms()]
        # map = []
        # for a in mol.GetAtoms():
        #     symm_eqs = [a.GetIdx()]
        #     for b in mol.GetAtoms():
        #         if b.GetIdx() != a.GetIdx() and a.GetProp('_CIPRank') == b.GetProp('_CIPRank'):
        #             symm_eqs.append(b.GetIdx())
        #     map.append(symm_eqs)
        # print(map)
        # combs = []
        # for i in range(len(map)):
        #     possible_ids = map[i]
        #     for pid in possible_ids:
        #         if pid not in combs:
        #             combs.append(pid)
        #             break
        # import itertools
        # map = list(itertools.product(*map))
        # print(len(map))
        # print([a.GetSymbol() for a in mol.GetAtoms()])
        # print(map)
        return GetBestRMS(mol_prob, mol, prbId=j, refId=i,)# maxMatches=2000000)


    @staticmethod
    def GetBestRMSMatrix(mol):
        """
        the problem of AllChem.GetConformerRMSMatrix(m) is that it uses AlignMolConformers
        this assumes atom order is the same, which may not be the case if the molecule has symmetry
        """
        from rdkit.Chem.rdMolAlign import GetBestRMS
        mol_prob = Chem.Mol(mol)
        rmsvals = []
        confIds = [conf.GetId() for conf in mol.GetConformers()]
        for i in range(1, len(confIds)):
            rmsvals.append(
                GetBestRMS(mol_prob, mol, prbId=i, refId=0)
            )
        # loop over the conformations (except the reference one)
        cmat = []
        for i in range(1, len(confIds)):
            cmat.append(rmsvals[i - 1])
            for j in range(1, i):
                cmat.append(GetBestRMS(mol_prob, mol, prbId=j, refId=i))
        return np.array(cmat, dtype=float)

    @staticmethod
    def cluster_by_rmsd(m: Chem.Mol, energies, rmsd_threshold: float):
        """
        conformer in m should be assigned consecutive ids with delta=1

        :param energies:
        :param m:
        :param rmsd_threshold:
        :return:
        """
        rmsd_condensed = AllChem.GetConformerRMSMatrix(m)
        # rmsd_condensed = ConfGen.GetBestRMSMatrix(m)
        from rdkit.ML.Cluster.Butina import ClusterData
        clusters = ClusterData(rmsd_condensed, len(m.GetConformers()), distThresh=rmsd_threshold, isDistData=True)
        m_after_prune = Chem.Mol(m)
        m_after_prune.RemoveAllConformers()
        n_energies = []
        for c in clusters:
            conf = m.GetConformer(c[0])
            n_energies.append(energies[c[0]])
            m_after_prune.AddConformer(conf, assignId=True)
        return m_after_prune, n_energies

    def genconfs(self, flavor="prune", write_confs=False, prune_max_conformer=20, rmsd_threshold=0.5):
        m_prune_by_eng, energies = self.genconfs_prune_by_energy()
        if flavor == "prune":
            m_after, e_after = self.prune_by_rmsd(m_prune_by_eng, energies=energies, max_conformers=prune_max_conformer,
                                                  rmsd_threshold=rmsd_threshold)
        elif flavor == 'cluster':
            m_after, e_after = self.cluster_by_rmsd(m_prune_by_eng, energies=energies,
                                                    rmsd_threshold=rmsd_threshold)
        else:
            raise ConfGenError('flavor {} not understood'.format(flavor))
        if write_confs:
            atom_list = [a.GetSymbol() for a in m_after.GetAtoms()]
            for c in m_after.GetConformers():
                cid = c.GetId()
                energy = e_after[cid]
                comment = '{} energy: {:.4f} kcal/mol'.format(self.ff, energy)
                RdFunc.conf2xyz(c, "{}_{}.xyz".format(self.mol_name, cid), atom_list, comment)
        return m_after, e_after

    @staticmethod
    def show_confs_diversity(m: Chem.Mol, energies: [float]):
        # TODO find a way to display the diversity in/clusters of the conformers generated
        #  m should be a rdkit molecule with conformers generated and confomer ids assigned
        #  s.t. energies[conformer_id] = energy
        #  possible descriptor of diversity:
        #       energies,
        #       shape descriptors, e.g. AllChem.ShapeTanimotoDist
        #       rmsd spread, e.g. https://www.rdkit.org/docs/Cookbook.html#rmsd-calculation-between-n-molecules
        #       dimension reduction method for rmsd matrix, e.g. http://sketchmap.org/
        #  notice we can take the advantage that all conformers share the same molecule, i.e. atom order does not matter
        pass

