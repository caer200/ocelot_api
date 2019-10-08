from pymatgen.io.vasp.sets import DictSet, loadfn
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core.structure import Structure, Molecule
from ocelot.routines.pbc import CIFparser
from ocelot.schema.config import Config
from ocelot.task.idchg import IdChg
from ocelot.routines.fileop import createdir, movefile, copyfile, removefile
import numpy as np

"""
calculate lattice energy for both neutral and ionic crystals

ionic crystals: see Chen2017 (Phys. Chem. Chem. Phys., 2017, 19, 4114)
    - encut = 1000 eV
    - monk kmesh
    - reciprocal density 64
    - DFT-D3 code for single point 3-body dispersion
    - HSE single point
    - EUGEN for Madelung constant
    
1. pbc relaxation with pbe-d3bj
2. pbc singlepoint with hse06-d3bj
3. molecule relaxation with pbe-d3bj with charge cal by mopac lewis
4. molecule singlepoint with hse06-d3bj

notes:
- a1 = 0.383, a2 = 5.685, s6=1.000, s8 = 2.310 for HSE06-D3BJ as in J. Phys. Chem. C 2014, 118, 7615−7621


"""


def boxstructure(pmgmol):
    cartmat = pmgmol.cart_coords
    maxdelta = np.amax(cartmat, 0) - np.amin(cartmat, 0)
    a, b, c = maxdelta + 20
    return pmgmol.get_boxed_structure(a, b, c)


def load_yaml_config(fn):
    config = loadfn(str(("%s.yaml" % fn)))
    config["INCAR"].update(loadfn("VASPIncarBase.yaml"))
    return config


def poscar2structure(fn):
    with open(fn, 'r') as f:
        s = f.read()
    return Structure.from_str(s, 'poscar')


MOPACCMD = 'singularity exec /home/ai/bin/mopac.sif /opt/mopac/MOPAC2016.exe'

onekpoints = Kpoints()

class LatEng:

    def __init__(self, jobname, config, rootdir=None):
        """
        - rootdir
            - 0_mopac
            - 1_pbc_relax
            - 2_pbc_sp
            - 3_mol_relax
                - mol-1
                - mol-2
            - 4_mol_sp
                - mol-1
                - mol-2

        :param jobname:
        :param config:
        :param rootdir:
        """
        if rootdir is None:
            rootdir = jobname
        self.rootdir = rootdir
        mopacwdir = self.rootdir + '/' + '0_mopac'
        createdir(self.rootdir)
        createdir(mopacwdir)
        self.jobname = jobname
        self.config = config
        self.structure = config.unwrap_structure
        self.pmgmols = config.mols
        self.mol_charges = []

        for i in range(len(self.pmgmols)):
            mol = self.pmgmols[i]
            idchgname = '{}_mol-{}'.format(self.jobname, i)
            idchg = IdChg(mol, idchgname, MOPACCMD + ' ' + idchgname)
            idchg.write_inp(mopacwdir)
            d = idchg.run(mopacwdir)
            total_charge = d['total_charge']
            self.mol_charges.append(total_charge)

        self.unique_mols = []
        self.unique_charges = []
        for i in range(len(self.pmgmols)):
            mol = self.pmgmols[i]
            charge = self.mol_charges[i]
            compos = [m.composition for m in self.unique_mols]
            if mol.composition not in compos:
                self.unique_mols.append(mol)
                self.unique_charges.append(charge)

    def gen_pbc_rlx(self, wdir):
        config = load_yaml_config('lateng_pbe_d3bj_relax')
        vasp_set = DictSet(self.structure, config)
        vasp_set.write_input(wdir)

    def gen_pbc_sp(self, wdir):
        config = load_yaml_config('lateng_hse_d3bj_sp')
        vasp_set = DictSet(self.structure, config)
        vasp_set.write_input(wdir)

    def gen_mol_rlx(self, pmgmol, charge, wdir):
        config = load_yaml_config('lateng_pbe_d3bj_relax')
        boxs = boxstructure(pmgmol)
        boxs._charge = charge
        vasp_set = DictSet(boxs, config, use_structure_charge=True, user_kpoints_settings=onekpoints)
        vasp_set.write_input(wdir)

    def gen_mol_sp(self, pmgmol, charge, wdir):
        config = load_yaml_config('lateng_hse_d3bj_sp')
        boxs = boxstructure(pmgmol)
        boxs._charge = charge
        vasp_set = DictSet(boxs, config, use_structure_charge=True, user_kpoints_settings=onekpoints)
        vasp_set.write_input(wdir)

    def gen_all(self):
        wdir_1_pbc_relax = self.rootdir + '/1_pbc_relax'
        createdir(wdir_1_pbc_relax)
        self.gen_pbc_rlx(wdir_1_pbc_relax)

        wdir_2_pbc_sp = self.rootdir + '/2_pbc_sp'
        createdir(wdir_2_pbc_sp)
        self.gen_pbc_sp(wdir_2_pbc_sp)
        removefile(wdir_2_pbc_sp + '/POSCAR')

        wdir_3_mol_relax = self.rootdir + '/3_mol_relax'
        createdir(wdir_3_mol_relax)

        wdir_4_mol_sp = self.rootdir + '/4_mol_sp'
        createdir(wdir_4_mol_sp)

        for i in range(len(self.unique_charges)):
            wdir_3_i = wdir_3_mol_relax + '/mol-{}'.format(i)
            mol = self.unique_mols[i]
            charge = self.unique_charges[i]
            self.gen_mol_rlx(mol, charge, wdir_3_i)

            wdir_4_i = wdir_4_mol_sp + '/mol-{}'.format(i)
            self.gen_mol_sp(mol, charge, wdir_4_i)
            removefile(wdir_4_i + '/POSCAR')


    @classmethod
    def from_rawcif(cls, ciffn, rootdir=None):
        name = ciffn.split('/')[-1]
        name = name[:-4]
        with open(ciffn, 'r') as f:  # read cif file into string
            fs = f.read()

        cp = CIFparser.from_cifstring(fs)  # create a parser object from cif string
        clean_cif_strings = cp.get_clean_cifs_stringlist()  # a list of cif strings without disorder
        config = Config.from_cifstring(clean_cif_strings[0])
        return cls(name, config, rootdir)


import glob
for i in glob.glob('./ptsalts/*.cif'):
    jobname = i.split('/')[-1].split('.')[0]
    lateng = LatEng.from_rawcif(i)
    lateng.gen_all()


# lateng_pbc_relax = load_yaml_config('lateng_pbc_relax')
# lateng_pbc_sp = load_yaml_config('lateng_pbc_sp')
#
# structure = Structure.from_file('POSCAR')
#
# latset1 = DictSet(structure, lateng_pbc_relax)
# latset1.write_input('./')
