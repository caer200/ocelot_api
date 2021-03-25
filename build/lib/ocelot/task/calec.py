import re
import sys
from collections import OrderedDict

import numpy as np
import numpy.matlib
from numba import jit
from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput
from scipy.linalg import blas
from scipy.linalg import lapack

'''
electronic coupling calculator with

1. gaussian, code comes from Sean M. Ryno, see original code at /legacy/sean-electronic_coupling.py for citation suggestions

2. zindo, under dev
'''

SUGGESTED_route_parameters = {
    'nosymm': '',
    'pop': 'full',
    'iop': {'3/33': '1', '3/59': '8'},
    'int': {'acc2e': '12'},
}
SUGGESTED_route_parameters = OrderedDict(SUGGESTED_route_parameters)


class ElectronicCoupling:
    """
    class to calculate electronic couplings from Gaussian w. fmo formalism

    micro wf:

    1. get 2 monomer structure in cart coords, label them mono_a and mono_b

    2. concatenate 2 input structure to give dimer structure, fmo requires site sequence remains intact

    3. run calculations to give mono_a_fchk, mono_b_fchk, dimer_fchk, dimer_log.

    suggested keywords `nosymm pop=full iop(3/33=1,3/59=8) int(acc2e=12)`

    4. four files in 3. will be used to give ec, the results is a list of lists,

    fields are ['#mo_a', '#mo_b', 'e_mo_a', 'e_mo_b', 'J_12', 'S_12', 'e^eff_1', 'e^eff_2', 'J^eff_12', 'dE_12']

    all in eV
    """

    def __init__(self, sysname, mono_a, mono_b, caltype='gauss'):
        """
        setup ec calculations
        :param sysname:
        :param mono_a: pmg mol obj
        :param mono_b: pmg mol obj
        :param caltype: 'gauss' or 'zindo'
        """

        self.sysname = sysname
        self.mono_a = mono_a
        self.mono_b = mono_b
        self.caltype = caltype

    def input_gen_gauss(
            self, charge=None, spin_multiplicity=None, title=None, functional='HF', basis_set='6-31G(d)',
            route_parameters=SUGGESTED_route_parameters, input_parameters=None, link0_parameters=None, dieze_tag='#P',
            gen_basis=None
    ):

        ginobj_mono_a = GaussianInput(self.mono_a, charge, spin_multiplicity, title, functional, basis_set,
                                      route_parameters, input_parameters, link0_parameters, dieze_tag, gen_basis)

        ginobj_mono_b = GaussianInput(self.mono_a, charge, spin_multiplicity, title, functional, basis_set,
                                      route_parameters, input_parameters, link0_parameters, dieze_tag, gen_basis)

        dimer_sites = self.mono_a.sites + self.mono_b.sites
        dimer_pmgmol = Molecule.from_sites(dimer_sites)
        ginobj_dimer = GaussianInput(dimer_pmgmol, charge, spin_multiplicity, title, functional, basis_set,
                                     route_parameters, input_parameters, link0_parameters, dieze_tag, gen_basis)
        return ginobj_mono_a.to_string(cart_coords=True), ginobj_mono_b.to_string(
            cart_coords=True), ginobj_dimer.to_string(cart_coords=True)

    @staticmethod
    def calculate_coupling_gauss(mono_a_fchk_string, mono_b_fchk_string, dimer_fchk_string, dimer_log_string,
                                 homo1=0, lumo1=0, homo2=0, lumo2=0):
        monomer_a_fchk_lines = mono_a_fchk_string.split('\n')
        monomer_b_fchk_lines = mono_b_fchk_string.split('\n')
        dimer_fchk_lines = dimer_fchk_string.split('\n')
        dimer_log_lines = dimer_log_string.split('\n')
        # Sanity check of homos and lumos
        if (homo1 > lumo1) or (homo2 > lumo2):
            sys.exit("E: Homo and lumo cannot be same of overlap.")

        # Get dimer values
        nbf_d = ElectronicCoupling.g09_get_basis_functions(dimer_fchk_lines)
        e_d = ElectronicCoupling.g09_get_alpha_orbital_energies(dimer_fchk_lines, nbf_d)
        c_d = ElectronicCoupling.g09_get_mo_coefficients(dimer_fchk_lines, nbf_d)

        s_d = ElectronicCoupling.get_overlap(dimer_log_string, nbf_d)
        f_d = ElectronicCoupling.g09_calculate_fock(c_d, s_d, nbf_d, e_d)

        # Get monomer A values
        nbf_monomer_a = ElectronicCoupling.g09_get_basis_functions(monomer_a_fchk_lines)
        homo_monomer_a, lumo_monomer_a = ElectronicCoupling.g09_get_homo_lumo(monomer_a_fchk_lines)
        homo_monomer_a += homo1
        lumo_monomer_a += lumo1
        c_monomer_a = ElectronicCoupling.g09_get_mo_coefficients(monomer_a_fchk_lines, nbf_monomer_a)

        # Get monomer B values
        nbf_monomer_b = ElectronicCoupling.g09_get_basis_functions(monomer_b_fchk_lines)
        homo_monomer_b, lumo_monomer_b = ElectronicCoupling.g09_get_homo_lumo(monomer_b_fchk_lines)
        homo_monomer_b += homo2
        lumo_monomer_b += lumo2
        c_monomer_b = ElectronicCoupling.g09_get_mo_coefficients(monomer_b_fchk_lines, nbf_monomer_b)

        # Calculate Coupling
        cfrag = ElectronicCoupling.g09_build_cfrag(nbf_d, nbf_monomer_a, nbf_monomer_b, c_monomer_a, c_monomer_b)
        out = ElectronicCoupling.g09_calculate_coupling(homo_monomer_a,
                                                        lumo_monomer_a, homo_monomer_b, lumo_monomer_b, nbf_monomer_a,
                                                        f_d, s_d, cfrag)
        return out

    @staticmethod
    @jit
    def g09_calculate_coupling(homo_monomer_a, lumo_monomer_a, homo_monomer_b, lumo_monomer_b, nbf_monomer_a,
                               f_d, s_d, cfrag):
        """

        :param homo_monomer_a:
        :param lumo_monomer_a:
        :param homo_monomer_b:
        :param lumo_monomer_b:
        :param nbf_monomer_a:
        :param f_d:
        :param s_d:
        :param cfrag:
        :return: "MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, eV", "S_12",
                "e^eff_1, eV", "e^eff_2, eV", "J^eff_12, eV", "dE_12, eV"
        """
        irange = tuple(range(homo_monomer_a - 1, lumo_monomer_a))
        jrange = tuple(range(homo_monomer_b - 1 + nbf_monomer_a, lumo_monomer_b + nbf_monomer_a))
        out = []
        for i in irange:
            for j in jrange:
                ctmp = cfrag[i, :] * f_d
                e1 = blas.ddot(cfrag[i, :], ctmp) * 27.2116
                ctmp = cfrag[j, :] * f_d
                e2 = blas.ddot(cfrag[j, :], ctmp) * 27.2116
                ctmp = cfrag[j, :] * f_d
                j12 = blas.ddot(cfrag[i, :], ctmp) * 27.2116
                ctmp = cfrag[j, :] * s_d
                s12 = blas.ddot(cfrag[i, :], ctmp)
                ee1 = 0.5 * ((e1 + e2) - 2.0 * j12 * s12 + (e1 - e2) * np.sqrt(1.0 - s12 * s12)) / (1.0 - s12 * s12)
                ee2 = 0.5 * ((e1 + e2) - 2.0 * j12 * s12 - (e1 - e2) * np.sqrt(1.0 - s12 * s12)) / (1.0 - s12 * s12)
                je12 = (j12 - 0.5 * (e1 + e2) * s12) / (1.0 - s12 * s12)
                de = np.sqrt((ee1 - ee2) ** 2 + (2 * je12) ** 2)
                out.append((i + 1, j + 1 - nbf_monomer_a, e1, e2, j12, s12, ee1, ee2,
                            je12, de))
        return out

    @staticmethod
    def g09_get_basis_functions(fchk_lines):
        # Get number of basis functions
        for line in fchk_lines:
            if "Number of independent functions" in line:
                nbf = int(line.split()[5])
                return nbf

    @staticmethod
    def g09_get_homo_lumo(fchk_lines):
        # Get HOMO and LUMO levels
        for line in fchk_lines:
            if "Number of alpha electrons" in line:
                homo = int(line.split()[5])
                lumo = homo + 1
                return homo, lumo

    @staticmethod
    def g09_get_alpha_orbital_energies(fchk_lines, nbf):
        # Get alpha orbital energies
        e = np.matlib.zeros(nbf)
        line_index = -1
        for line in fchk_lines:
            if "Alpha Orbital Energies" in line:
                line_index = int(fchk_lines.index(line)) + 1
                break
        if line_index == -1:
            print('Alpha Orbital Enetrgies not present')
            sys.exit(1)
        count = 0
        icount = 0
        while count < nbf:
            mo_line = fchk_lines[line_index + icount].split()
            icount += 1
            for i in range(len(mo_line)):
                e[0, count] = float(mo_line[i])
                count += 1
        return e

    @staticmethod
    def g09_get_mo_coefficients(fchk_lines, nbf):
        # Get alpha MO coefficients
        c_temp = np.matlib.zeros(nbf * nbf)
        line_index = -1
        for line in fchk_lines:
            if "Alpha MO coefficients" in line:
                line_index = int(fchk_lines.index(line)) + 1
                break

        if line_index == -1:
            print('Alpha MO coefficients not present')
            sys.exit(1)

        count = 0
        icount = 0
        while count < nbf * nbf:
            mo_line = fchk_lines[line_index + icount].split()
            icount += 1
            for i in range(len(mo_line)):
                c_temp[0, count] = float(mo_line[i])
                count += 1
        c = c_temp.reshape((nbf, nbf))
        return c

    @staticmethod
    @jit
    def g09_calculate_fock(c, s, nbf, e):
        # Calculate the Fock matrix
        sc_temp = blas.dgemm(1.0, np.array(c), np.array(s), 1.0)
        sc = np.asmatrix(sc_temp)
        sce = np.matlib.zeros((nbf, nbf))
        for i in range(nbf):
            sce[i, :] = e[0, i] * sc[i, :]
        c_lu, ipiv, info = lapack.dgetrf(c)
        ic_lu, info = lapack.dgetri(c_lu, ipiv)
        f_temp = blas.dgemm(1.0, np.array(ic_lu), np.array(sce), 1.0)
        f = np.asmatrix(f_temp)
        return f

    @staticmethod
    @jit
    def g09_build_cfrag(nbf_d, nbf_mon1, nbf_mon2, c_mon1, c_mon2):
        # Build the block matrix of monomer MO coefficients
        cfrag = np.matlib.zeros((nbf_d, nbf_d))
        for i in range(nbf_mon1):
            for j in range(nbf_mon1):
                cfrag[j, i] = c_mon1[j, i]

        for i in range(nbf_mon2):
            for j in range(nbf_mon2):
                cfrag[j + nbf_mon1, i + nbf_mon1] = c_mon2[j, i]

        return cfrag

    @staticmethod
    def get_overlap(g09logstring, nbf):
        """
        Extracts the overlap matrix from a Gaussian logfile string.
        Returns a numpy matrix.
        """
        data = g09logstring
        overlap_matrix = np.zeros((nbf, nbf))
        # grab all text between  "Overlap ***" and "*** Kinetic"
        raw_overlap_string = re.findall(r'Overlap \*\*\*(.*?)\*\*\* Kinetic', data, re.DOTALL)
        raw_overlap_string = raw_overlap_string[0].replace('D', 'E')
        raw_overlap_elements = raw_overlap_string.split()
        matrix_elements = []
        for overlap_value in raw_overlap_elements:
            if 'E' in overlap_value:
                matrix_elements.append(overlap_value)
        overlap = ElectronicCoupling.create_matrix(overlap_matrix, matrix_elements, nbf)
        overlap = ElectronicCoupling.make_symmetric(overlap)
        overlap = np.array(overlap)

        return overlap

    @staticmethod
    def create_matrix(matrix, elements, nbf):
        """
        create lower triangular matrix from list of matrix elements indexed like so:

            [[0,0,0,...,0],
            [1,2,0,...,0],
            [3,4,5,...,0]]

        nbf is number of basis functions

        elements is a list of matrix elements indexed like above, e.g. [0,1,2,3,...]

        Gaussian prints every 5 columns, so the mod5 accounts for this

        """
        count = 0  # count is our index
        # fill the main block, leaving remainder for triangle fill
        for i in range(0, nbf - nbf % 5, 5):
            matrix, count = ElectronicCoupling.triangle_fill(matrix, i, i + 5, i, count, elements)
            matrix, count = ElectronicCoupling.block_fill(matrix, i + 5, nbf, i, i + 5, count, elements)
        # finish filling the last triangle bit
        matrix, count = ElectronicCoupling.triangle_fill(matrix, nbf - nbf % 5, nbf, nbf - nbf % 5, count, elements)

        return matrix

    @staticmethod
    @jit
    def make_symmetric(matrix):
        return matrix + matrix.T - np.diag(np.diag(matrix))

    @staticmethod
    @jit
    def triangle_fill(matrix, istrt, iend, jstrt, count, element):
        for i in range(istrt, iend):
            for j in range(jstrt, i + 1):
                matrix[i, j] = element[count]
                count += 1
        return matrix, count

    @staticmethod
    @jit
    def block_fill(matrix, istrt, nbf, jstrt, jend, count, element):
        for i in range(istrt, nbf):
            for j in range(jstrt, jend):
                matrix[i, j] = element[count]
                count += 1
        return matrix, count
