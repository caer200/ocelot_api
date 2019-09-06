import sys
import re
import numpy as np
import numpy.matlib
from scipy.linalg import blas
from scipy.linalg import lapack
from numba import jit
'''
                        Electronic Coupling Program                         
                               Sean M. Ryno                                 
                              June 13, 2017                                 
 Note: This program makes use of components of the gauss_parse.py           
       program of J. J. Goings                                              
 Based upon: E. F. Valeev, V. Coropceanu, D. A. da Silva Filho, S. Salman,  
             J. L. Bredas, J. Am. Chem. Soc. 128, 9882-9886 (2006)          
'''


__author__ = 'Sean M. Ryno'
__copyright__ = 'Copyright 2017, Sean M. Ryno'
__credits__ = ['Sean M. Ryno', 'J. J. Goings']
__license__ = 'GPL v3.0'
__version__ = '1.00'
__maintainer__ = 'Sean M. Ryno'
__email__ = 'sean.m.ryno@gmail.com'
__status__ = 'Development'


class ElectronicCoupling:
    """ class to calculate electronic couplings from Gaussian output files."""
    def __init__(self, sysname, homo1=0, lumo1=0, homo2=0, lumo2=0):
        self.sysname = sysname
        # Homo/Lumo values are in addition to those extacted. Example: for homo-1 to lumo+1, homo=-1, lumo=1
        self.homo1 = homo1
        self.lumo1 = lumo1
        self.homo2 = homo2
        self.lumo2 = lumo2

        with open(self.sysname + '_A.fchk', 'r') as f:
            self.monomer_a_fchk_lines = f.readlines()
        with open(self.sysname + '_B.fchk', 'r') as f:
            self.monomer_b_fchk_lines = f.readlines()
        with open(self.sysname + '_A_B.fchk', 'r') as f:
            self.dimer_fchk_lines = f.readlines()

        # Sanity check of homos and lumos
        if (self.homo1 > self.lumo1) or (self.homo2 > self.lumo2):
            print("Homo and lumo cannot be same of overlap.", file=sys.stderr)
            sys.exit(1)

        self.ahomo, self.alumo = self.g09_get_homo_lumo(self.monomer_a_fchk_lines)
        self.bhomo, self.blumo = self.g09_get_homo_lumo(self.monomer_b_fchk_lines)

        self.output = self.calculate_coupling()
        # when using default
        self.echh = self.output[0][4]
        self.echh_eff = self.output[0][8]
        self.ecll = self.output[3][4]
        self.ecll_eff = self.output[3][8]

    def calculate_coupling(self):

        # Get dimer values
        nbf_d = self.g09_get_basis_functions(self.dimer_fchk_lines)
        e_d = self.g09_get_alpha_orbital_energies(self.dimer_fchk_lines, nbf_d)
        c_d = self.g09_get_mo_coefficients(self.dimer_fchk_lines, nbf_d)
        s_d = self.get_overlap(self.sysname + '_A_B.log', nbf_d)
        f_d = self.g09_calculate_fock(c_d, s_d, nbf_d, e_d)

        # Get monomer A values
        nbf_monomer_a = self.g09_get_basis_functions(self.monomer_a_fchk_lines)
        homo_monomer_a, lumo_monomer_a = self.g09_get_homo_lumo(self.monomer_a_fchk_lines)
        homo_monomer_a += self.homo1
        lumo_monomer_a += self.lumo1
        c_monomer_a = self.g09_get_mo_coefficients(self.monomer_a_fchk_lines, nbf_monomer_a)

        # Get monomer B values
        nbf_monomer_b = self.g09_get_basis_functions(self.monomer_b_fchk_lines)
        homo_monomer_b, lumo_monomer_b = self.g09_get_homo_lumo(self.monomer_b_fchk_lines)
        homo_monomer_b += self.homo2
        lumo_monomer_b += self.lumo2
        c_monomer_b = self.g09_get_mo_coefficients(self.monomer_b_fchk_lines, nbf_monomer_b)

        # Calculate Coupling
        cfrag = self.g09_build_cfrag(nbf_d, nbf_monomer_a, nbf_monomer_b, c_monomer_a, c_monomer_b)
        out = self.g09_calculate_coupling(homo_monomer_a,
                                          lumo_monomer_a, homo_monomer_b, lumo_monomer_b, nbf_monomer_a,
                                          f_d, s_d, cfrag)
        return out

    @staticmethod
    @jit
    def g09_calculate_coupling(homo_monomer_a, lumo_monomer_a, homo_monomer_b, lumo_monomer_b, nbf_monomer_a,
                               f_d, s_d, cfrag):
        #            % ("MO 1", "MO 2", "e 1, eV", "e 2, eV", "J_12, meV", "S_12", "e^eff_1, eV", "e^eff_2, eV",
        #               "J^eff_12, meV",
        #               "J^eff_12, cm-1", "dE_12, meV"))
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
                out.append((i + 1, j + 1 - nbf_monomer_a, e1, e2, j12 * 1000.0, s12, ee1, ee2,
                            je12 * 1000.0, je12 * 8065.544, de * 1000.0))
        return out

    @staticmethod
    def g09_get_basis_functions(fchk_lines):
        # Get number of basis functions
        for line in fchk_lines:
            if "Number of basis functions" in line:
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

    def get_overlap(self, g09file, nbf):
        """
        Extracts the overlap matrix from a Gaussian logfile.
        Returns a numpy matrix.
        """
        with open(g09file, 'r') as f:
            data = f.read()
        overlap_matrix = np.zeros((nbf, nbf))
        # grab all text between  "Overlap ***" and "*** Kinetic"
        raw_overlap_string = re.findall(r'Overlap \*\*\*(.*?)\*\*\* Kinetic', data, re.DOTALL)
        raw_overlap_string = raw_overlap_string[0].replace('D', 'E')
        raw_overlap_elements = raw_overlap_string.split()
        matrix_elements = []
        for overlap_value in raw_overlap_elements:
            if 'E' in overlap_value:
                matrix_elements.append(overlap_value)
        overlap = self.create_matrix(overlap_matrix, matrix_elements, nbf)
        overlap = self.make_symmetric(overlap)
        overlap = np.matrix(overlap)

        return overlap

    def create_matrix(self, matrix, elements, nbf):
        """ create lower triangular matrix from list of matrix elements
        indexed like so:
            [[0,0,0,...,0],
            [1,2,0,...,0],
            [3,4,5,...,0]]
        nbf is number of basis functions
        elements is a list of matrix elements indexed like above, e.g.
            [0,1,2,3,...]
        Gaussian prints every 5 columns, so the mod5 accounts for this
        """
        count = 0  # count is our index
        # fill the main block, leaving remainder for triangle fill
        for i in range(0, nbf - nbf % 5, 5):
            matrix, count = self.triangle_fill(matrix, i, i + 5, i, count, elements)
            matrix, count = self.block_fill(matrix, i + 5, nbf, i, i + 5, count, elements)
        # finish filling the last triangle bit
        matrix, count = self.triangle_fill(matrix, nbf - nbf % 5, nbf, nbf - nbf % 5, count, elements)

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
