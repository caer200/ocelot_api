from pymatgen.io.vasp import Poscar, PotcarSingle, Potcar
import numpy as np
import math
import glob
import sys
import os
import warnings
import multiprocessing


####################

incar_single_point_scf = {  # vasp single point scf
    'SYSTEM': 'single point scf',
    'PREC': 'Accurate',
    # electronic
    'ENCUT': '0',
    'NELM': '200',
    'EDIFF': '1e-06',
    'ISMEAR': '0',
    'SIGMA': '0.05',
    'ISPIN': '1',
    'LORBIT': '0',
    'LREAL': 'Auto',
    'LASPH': 'True',
    'ALGO': 'Normal',
    # relaxation
    'NSW': '1',
    'IBRION': '2',
    'ISIF': '2',
    'EDIFFG': '0.01',
    'IVDW': '12',
    # file reading
    'ISTART': '0',
    'ICHARG': '0',
    'LCHARGE': 'True',
    'LWAVE': 'False',
    # parallel
    'NCORE': str(multiprocessing.cpu_count()),
}

incar_nscf = {  # vasp nscf band structure
    'SYSTEM': 'nscf',
    'PREC': 'Accurate',
    # electronic
    'ENCUT': '0',
    'NELM': '200',
    'EDIFF': '1e-06',
    'ISMEAR': '0',
    'SIGMA': '0.05',
    'ISPIN': '1',
    'LORBIT': '0',
    'LREAL': 'Auto',
    'LASPH': 'True',
    'ALGO': 'Normal',
    # relaxation
    'NSW': '1',
    'IBRION': '2',
    'ISIF': '2',
    'EDIFFG': '0.01',
    'IVDW': '12',
    # file reading
    'ISTART': '0',
    'ICHARG': '11',
    'LCHARGE': 'False',
    'LWAVE': 'False',
    # parallel
    'NCORE': str(multiprocessing.cpu_count()),
}

incar_full_relaxation = {
    'SYSTEM': 'full relaxation',
    'PREC': 'Accurate',
    # electronic
    'ENCUT': '0',
    'NELM': '200',
    'EDIFF': '1e-06',
    'ISMEAR': '0',
    'SIGMA': '0.05',
    'ISPIN': '1',
    'LORBIT': '0',
    'LREAL': 'Auto',
    'LASPH': 'True',
    'ALGO': 'Normal',
    # relaxation
    'NSW': '500',
    'IBRION': '2',
    'ISIF': '3',
    'EDIFFG': '0.01',
    'IVDW': '12',
    # file reading
    'ISTART': '0',
    'ICHARG': '0',
    'LCHARGE': 'False',
    'LWAVE': 'False',
    # parallel
    'NCORE': str(multiprocessing.cpu_count()),
}

incar_fix_box_relaxation = {
    'SYSTEM': 'fix box relaxation',
    'PREC': 'Accurate',
    # electronic
    'ENCUT': '0',
    'NELM': '200',
    'EDIFF': '1e-06',
    'ISMEAR': '0',
    'SIGMA': '0.05',
    'ISPIN': '1',
    'LORBIT': '0',
    'LREAL': 'Auto',
    'LASPH': 'True',
    'ALGO': 'Normal',
    # relaxation
    'NSW': '500',
    'IBRION': '2',
    'ISIF': '2',
    'EDIFFG': '0.01',
    'IVDW': '12',
    # file reading
    'ISTART': '0',
    'ICHARG': '0',
    'LCHARGE': 'False',
    'LWAVE': 'False',
    # parallel
    'NCORE': str(multiprocessing.cpu_count()),
}

kpoint_auto = {
    'SYSTEM': 'auto mesh',
    'nkpts': '0',
    'type': 'gamma',
    'grid': '0 0 0',
    'shift': '0 0 0'
}

template = dict(
    incar_single_point_scf=incar_single_point_scf,
    incar_fix_box_relaxation=incar_fix_box_relaxation,
    incar_full_relaxation=incar_full_relaxation,
    incar_nscf=incar_nscf,
    kpoint_auto=kpoint_auto,
)


class WriteInp:

    def __init__(self, jobtype, pawpath):
        """

        :param jobtype: fbrlx, arlx, scf, nscf
        :param pawpath:
        """
        if len(glob.glob('POSCAR')) != 1:
            warnings.warn('failed in writing inp for vasp! there is no POSCAR file at ' + os.getcwd())
            sys.exit(1)
        self.jobtype = jobtype
        self.pawpath = pawpath
        self.enmax = self.write_potcar(pawpath, 'POSCAR')
        self.write_kpoints()
        self.write_incar(self.jobtype, self.enmax)

    @staticmethod
    def write_potcar(pawpath, poscarfilename):
        """
        Function to create molecule POTCAR from individual element POTCARs.
        This allow ENMAX for a given atom to be extracted.
        :param pawpath: abs path to paw
        :param poscarfilename: Output molecule poscar filename.
        :return: ENMAX of POTCAR
        """
        poscar = Poscar.from_file(poscarfilename, check_for_POTCAR=False)  # Read POSCAR from file.
        potcar_list = []  # Initialize list of element POTCARs
        enmax_list = []  # Initialize list of element ENMAXs
        elements = poscar.site_symbols  # Extract elements from POSCAR
        for i in elements:
            single_potcar = PotcarSingle.from_file(pawpath + "/" + i + "/POTCAR")  # Extract element POTCARs
            potcar_list.append(single_potcar.data)
            enmax_list.append(single_potcar.ENMAX)  # Extract ENMAX for pseudopotential
        potcar_dictionary = dict(zip(elements, potcar_list))
        potcar = Potcar(symbols=elements, functional='PBE',
                        sym_potcar_map=potcar_dictionary)  # Create POTCAR for POSCAR
        potcar.write_file('POTCAR')  # Write POTCAR to given file
        enmax_max = max(enmax_list)  # Determine ENMAX
        return enmax_max

    @staticmethod
    def write_kpoints(density=20, temp_dict=template['kpoint_auto']):
        if len(glob.glob('*kpath*')) != 1:
            with open('POSCAR', 'r') as p:
                lines = p.readlines()[2:5]
            la = np.array(lines[0].strip().split(), dtype='float')
            lb = np.array(lines[1].strip().split(), dtype='float')
            lc = np.array(lines[2].strip().split(), dtype='float')
            a = np.linalg.norm(la)
            b = np.linalg.norm(lb)
            c = np.linalg.norm(lc)
            ga = math.ceil(density/a)
            gb = math.ceil(density/b)
            gc = math.ceil(density/c)
            temp_dict['grid'] = ' '.join([str(ga), str(gb), str(gc)])
            WriteInp.write_fdict(temp_dict, 'KPOINTS')
        else:
            os.system('mv ' + glob.glob('*kpath*')[0] + ' KPOINTS')

    @staticmethod
    def write_fdict(fdict, fname):
        if fname == 'INCAR':
            with open(fname, 'w') as f:
                for i in fdict.keys():
                    f.write(i + ' = ' + str(fdict[i]) + '\n')
        elif fname == 'KPOINTS':
            with open(fname, 'w') as f:
                for i in fdict.keys():
                    f.write(str(fdict[i]) + '\n')

    @staticmethod
    def write_incar(jtype, enmax):
        if jtype == 'fbrlx':
            incar_d = template['incar_fix_box_relaxation']
            incar_d['ENCUT'] = enmax
            WriteInp.write_fdict(incar_d, 'INCAR')
        elif jtype == 'arlx':
            incar_d = template['incar_full_relaxation']
            incar_d['ENCUT'] = enmax * 1.3
            WriteInp.write_fdict(incar_d, 'INCAR')
        elif jtype == 'scf':
            incar_d = template['incar_single_point_scf']
            incar_d['ENCUT'] = enmax
            WriteInp.write_fdict(incar_d, 'INCAR')
        elif jtype == 'nscf':
            incar_d = template['incar_nscf']
            incar_d['ENCUT'] = enmax
            WriteInp.write_fdict(incar_d, 'INCAR')
        else:
            warnings.warn('jobtype not recognized')
            sys.exit(1)
