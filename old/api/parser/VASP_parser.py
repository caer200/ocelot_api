from api.parser.Parser import Parser
from api.routines.geometry import frac2cart
import numpy as np


class VASP_parser(Parser):

    def __init__(self, filename):
        super().__init__(filename)
        if 'POSCAR' in self.filename.upper() or 'CONTCAR' in self.filename.upper():
            self.type = 'poscar'
        if 'OUTCAR' in self.filename.upper():
            self.type = 'outcar'
        if 'KPOINTS' in self.filename.upper():
            self.type = 'kpoints'
        if 'INCAR' in self.filename.upper():
            self.type = 'incar'
        if 'POTCAR' in self.filename.upper():
            self.type = 'potcar'
        if 'OSZICAR' in self.filename.upper():
            self.type = 'oszicar'
        if 'EIGENVAL' in self.filename.upper():
            self.type = 'eigenval'

    def to_cart(self, output):
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        system = lines[7][0].upper()
        if system == 'C' or system == 'K':
            with open(output, 'w') as f:
                f.write(''.join(lines))
        else:
            pbc = np.zeros((3, 3))
            for i in range(3):
                pbc[i] = [float(f) for f in lines[2+i].strip().split()]
            headerlines = lines[:8]
            fraclines = [l for l in lines[8:] if len(l.strip().split()) > 0]
            cartlines = []
            for i in range(len(fraclines)):
                fraccoords = frac2cart([float(f) for f in fraclines[i].strip().split()[:3]], pbc)
                cartlines += '{:.6f} {:.6f} {:.6f} \n'.format(*fraccoords)
            with open(output, 'w') as f:
                f.write(''.join(headerlines + cartlines))


