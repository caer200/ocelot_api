class Element:
    covalent_radii = {'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.73, 'N': 0.71, 'O': 0.66,
                      'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
                      'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
                      'Mn': 1.50, 'Fe': 1.42, 'Co': 1.38, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20,
                      'As': 1.19, 'Se': 1.20, 'Br': 1.20, 'Kr': 1.16, 'Rb': 2.20, 'Sr': 1.95, 'Y': 1.90, 'Zr': 1.75,
                      'Nb': 1.64, 'Mo': 1.54, 'Tc': 1.47, 'Ru': 1.46, 'Rh': 1.42, 'Pd': 1.39, 'Ag': 1.45, 'Cd': 1.44,
                      'In': 1.42, 'Sn': 1.39, 'Sb': 1.39, 'Te': 1.38, 'I': 1.39, 'Xe': 1.40, 'Cs': 2.44, 'Ba': 2.15,
                      'La': 2.07, 'Ce': 2.04, 'Pr': 2.03, 'Nd': 2.01, 'Pm': 1.99, 'Sm': 1.98, 'Eu': 1.98, 'Gd': 1.96,
                      'Tb': 1.94, 'Dy': 1.92, 'Ho': 1.92, 'Er': 1.89, 'Tm': 1.90, 'Yb': 1.87, 'Lu': 1.87, 'Hf': 1.75,
                      'Ta': 1.70, 'W': 1.62, 'Re': 1.51, 'Os': 1.44, 'Ir': 1.41, 'Pt': 1.36, 'Au': 1.36, 'Hg': 1.32,
                      'Tl': 1.45, 'Pb': 1.46, 'Bi': 1.48, 'Po': 1.40, 'At': 1.50, 'Rn': 1.50, 'Fr': 2.60, 'Ra': 2.21,
                      'Ac': 2.15, 'Th': 2.06, 'Pa': 2.00, 'U': 1.96, 'Np': 1.90, 'Pu': 1.87, 'Am': 1.80, 'Cm': 1.69}
    vdw_radii = {'H': 1.2, 'He': 1.4, 'Li': 1.82, 'Be': 1.53, 'B': 1.92, 'C': 1.7, 'N': 1.55, 'O': 1.52, 'F': 1.47,
                 'Ne': 1.54, 'Na': 2.27, 'Mg': 1.73, 'Al': 1.84, 'Si': 2.1, 'P': 1.8, 'S': 1.8, 'Cl': 1.75, 'Ar': 1.88,
                 'K': 2.75, 'Ca': 2.31, 'Sc': 2.11, 'Ti': None, 'V': None, 'Cr': None, 'Mn': None, 'Fe': None,
                 'Co': None, 'Ni': 1.63, 'Cu': 1.4, 'Zn': 1.39, 'Ga': 1.87, 'Ge': 2.11, 'As': 1.85, 'Se': 1.9,
                 'Br': 1.85, 'Kr': 0.88, 'Rb': 3.03, 'Sr': 2.49, 'Y': None, 'Zr': None, 'Nb': None, 'Mo': None,
                 'Tc': None, 'Ru': None, 'Rh': None, 'Pd': 1.63, 'Ag': 1.72, 'Cd': 1.58, 'In': 1.93, 'Sn': 2.17,
                 'Sb': 2.06, 'Te': 2.06, 'I': 1.98, 'Xe': 1.08, 'Cs': 3.43, 'Ba': 2.68, 'La': None, 'Ce': None,
                 'Pr': None, 'Nd': None, 'Pm': None, 'Sm': None, 'Eu': None, 'Gd': None, 'Tb': None, 'Dy': None,
                 'Ho': None, 'Er': None, 'Tm': None, 'Yb': None, 'Lu': None, 'Hf': None, 'Ta': None, 'W': None,
                 'Re': None, 'Os': None, 'Ir': None, 'Pt': 1.75, 'Au': 1.66, 'Hg': 1.55, 'Tl': 1.96, 'Pb': 2.02,
                 'Bi': 2.07, 'Po': 1.97, 'At': 1.27, 'Rn': 1.2, 'Fr': None, 'Ra': None, 'Ac': None, 'Th': None,
                 'Pa': None, 'U': None, 'Np': None, 'Pu': None, 'Am': None, 'Cm': None, 'Bk': None, 'Cf': None,
                 'Es': None, 'Fm': None, 'Md': None, 'No': None, 'Lr': None, 'Rf': None, 'Db': None, 'Sg': None,
                 'Bh': None, 'Hs': None, 'Mt': None, 'Ds': None, 'Rg': None, 'Cn': None, 'Nh': None, 'Fl': None,
                 'Mc': None, 'Lv': None, 'Ts': None, 'Og': None, }
    atomic_numbers = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                      'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21,
                      'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
                      'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41,
                      'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
                      'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60,
                      'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69,
                      'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
                      'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88,
                      'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
                      'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107,
                      'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
                      'Ts': 117, 'Og': 118, }

    # max bonds that an atom can have, a short dict should be enough for organics...
    valence_dict = {'O': 2, 'S': 2, 'N': 3, 'P': 3, 'H': 1, 'C': 4, 'F': 1, 'Cl': 1, 'Br': 1, }

    def __init__(self, name):
        self.name = name

    def __eq__(self, other):
        if self.name == other.name:
            return 1
        return 0

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.name)

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    @property
    def covrad(self):
        try:
            return self.covalent_radii[self.name]
        except KeyError:
            return None

    @property
    def vanrad(self):
        try:
            return self.vdw_radii[self.name]
        except KeyError:
            return None

    @property
    def atomic_number(self):
        try:
            return self.atomic_numbers[self.name]
        except KeyError:
            return None

    @property
    def valence(self):
        try:
            return self.valence_dict[self.name]
        except KeyError:
            return None
