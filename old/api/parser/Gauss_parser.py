# import re
# import sys
# from collections import OrderedDict
#
# GaussPatterns = OrderedDict(
#     start=re.compile(r" \(Enter \S+l101\.exe\)"), route=re.compile(r" #[pPnNtT]*.*"),
#     link0=re.compile(r"^\s(%.+)\s*=\s*(.+)"),
#     charge_mul=re.compile(r"Charge\s+=\s*([-\d]+)\s+Multiplicity\s+=\s*(\d+)"),
#     num_basis_func=re.compile(r"([0-9]+)\s+basis functions"),
#     num_elec=re.compile(r"(\d+)\s+alpha electrons\s+(\d+)\s+beta electrons"),
#     pcm=re.compile(r"Polarizable Continuum Model"), stat_type=re.compile(r"imaginary frequencies"),
#     scf=re.compile(r"E\(.*\)\s*=\s*([-.\d]+)\s+"), mp2=re.compile(r"EUMP2\s*=\s*(.*)"),
#     oniom=re.compile(r"ONIOM:\s+extrapolated energy\s*=\s*(.*)"), termination=re.compile(r"(Normal|Error) termination"),
#     error=re.compile(r"(! Non-Optimized Parameters !|Convergence failure)"),
#     mulliken=re.compile(r"^\s*(Mulliken charges|Mulliken atomic charges)"),
#     mulliken_charge=re.compile(r'^\s+(\d+)\s+([A-Z][a-z]?)\s*(\S*)'),
#     end_mulliken=re.compile(r'(Sum of Mulliken )(.*)(charges)\s*=\s*(\D)'),
#     std_orientation=re.compile(r"Standard orientation"), input_orientation=re.compile(r"Input orientation"),
#     end=re.compile(r"--+"), orbital=re.compile(r"(Alpha|Beta)\s*\S+\s*eigenvalues --(.*)"),
#     thermo=re.compile(r"(Zero-point|Thermal) correction(.*)=\s+([\d.-]+)"),
#     forces_on=re.compile(r"Center\s+Atomic\s+Forces\s+\(Hartrees/Bohr\)"),
#     forces_off=re.compile(r"Cartesian\s+Forces:\s+Max.*RMS.*"),
#     forces=re.compile(r"\s+(\d+)\s+(\d+)\s+([0-9.-]+)\s+([0-9.-]+)\s+([0-9.-]+)"),
#     freq_on=re.compile(r"Harmonic\sfrequencies\s+\(cm\*\*-1\),\sIR\sintensities.*Raman.*"),
#     freq=re.compile(r"Frequencies\s--\s+(.*)"),
#     normal_mode=re.compile(r"\s+(\d+)\s+(\d+)\s+([0-9.-]{4,5})\s+([0-9.-]{4,5}).*"),
#     mo_coeff=re.compile(r"Molecular Orbital Coefficients:"),
#     mo_coeff_name=re.compile(r"\d+\s((\d+|\s+)\s+([a-zA-Z]{1,2}|\s+))\s+(\d+\S+)"),
#     hessian=re.compile(r"Force constants in Cartesian coordinates:"), resume=re.compile(r"^\s1\\1\\GINC-\S*"),
#     resume_end=re.compile(r"^\s.*\\\\@"), bond_order=re.compile(r"Wiberg bond index matrix in the NAO basis:"),
# )
#
#
# class Logparser:
#
#     def __init__(self, filename):
#         self.filename = filename
#         if self.filename[-4:] != '.log':
#             sys.exit('trying to parse a file with wrong extension!')
#         self.data = OrderedDict()
#
#     def parse(self):
#         with open(self.filename, 'r') as f:
#             for line in f:

