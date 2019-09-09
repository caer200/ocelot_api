import sys
import random
import os
import re
import subprocess
import argparse
from scipy import optimize
import shutil
import time
from pymatgen.io.gaussian import GaussianInput
from collections import OrderedDict
import os

"""
originally written by Sean M. Ryno, Cheng Zhong, Haitao Sun, see /legacy/aw_tuning.py
for a certain mol, get tuned w

default keywords for gaussian
'scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)'

dev-ing
"""

# SUGGESTED_route_parameters = {
#     'scf': {
#         'xqc': '',
#         'fermi': '',
#         'noincfock': '',
#         'novaracc': '',
#         'ndamp': '35',
#         'conver': '6',
#         'vshift': '500'
#     }
# }
# SUGGESTED_route_parameters = OrderedDict(SUGGESTED_route_parameters)
#
# def extractEHL(logfilename):
#     # Extract required information from log files
#     with open(logfilename, 'r') as f:
#         LogFile = f.readlines()
#     normal = 0
#     for line in LogFile:
#         if re.search('Normal termination', line):
#             normal = 1
#         if re.search('SCF Done', line):
#             Energy = line.split()[4]
#         if re.search('Alpha  occ. eigenvalue', line):
#             AlphaHOMO = line.split()[-1]
#             occFlag = 1
#         if re.search('Beta  occ. eigenvalues', line):
#             BetaHOMO = line.split()[-1]
#             occFlag = 1
#         if re.search('Alpha virt. eigenvalues', line) and (occFlag == 1):
#             AlphaLUMO = line.split()[4]
#             occFlag = 0
#         if re.search('Beta virt. eigenvalues', line) and (occFlag == 1):
#             BetaLUMO = line.split()[4]
#             occFlag = 0
#
#     if normal == 0:
#         return None
#     try:
#         BetaHOMO
#     except NameError:
#         BetaHOMO = AlphaHOMO
#     try:
#         BetaLUMO
#     except NameError:
#         BetaLUMO = AlphaLUMO
#
#     if float(BetaHOMO) > float(AlphaHOMO):
#         HOMO = float(BetaHOMO)
#     else:
#         HOMO = float(AlphaHOMO)
#     if float(BetaLUMO) < float(AlphaLUMO):
#         LUMO = float(BetaLUMO)
#     else:
#         LUMO = float(AlphaLUMO)
#
#     return [float(Energy) * 27.211396132, HOMO * 27.211396132, LUMO * 27.211396132]
#
#
# class GaussTuning:
#
#     def __init__(
#             self, sysname, mol, wdir, gaubin, neutral_charge=0, neutralspin=1, tuning_what='w', tuning_scheme='gap',
#             functional='lc-whpbe', basis_set='6-311g(d,p)', route_parameters=SUGGESTED_route_parameters,
#             link0_parameters=None, input_parameters=None, w_max=0.5, w_min=0.05, conver=4,
#     ):
#         """
#         :param mol: pmg mol
#         :param wdir: absolute path for working dir
#         :param tuning_what: 'w' 'a'
#         :param tuning_scheme: 'gap' 'jh' 'jl'
#         """
#         self.sysname = sysname
#         self.gaubin = gaubin
#
#         self.mol = mol
#         self.tuning_what = tuning_what
#         self.tuning_scheme = tuning_scheme
#         self.wdir = wdir
#         self.cycle = 0
#         self.conver = conver
#         self.omega = random.uniform(0, 1)
#         self.n_charge = neutral_charge
#         self.n_spin = neutralspin
#
#         self.functional = functional
#         self.basis_set = basis_set
#         self.route_parameters = route_parameters
#         self.link0_parameters = link0_parameters
#         self.input_parameters = input_parameters
#
#         self.w_max = w_max
#         self.w_min = w_min
#
#         self.c_charge = int(self.n_charge) + 1
#         self.a_charge = int(self.n_charge) - 1
#
#
#         if self.n_spin == 1:
#             self.c_spin = 2
#             self.a_spin = 2
#         else:
#         # elif self.n_spin > 1:
#             self.c_spin = self.n_spin - 1
#             self.a_spin = self.n_spin - 1
#
#         self.logstring = '# n_charge={} c_charge={} a_charge={} \n'.format(self.n_charge, self.c_charge, self.a_charge)
#         self.logstring += '# n_spin={} c_spin={} a_spin={} \n'.format(self.n_spin, self.c_spin, self.a_spin)
#         self.logstring += '# current C\t '
#
#     @property
#     def omega_iopstring(self):
#         return str(int(float(self.omega) * 10000)).zfill(5) + '00000)'
#
#     def getc(self, aorw):
#         """
#         write
#         :return:
#         """
#         self.omega = round(self.omega, self.conver)
#         self.route_parameters['iop'] = {
#             '3/107': self.omega_iopstring,
#             '3/108': self.omega_iopstring
#         }
#         n_input_fn = 'n_omega_{}.com'.format(self.omega_iopstring)
#         a_input_fn = 'a_omega_{}.com'.format(self.omega_iopstring)
#         c_input_fn = 'c_omega_{}.com'.format(self.omega_iopstring)
#
#
#
#
#
#
#
#
#     def run(self, guassbin, functional, basis_set, route_parameters, link0_parameters, input_parameters=None):
#         os.makedirs(self.wdir, exist_ok=True)
#         os.chdir(self.wdir)
#         iiter = 0
#         gin = GaussianInput(
#             self.mol, self.n_charge, self.n_spin, self.sysname, functional='lc-whpbe', basis_set='6-31g(d)',
#             route_parameters=route_parameters, input_parameters=input_parameters, link0_parameters=None,
#         )
#         gin.write_file(, cart_coords=True)
#
#     def runGaussian(self):
#         print(' ', file=sys.stdout)
#
#         # Write Alpha and Omega values and perform a sanity check for the values.
#         tuningInput = ''
#         tuningNameSuffix = ''
#         if self.omega != '-----':
#             if (float(self.omega) < 10) and (float(self.omega) >= 0):
#                 self.omega = round(self.omega, self.conver)
#                 tuningInput = tuningInput + 'IOP(3/107=' + str(int(float(self.omega) * 10000)).zfill(
#                     5) + '00000) IOP(3/108=' + str(int(float(self.omega) * 10000)).zfill(5) + '00000) '
#             else:
#                 print("Omega", self.omega, "must be a number between 0 and 10", file=sys.stderr)
#                 print("Exiting...", file=sys.stderr)
#                 sys.exit(1)
#             tuningNameSuffix = tuningNameSuffix + 'w' + str(int(float(self.omega) * 10000)).zfill(5)
#
#         # Move old checkpoint files and set new names
#         try:
#             oldNF, oldCF, oldAF = self.neutralFilename, self.cationFilename, self.anionFilename
#         except:
#             oldNF, oldCF, oldAF = '', '', ''
#
#         self.neutralFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
#                                                                                                     '_') + '_' + tuningNameSuffix + '_n'
#         self.cationFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
#                                                                                                    '_') + '_' + tuningNameSuffix + '_c'
#         self.anionFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
#                                                                                                   '_') + '_' + tuningNameSuffix + '_a'
#         self.tuningNameSuffix = tuningNameSuffix
#
#         # Search the InfoTable, if data is already there extract it from the table
#         # If alpha and omega are the same as the previous cycle do not increase cycle number
#         try:
#             infoTable = open(self.infoTable, 'r')
#             cycle = self.cycle
#             lastLine = 1
#             for line in infoTable:
#                 if self.cycle == cycle - 1:
#                     lastLine = 0
#                     break
#                 if (str(self.alpha) == line.split()[0]) and (str(self.omega) == line.split()[1]):
#                     print("alpha:", str(self.alpha), ' omega:', str(self.omega), ' found in ', self.infoTable, sep='',
#                           file=sys.stdout)
#                     self.C, self.J2, self.Jn2, self.O2, self.gap, self.IP, self.EA, self.neutralE, self.N_HOMO, self.N_LUMO, self.cationE, self.C_HOMO, self.C_LUMO, self.anionE, self.A_HOMO, self.A_LUMO = [
#                         float(i) for i in line.split()[2:]]
#                     cycle = cycle + 1
#             if lastLine == 0:
#                 self.cycle = cycle
#                 print("-----cycle", str(self.cycle).zfill(2), "-----", tuningNameSuffix, "  ", self.criterion, ":",
#                       str(self.C), "  gap:", str(self.gap), sep='', file=sys.stdout)
#                 return self.C
#             if self.cycle == cycle - 1:
#                 return self.C
#             infoTable.close()
#         except IOError:
#             pass
#
#         # Create a working directory to perform the tuning in
#         if not os.path.exists(self.workingDIR):
#             os.makedirs(self.workingDIR)
#         os.chdir(self.workingDIR)
#
#         # Create Gaussian input file
#         self.genGJF(tuningInput, tuningNameSuffix)
#
#         # Call Gaussian
#         if len(self.extractEHL(self.neutralFilename)) != 3:
#             print("Running: ", self.neutralFilename, ".gjf...",
#                   time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
#             if os.path.isfile(oldNF + '.chk'):
#                 shutil.copy(oldNF + '.chk', self.neutralFilename + '.chk')
#             else:
#                 self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
#             runNeutral = subprocess.Popen([self.program, self.neutralFilename + '.gjf'])
#             runNeutral.wait()
#         else:
#             print("logfile ", self.neutralFilename, ".log found, directly extracting information.", sep='',
#                   file=sys.stdout)
#         try:
#             neutralE, N_HOMO, N_LUMO = self.extractEHL(self.neutralFilename)
#         except ValueError:
#             print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.neutralFilename,
#                   ".log. Please check this files.", sep='', file=sys.stderr)
#             print(
#                 "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
#                 file=sys.stderr)
#             sys.exit(1)
#         if self.criterion != "Jl":
#             if len(self.extractEHL(self.cationFilename)) != 3:
#                 print("Running ", self.cationFilename, ".gjf...",
#                       time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
#                 if os.path.isfile(oldCF + ".chk"):
#                     shutil.copy(oldCF + ".chk", self.cationFilename + ".chk")
#                 else:
#                     self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
#                 runCation = subprocess.Popen([self.program, self.cationFilename + ".gjf"])
#                 runCation.wait()
#             else:
#                 print("logfile ", self.cationFilename, ".log found, directly extracting information.", sep='',
#                       file=sys.stdout)
#             try:
#                 cationE, C_HOMO, C_LUMO = self.extractEHL(self.cationFilename)
#             except ValueError:
#                 print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.cationFilename,
#                       ".log. Please check this files.", sep='', file=sys.stderr)
#                 print(
#                     "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
#                     file=sys.stderr)
#                 sys.exit(1)
#         else:
#             cationE, C_HOMO, C_LUMO = 0, 0, 0
#         if self.criterion != "Jh":
#             if len(self.extractEHL(self.anionFilename)) != 3:
#                 print("Running ", self.anionFilename, ".gjf...",
#                       time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
#                 if os.path.isfile(oldAF + ".chk"):
#                     shutil.copy(oldAF + ".chk", self.anionFilename + ".chk")
#                 else:
#                     self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
#                 runAnion = subprocess.Popen([self.program, self.anionFilename + ".gjf"])
#                 runAnion.wait()
#             else:
#                 print("logfile ", self.anionFilename, ".log found, directly extracting information.", sep='',
#                       file=sys.stdout)
#             try:
#                 anionE, A_HOMO, A_LUMO = self.extractEHL(self.anionFilename)
#             except ValueError:
#                 print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.anionFilename,
#                       ".log. Please check this files.", sep='', file=sys.stderr)
#                 print(
#                     "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
#                     file=sys.stderr)
#                 sys.exit(1)
#         else:
#             anionE, A_HOMO, A_LUMO = 0, 0, 0
#
#         os.chdir(self.currentDIR)
#         if (self.criterion != "Jh") and (self.criterion != "Jl"):
#             IP = neutralE - cationE
#             EA = neutralE - anionE
#             J2 = (N_HOMO - IP) ** 2 + (A_HOMO + EA) ** 2
#             J = J2 ** 0.5
#             gap = N_LUMO - N_HOMO
#             Jh = abs(N_HOMO - IP)
#             Jl = abs(N_LUMO + EA)
#             Jn2 = Jh ** 2 + Jl ** 2
#             Jn = Jn2 ** 0.5
#             O2 = (N_HOMO - C_LUMO) ** 2 + (N_LUMO - A_HOMO) ** 2
#             C = eval(self.criterion)
#         else:
#             IP = 0
#             EA = 0
#             J2 = 0
#             J = J2 ** 0.5
#             Jl = 0
#             Jn2 = 0
#             Jn = Jn2 ** 0.5
#             O2 = 0
#         if self.criterion == "Jh":
#             IP = neutralE - cationE
#             gap = N_LUMO - N_HOMO
#             Jh = abs(N_HOMO - IP)
#             C = eval(self.criterion)
#         if self.criterion == "Jl":
#             EA = neutralE - anionE
#             gap = N_LUMO - N_HOMO
#             Jl = abs(N_LUMO + EA)
#             C = eval(self.criterion)
#
#         # Assign self values based on new run
#         self.C, self.J, self.J2, self.Jn, self.Jn2, self.O2, self.gap, self.IP, self.EA, self.neutralE, self.N_HOMO, self.N_LUMO, self.cationE, self.C_HOMO, self.C_LUMO, self.anionE, self.A_HOMO, self.A_LUMO = C, J, J2, Jn, Jn2, O2, gap, IP, EA, neutralE, N_HOMO, N_LUMO, cationE, C_HOMO, C_LUMO, anionE, A_HOMO, A_LUMO
#         lenNE = str(len('{:.6f}'.format(self.neutralE)) + 2)
#         lenCE = str(len('{:.6f}'.format(self.cationE)) + 2)
#         lenAE = str(len('{:.6f}'.format(self.anionE)) + 2)
#
#         # Make InfoTable
#         try:
#             infoTable = open(self.infoTable, 'r')
#         except IOError:
#             infoTable = open(self.infoTable, 'w')
#             infoTable.write((
#                     '{:<10}{:<10}{:<11}{:<11}{:<11}{:<11}{:<11}{:<11}{:<11}{:<' + lenNE + '}{:<11}{:<11}{:<' + lenCE + '}{:<11}{:<11}{:<' + lenAE + '}{:<11}{:<11}\n').format(
#                 'Alpha', 'Omega', 'opt:' + self.criterion, 'J2', 'Jn2', 'O2', 'gap', 'IP', 'EA', 'neutralE', 'N_HOMO',
#                 'N_LUMO', 'cationE', 'C_HOMO', 'C_LUMO', 'anionE', 'A_HOMO', 'A_LUMO'))
#         infoTable.close()
#         infoTable = open(self.infoTable, 'a')
#         infoTable.write((
#                 '{:<10}{:<10}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<' + lenNE + '.6f}{:<11.6f}{:<11.6f}{:<' + lenCE + '.6f}{:<11.6f}{:<11.6f}{:<' + lenAE + '.6f}{:<11.6f}{:<11.6f}\n').format(
#             self.alpha, self.omega, C, J2, Jn2, O2, gap, IP, EA, neutralE, N_HOMO, N_LUMO, cationE, C_HOMO, C_LUMO,
#             anionE, A_HOMO, A_LUMO))
#         infoTable.close()
#         infoTable = open(self.infoTable, 'r')
#         self.cycle = self.cycle + 1
#
#         print("-----cycle", str(self.cycle).zfill(2), "-----", tuningNameSuffix, "  ", self.criterion, ":", str(C),
#               "  gap:", str(gap), sep='', file=sys.stdout)
#
#         return C
#
#
# # Setup Parser and help files
# parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
#                                  description="A script used with Gaussian09/16 to auto optimize omega and alpha for long-range corrected DFT methodologies. This program would not be possible without the original version by Cheng Zhong.")
# parser.add_argument("-g", "--gauss", help="Set the executable name to be used. Default: g09", type=str, default='g09')
# parser.add_argument("-l", "--level",
#                     help="Set the theoretical level. Default: read from gjf file. If not found, use LC-wPBE/6-31g(d)",
#                     type=str, default='')
# parser.add_argument("-m", "--memory",
#                     help="Set the memory used by Gaussian (with unit).  Default: read from gjf file. If not found, use: 500MB",
#                     type=str, default='')
# parser.add_argument("-n", "--nproc",
#                     help="Set the number of cores used by Gaussian. Default: read from gjf file. If not found, use 1",
#                     type=str, default='')
# parser.add_argument("-p", "--precision",
#                     help="Set the precision of omega. Default: 4, means omega will converge to 0.0001 ", type=int,
#                     default='4')
# parser.add_argument("-P", "--alphaPrecision", help="Set the precision of alpha. Default: equal to precision of omega",
#                     type=int)
# parser.add_argument("-b", "--bound", help="Set the boundary for omega, e.g. -b 0.05,0.3. Default: 0.05,0.5", type=str,
#                     default='0.05,0.5')
# parser.add_argument("-B", "--alphaBound", help="Set the boundary for alpha, e.g. -B 0.05,0.3. Default: 0.0,0.5",
#                     type=str, default='0.0,0.5')
# parser.add_argument("-k", "--keywords",
#                     help="Set additional Gaussian keywords. Default: \"scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)\"",
#                     type=str, default='')
# parser.add_argument("-c", "--criterion",
#                     help="Set the optimization criterion (the value to minimize), available options are:\nJ2---((HOMO-IP)^2+(A_HOMO+EA)^2)\nJh---(HOMO-IP)\nJl(LUMO+EA)\nJn2---((HOMO-IP)^2+(LUMO+EA)^2)\nO2---((A_HOMO-LUMO)^2+(C_LUMO-HOMO)^2)\nor any other valid experssion, e.g. -c \"(A_HOMO-N_LUMO)**2\"\nIf multiple criterion is used, seperate them by comma. The order and number of criterions must be the same to that of structures. Default: J2",
#                     type=str, default='J2')
# parser.add_argument("-t", "--tuning", help="Set the tuning scheme. Available options are: w(omega),a(alpha),aw(alpha and omega),s(scan alpha or omege, must be used with -a or -w option)\nPut m in front of above options (mw,ma,maw,ms) means tuning with multiple structure,\n\
# e.g. \"-t mw\" means optimize the omega to minimize sum of the criterion of all the input structure.\nIf different criterion is used for different structure, use -c option to designate criterion explicitly for each structure, \ne.g. \"-t ma -c Jh,Jl\" means optimize alpha to minimize sum of Jh of structure 1 and Jl of structure 2.  Default: w",
#                     type=str, default='w')
# parser.add_argument("-a", "--alpha",
#                     help="Set the value of alpha. Can be a single value or multiple value seperated by space, e.g. -a '0.1 0.3 0.5 0.7' Default: None",
#                     type=str)
# parser.add_argument("-w", "--omega",
#                     help="Set the value of omega. Can be a single value or multiple value seperated by space. Default: None",
#                     type=str)
# parser.add_argument("gaussianInput", nargs='*')
# args = parser.parse_args()
#
# # Set values based on parser input
# program = args.gauss
# level = args.level
# memory = args.memory
# nproc = args.nproc
# bound = [float(i) for i in args.bound.split(',')]
# alphaBound = [float(i) for i in args.alphaBound.split(',')]
# keywords = args.keywords
# tuning = args.tuning
# criterions = args.criterion.split(',')
#
# # Do check to see if multiple files are being requested to run and run them
# # Tuning for single structures
# if 'm' not in args.tuning:
#     if len(criterions) > 1:
#         print("ERROR: Multiple criterion have been requested, but this must be used with the '-t m; option.",
#               file=sys.stderr)
#         print("Exiting...", file=sys.stderr)
#         sys.exit(1)
#     for g09input in args.gaussianInput:
#         tuningObj = tuningGaussian(g09input, nproc=nproc, mem=memory, level=level, addKeyWords=keywords,
#                                    criterion=args.criterion, program=program)
#         tuningObj.conver = args.precision
#         if not args.alphaPrecision:
#             args.alphaPrecision = args.precision
#         tuningObj.alphaConver = args.alphaPrecision
#
#         # Now run the requested tuning procedure and optimize as requested
#         # Omega tuning
#         if args.tuning == 'w':
#             def wtune(omega, wtune=tuningObj):
#                 wtune.omega = omega
#                 return wtune.runGaussian()
#
#
#             if args.alpha:
#                 for alpha in args.alpha.split():
#                     tuningObj.alpha = float(alpha)
#                     res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)
#             else:
#                 res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)
#
#         # Alpha tuning
#         if args.tuning == 'a':
#             def atune(alpha, atune=tuningObj):
#                 atune.alpha = alpha
#                 return atune.runGaussian()
#
#
#             if args.omega:
#                 for omega in args.omega.split():
#                     tuningObj.omega = float(omega)
#                     res = optimize.fminbound(atune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)
#             else:
#                 res = optimize.fminbound(atune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)
#
#         # Alpha and Omega tuning
#         if args.tuning == 'aw':
#             OUTCYCLE = 1
#
#
#             def awtune(alpha, awtune=tuningObj):
#                 global OUTCYCLE
#
#                 def wtune(omega, wtune=tuningObj):
#                     global OUTCYCLE
#                     wtune.omega = omega
#                     return wtune.runGaussian()
#
#                 awtune.alpha = alpha
#                 res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)
#                 awtune.omega = res
#                 print("=====OUTCYCLE", str(OUTCYCLE).zfill(2), "=====", awtune.tuningNameSuffix, "  ", awtune.criterion,
#                       ":", str(awtune.C), sep='', file=sys.stdout)
#                 OUTCYCLE = OUTCYCLE + 1
#                 return awtune.runGaussian()
#
#
#             res = optimize.fminbound(awtune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)
#
#         # Scan of Alpha or Omega
#         if args.tuning == 's':
#             if args.alpha:
#                 for alpha in args.alpha.split():
#                     tuningObj.alpha = float(alpha)
#                     if args.omega:
#                         for omega in args.omega.split():
#                             tuningObj.omega = float(omega)
#                             tuningObj.runGaussian()
#                     else:
#                         tuningObj.runGaussian()
#             elif args.omega:
#                 for omega in args.omega.split():
#                     tuningObj.omega = float(omega)
#                     tuningObj.runGaussian()
#             else:
#                 print("You must set either alpha (-a) or omega (-w) to run a scan job.", file=sys.stderr)
#                 print("Exiting...", file=sys.stderr)
#                 sys.exit(1)
#
#         shutil.copy(tuningObj.workingDIR + '/' + tuningObj.neutralFilename + '.gjf',
#                     tuningObj.currentDIR + '/' + tuningObj.originFilename + '_' + re.sub('[+)(,]', '',
#                                                                                          tuningObj.level).replace('/',
#                                                                                                                   '_') + '_' + tuningObj.criterion + '_tuned.gjf')
