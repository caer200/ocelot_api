import sys
import os
import re
import subprocess
import argparse
from scipy import optimize
import shutil
import time
from pymatgen.io.gaussian import GaussianInput
"""
originally written by Sean M. Ryno, Cheng Zhong, Haitao Sun
for a certain mol, get tuned a or w
"""

class GaussTuning:

    def __init__(self, omol, gparams):

        self.omol = omol
        self.inp = GaussianInput(self.omol.to_pymatgen_mol(), charge=omol.charge, spin_multiplicity=omol.multiplicity)











class tuningGaussian():
    def __init__(self, filename, nproc='', mem='', level='', addKeyWords='', criterion='', program=''):
        suffix = filename.split('.').pop().lower()
        if (suffix == 'gjf') or (suffix == 'com'):
            self.filename = filename
        else:
            print("Filename should end with .gjf or .com", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)

        print("Filename:", filename, file=sys.stdout)

        try:
            self.inputFile = open(filename, 'r')
        except IOError:
            print("Cannot open file:", filename, file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)

        # Read nproc, mem, methods, and coordinates from input file.
        self.coord = []
        readCoord = 0
        keywordLine = 0

        # Read Input file first
        for line in self.inputFile:
            if re.search('nproc', line):
                nprocInFile = line.split('=')[1].strip('\n')
            if re.search('nprocshared', line):
                nprocInFile = line.split('=')[1].strip('\n')
            if re.search('CPU', line):
                nprocInFile = line.split('=')[1].strip('\n')
            if re.search('mem', line):
                memInFile = line.split('=')[1].strip('\n')
            if re.search('#', line):
                keywordLine = keywordLine + 1
                if keywordLine > 1:
                    print("Only a single keyword line is allowed and must begin with functional/basisset.",
                          file=sys.stderr)
                    print("Exiting...", file=sys.stderr)
                    sys.exit(1)
                try:
                    levelInFile = line.split(None, 2)[1].strip('\n')
                except IndexError:
                    levelInFile = ''
                try:
                    addKeyWordsInFile = line.split(None, 2)[2].strip('\n')
                except IndexError:
                    addKeyWordsInFile = ''
            if re.search(r'^-?[0-9] +[0-9] *$', line):
                self.charge = line.split()[0]
                self.spin = line.split()[1]
                readCoord = 1
                continue
            if readCoord == 1:
                self.coord.append(line)
            else:
                readCoord = 0

        # Set number of processors to use. Default: 1
        if nproc != '':
            self.nrpoc = nproc
        elif nprocInFile != '':
            self.nproc = nprocInFile
        else:
            self.nproc = '1'

        # Set amount of memory to use. Default: 500MB
        if mem != '':
            self.mem = mem
        elif memInFile != '':
            self.mem = memInFile
        else:
            self.mem = '500MB'

        # Set the level of the calculations. Default: LC-wPBE/6-31G(d)
        if level != '':
            self.level = level
        elif levelInFile != '':
            self.level = levelInFile
        else:
            self.level = 'LC-wPBE/6-31G(d)'

        # Include any additional keywords that have been specified
        if addKeyWords != '':
            self.addKeyWords = addKeyWords
        elif addKeyWordsInFile.strip() != '':
            self.addKeyWords = addKeyWordsInFile
        else:
            self.addKeyWords = 'scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)'

        # Set the program to be used
        if program != '':
            self.program = program
        else:
            self.program = 'g09'

        # Set other variables needed for the tuning procedure
        self.currentDIR = os.getcwd()
        self.criterion = criterion
        self.workingDIR = self.filename.rpartition('.')[0]
        self.originFilename = ''.join(re.split('[/.]', self.filename)[-2:-1])
        self.startTime = time.time()
        self.cycle = 0
        self.conver = 4
        self.alphaConver = 4
        self.omega = '-----'
        self.alpha = '-----'
        self.infoTable = ''.join(re.split('[/.]', self.filename)[-2:-1]) + '_' + re.sub('[+)(,]', '',
                                                                                        self.level).replace('/',
                                                                                                            '_') + '_' + self.criterion + '_awDATA'

    # Generate Gaussian input file, run Gaussian, and extract the necessary data.
    def genGJF(self, tuningInput, tuningNameSuffix, noGuess=0):

        # Setup initial spin and charges
        neutralCharge = self.charge
        neutralSpin = self.spin
        cationCharge = str(int(self.charge) + 1)
        anionCharge = str(int(self.charge) - 1)
        if int(neutralSpin) == 1:
            cationSpin = '2'
            anionSpin = '2'
        if int(neutralSpin) > 1:
            cationSpin = str(int(neutralSpin) - 1)
            anionSpin = str(int(neutralSpin) - 1)

        # Determine if this is the first run or subsequent run.
        if (self.cycle >= 1) and (noGuess != 1):
            ifguess = 'guess=read'
        else:
            ifguess = ''

        # Write neutral input file
        methods = self.level + ' ' + self.addKeyWords + ' ' + ifguess
        neutralInput = open(self.neutralFilename + '.gjf', 'w')
        neutralInput.write('%chk=' + self.neutralFilename + '.chk\n')
        neutralInput.write('%nprocshared=' + self.nproc + '\n')
        neutralInput.write('%mem=' + self.mem + '\n')
        neutralInput.write('#p ' + methods + ' ' + tuningInput + '\n')
        neutralInput.write('\n')
        neutralInput.write('Tuning input file generated by awtune with parameters: ' + tuningNameSuffix + '\n')
        neutralInput.write('\n')
        neutralInput.write(neutralCharge + ' ' + neutralSpin + '\n')
        for coords in self.coord:
            neutralInput.write(coords)
        neutralInput.write('\n')

        # Write cation input file
        cationInput = open(self.cationFilename + '.gjf', 'w')
        cationInput.write('%chk=' + self.cationFilename + '.chk\n')
        cationInput.write('%nprocshared=' + self.nproc + '\n')
        cationInput.write('%mem=' + self.mem + '\n')
        cationInput.write('#p ' + methods + ' ' + tuningInput + '\n')
        cationInput.write('\n')
        cationInput.write('Tuning input file generated by awtune with parameters: ' + tuningNameSuffix + '\n')
        cationInput.write('\n')
        cationInput.write(cationCharge + ' ' + cationSpin + '\n')
        for coords in self.coord:
            cationInput.write(coords)
        cationInput.write('\n')

        # Write anion input file
        anionInput = open(self.anionFilename + '.gjf', 'w')
        anionInput.write('%chk=' + self.anionFilename + '.chk\n')
        anionInput.write('%nprocshared=' + self.nproc + '\n')
        anionInput.write('%mem=' + self.mem + '\n')
        anionInput.write('#p ' + methods + ' ' + tuningInput + '\n')
        anionInput.write('\n')
        anionInput.write('Tuning input file generated by awtune with parameters: ' + tuningNameSuffix + '\n')
        anionInput.write('\n')
        anionInput.write(anionCharge + ' ' + anionSpin + '\n')
        for coords in self.coord:
            anionInput.write(coords)
        anionInput.write('\n')

        # Write closing line breaks and close files
        neutralInput.write('\n')
        cationInput.write('\n')
        anionInput.write('\n')
        neutralInput.close()
        cationInput.close()
        anionInput.close()

    def extractEHL(self, filename):
        # Extract required information from log files
        try:
            LogFile = open(filename + '.log', 'r')
        except IOError:
            # print("There was a problem reading: ", filename + '.log', file=sys.stderr)
            return ''

        normal = 0
        for line in LogFile:
            if re.search('Normal termination', line):
                normal = 1
            if re.search('SCF Done', line):
                Energy = line.split()[4]
            if re.search('Alpha  occ. eigenvalue', line):
                AlphaHOMO = line.split()[-1]
                occFlag = 1
            if re.search('Beta  occ. eigenvalues', line):
                BetaHOMO = line.split()[-1]
                occFlag = 1
            if re.search('Alpha virt. eigenvalues', line) and (occFlag == 1):
                AlphaLUMO = line.split()[4]
                occFlag = 0
            if re.search('Beta virt. eigenvalues', line) and (occFlag == 1):
                BetaLUMO = line.split()[4]
                occFlag = 0

        if normal == 0:
            print("Job:", filename, " Exited with Abnormal Termination.", file=sys.stderr)
            return ''
        try:
            BetaHOMO
        except NameError:
            BetaHOMO = AlphaHOMO
        try:
            BetaLUMO
        except NameError:
            BetaLUMO = AlphaLUMO

        if float(BetaHOMO) > float(AlphaHOMO):
            HOMO = float(BetaHOMO)
        else:
            HOMO = float(AlphaHOMO)
        if float(BetaLUMO) < float(AlphaLUMO):
            LUMO = float(BetaLUMO)
        else:
            LUMO = float(AlphaLUMO)

        return [float(Energy) * 27.211396132, HOMO * 27.211396132, LUMO * 27.211396132]

    def runGaussian(self):
        print(' ', file=sys.stdout)

        # Write Alpha and Omega values and perform a sanity check for the values.
        tuningInput = ''
        tuningNameSuffix = ''
        if self.omega != '-----':
            if (float(self.omega) < 10) and (float(self.omega) >= 0):
                self.omega = round(self.omega, self.conver)
                tuningInput = tuningInput + 'IOP(3/107=' + str(int(float(self.omega) * 10000)).zfill(
                    5) + '00000) IOP(3/108=' + str(int(float(self.omega) * 10000)).zfill(5) + '00000) '
            else:
                print("Omega", self.omega, "must be a number between 0 and 10", file=sys.stderr)
                print("Exiting...", file=sys.stderr)
                sys.exit(1)
            tuningNameSuffix = tuningNameSuffix + 'w' + str(int(float(self.omega) * 10000)).zfill(5)

        if self.alpha != '-----':
            if (float(self.alpha) <= 1) and (float(self.alpha) >= 0):
                self.alpha = round(self.alpha, self.alphaConver)
                beta = 1 - float(self.alpha)
                tuningInput = tuningInput + 'IOP(3/130=' + str(int(self.alpha * 10000)).zfill(5) + ') IOP(3/131=' + str(
                    int(self.alpha * 10000)).zfill(5) + ') '
                tuningInput = tuningInput + 'IOP(3/119=' + str(int(beta * 10000)).zfill(5) + '00000) IOP(3/120=' + str(
                    int(beta * 10000)).zfill(5) + '00000) '
            else:
                print("Alpha", self.alpha, "must be a number between 0 and 1", file=sys.stderr)
                print("Exiting...", file=sys.stderr)
                sys.exit(1)
            tuningNameSuffix = 'a' + str(int(float(self.alpha) * 10000)).zfill(5) + tuningNameSuffix

        # Move old checkpoint files and set new names
        try:
            oldNF, oldCF, oldAF = self.neutralFilename, self.cationFilename, self.anionFilename
        except:
            oldNF, oldCF, oldAF = '', '', ''

        self.neutralFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
                                                                                                    '_') + '_' + tuningNameSuffix + '_n'
        self.cationFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
                                                                                                   '_') + '_' + tuningNameSuffix + '_c'
        self.anionFilename = self.originFilename + '_' + re.sub('[+)(,]', '', self.level).replace('/',
                                                                                                  '_') + '_' + tuningNameSuffix + '_a'
        self.tuningNameSuffix = tuningNameSuffix

        # Search the InfoTable, if data is already there extract it from the table
        # If alpha and omega are the same as the previous cycle do not increase cycle number
        try:
            infoTable = open(self.infoTable, 'r')
            cycle = self.cycle
            lastLine = 1
            for line in infoTable:
                if self.cycle == cycle - 1:
                    lastLine = 0
                    break
                if (str(self.alpha) == line.split()[0]) and (str(self.omega) == line.split()[1]):
                    print("alpha:", str(self.alpha), ' omega:', str(self.omega), ' found in ', self.infoTable, sep='',
                          file=sys.stdout)
                    self.C, self.J2, self.Jn2, self.O2, self.gap, self.IP, self.EA, self.neutralE, self.N_HOMO, self.N_LUMO, self.cationE, self.C_HOMO, self.C_LUMO, self.anionE, self.A_HOMO, self.A_LUMO = [
                        float(i) for i in line.split()[2:]]
                    cycle = cycle + 1
            if lastLine == 0:
                self.cycle = cycle
                print("-----cycle", str(self.cycle).zfill(2), "-----", tuningNameSuffix, "  ", self.criterion, ":",
                      str(self.C), "  gap:", str(self.gap), sep='', file=sys.stdout)
                return self.C
            if self.cycle == cycle - 1:
                return self.C
            infoTable.close()
        except IOError:
            pass

        # Create a working directory to perform the tuning in
        if not os.path.exists(self.workingDIR):
            os.makedirs(self.workingDIR)
        os.chdir(self.workingDIR)

        # Create Gaussian input file
        self.genGJF(tuningInput, tuningNameSuffix)

        # Call Gaussian
        if len(self.extractEHL(self.neutralFilename)) != 3:
            print("Running: ", self.neutralFilename, ".gjf...",
                  time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
            if os.path.isfile(oldNF + '.chk'):
                shutil.copy(oldNF + '.chk', self.neutralFilename + '.chk')
            else:
                self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
            runNeutral = subprocess.Popen([self.program, self.neutralFilename + '.gjf'])
            runNeutral.wait()
        else:
            print("logfile ", self.neutralFilename, ".log found, directly extracting information.", sep='',
                  file=sys.stdout)
        try:
            neutralE, N_HOMO, N_LUMO = self.extractEHL(self.neutralFilename)
        except ValueError:
            print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.neutralFilename,
                  ".log. Please check this files.", sep='', file=sys.stderr)
            print(
                "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
                file=sys.stderr)
            sys.exit(1)
        if self.criterion != "Jl":
            if len(self.extractEHL(self.cationFilename)) != 3:
                print("Running ", self.cationFilename, ".gjf...",
                      time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
                if os.path.isfile(oldCF + ".chk"):
                    shutil.copy(oldCF + ".chk", self.cationFilename + ".chk")
                else:
                    self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
                runCation = subprocess.Popen([self.program, self.cationFilename + ".gjf"])
                runCation.wait()
            else:
                print("logfile ", self.cationFilename, ".log found, directly extracting information.", sep='',
                      file=sys.stdout)
            try:
                cationE, C_HOMO, C_LUMO = self.extractEHL(self.cationFilename)
            except ValueError:
                print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.cationFilename,
                      ".log. Please check this files.", sep='', file=sys.stderr)
                print(
                    "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
                    file=sys.stderr)
                sys.exit(1)
        else:
            cationE, C_HOMO, C_LUMO = 0, 0, 0
        if self.criterion != "Jh":
            if len(self.extractEHL(self.anionFilename)) != 3:
                print("Running ", self.anionFilename, ".gjf...",
                      time.strftime('%H:%M:%S', time.gmtime(time.time() - self.startTime)), sep='', file=sys.stdout)
                if os.path.isfile(oldAF + ".chk"):
                    shutil.copy(oldAF + ".chk", self.anionFilename + ".chk")
                else:
                    self.genGJF(tuningInput, tuningNameSuffix, noGuess=1)
                runAnion = subprocess.Popen([self.program, self.anionFilename + ".gjf"])
                runAnion.wait()
            else:
                print("logfile ", self.anionFilename, ".log found, directly extracting information.", sep='',
                      file=sys.stdout)
            try:
                anionE, A_HOMO, A_LUMO = self.extractEHL(self.anionFilename)
            except ValueError:
                print("ERROR: Cannot extract information from ", self.workingDIR, "/", self.anionFilename,
                      ".log. Please check this files.", sep='', file=sys.stderr)
                print(
                    "It is likely that the job either failed or that your initial Gaussian input file is DOS formatted instead of UNIX formatted.",
                    file=sys.stderr)
                sys.exit(1)
        else:
            anionE, A_HOMO, A_LUMO = 0, 0, 0

        os.chdir(self.currentDIR)
        if (self.criterion != "Jh") and (self.criterion != "Jl"):
            IP = neutralE - cationE
            EA = neutralE - anionE
            J2 = (N_HOMO - IP) ** 2 + (A_HOMO + EA) ** 2
            J = J2 ** 0.5
            gap = N_LUMO - N_HOMO
            Jh = abs(N_HOMO - IP)
            Jl = abs(N_LUMO + EA)
            Jn2 = Jh ** 2 + Jl ** 2
            Jn = Jn2 ** 0.5
            O2 = (N_HOMO - C_LUMO) ** 2 + (N_LUMO - A_HOMO) ** 2
            C = eval(self.criterion)
        else:
            IP = 0
            EA = 0
            J2 = 0
            J = J2 ** 0.5
            Jl = 0
            Jn2 = 0
            Jn = Jn2 ** 0.5
            O2 = 0
        if self.criterion == "Jh":
            IP = neutralE - cationE
            gap = N_LUMO - N_HOMO
            Jh = abs(N_HOMO - IP)
            C = eval(self.criterion)
        if self.criterion == "Jl":
            EA = neutralE - anionE
            gap = N_LUMO - N_HOMO
            Jl = abs(N_LUMO + EA)
            C = eval(self.criterion)

        # Assign self values based on new run
        self.C, self.J, self.J2, self.Jn, self.Jn2, self.O2, self.gap, self.IP, self.EA, self.neutralE, self.N_HOMO, self.N_LUMO, self.cationE, self.C_HOMO, self.C_LUMO, self.anionE, self.A_HOMO, self.A_LUMO = C, J, J2, Jn, Jn2, O2, gap, IP, EA, neutralE, N_HOMO, N_LUMO, cationE, C_HOMO, C_LUMO, anionE, A_HOMO, A_LUMO
        lenNE = str(len('{:.6f}'.format(self.neutralE)) + 2)
        lenCE = str(len('{:.6f}'.format(self.cationE)) + 2)
        lenAE = str(len('{:.6f}'.format(self.anionE)) + 2)

        # Make InfoTable
        try:
            infoTable = open(self.infoTable, 'r')
        except IOError:
            infoTable = open(self.infoTable, 'w')
            infoTable.write((
                '{:<10}{:<10}{:<11}{:<11}{:<11}{:<11}{:<11}{:<11}{:<11}{:<' + lenNE + '}{:<11}{:<11}{:<' + lenCE + '}{:<11}{:<11}{:<' + lenAE + '}{:<11}{:<11}\n').format(
                'Alpha', 'Omega', 'opt:' + self.criterion, 'J2', 'Jn2', 'O2', 'gap', 'IP', 'EA', 'neutralE', 'N_HOMO',
                'N_LUMO', 'cationE', 'C_HOMO', 'C_LUMO', 'anionE', 'A_HOMO', 'A_LUMO'))
        infoTable.close()
        infoTable = open(self.infoTable, 'a')
        infoTable.write((
            '{:<10}{:<10}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<11.6f}{:<' + lenNE + '.6f}{:<11.6f}{:<11.6f}{:<' + lenCE + '.6f}{:<11.6f}{:<11.6f}{:<' + lenAE + '.6f}{:<11.6f}{:<11.6f}\n').format(
            self.alpha, self.omega, C, J2, Jn2, O2, gap, IP, EA, neutralE, N_HOMO, N_LUMO, cationE, C_HOMO, C_LUMO,
            anionE, A_HOMO, A_LUMO))
        infoTable.close()
        infoTable = open(self.infoTable, 'r')
        self.cycle = self.cycle + 1

        print("-----cycle", str(self.cycle).zfill(2), "-----", tuningNameSuffix, "  ", self.criterion, ":", str(C),
              "  gap:", str(gap), sep='', file=sys.stdout)

        return C

    def __str__(self):
        return self.filename + " charge:" + self.charge + " spin:" + self.spin + " nproc=" + self.nproc + " mem=" + self.mem + " level:" + self.level


# Setup Parser and help files
parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                 description="A script used with Gaussian09/16 to auto optimize omega and alpha for long-range corrected DFT methodologies. This program would not be possible without the original version by Cheng Zhong.")
parser.add_argument("-g", "--gauss", help="Set the executable name to be used. Default: g09", type=str, default='g09')
parser.add_argument("-l", "--level",
                    help="Set the theoretical level. Default: read from gjf file. If not found, use LC-wPBE/6-31g(d)",
                    type=str, default='')
parser.add_argument("-m", "--memory",
                    help="Set the memory used by Gaussian (with unit).  Default: read from gjf file. If not found, use: 500MB",
                    type=str, default='')
parser.add_argument("-n", "--nproc",
                    help="Set the number of cores used by Gaussian. Default: read from gjf file. If not found, use 1",
                    type=str, default='')
parser.add_argument("-p", "--precision",
                    help="Set the precision of omega. Default: 4, means omega will converge to 0.0001 ", type=int,
                    default='4')
parser.add_argument("-P", "--alphaPrecision", help="Set the precision of alpha. Default: equal to precision of omega",
                    type=int)
parser.add_argument("-b", "--bound", help="Set the boundary for omega, e.g. -b 0.05,0.3. Default: 0.05,0.5", type=str,
                    default='0.05,0.5')
parser.add_argument("-B", "--alphaBound", help="Set the boundary for alpha, e.g. -B 0.05,0.3. Default: 0.0,0.5",
                    type=str, default='0.0,0.5')
parser.add_argument("-k", "--keywords",
                    help="Set additional Gaussian keywords. Default: \"scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)\"",
                    type=str, default='')
parser.add_argument("-c", "--criterion",
                    help="Set the optimization criterion (the value to minimize), available options are:\nJ2---((HOMO-IP)^2+(A_HOMO+EA)^2)\nJh---(HOMO-IP)\nJl(LUMO+EA)\nJn2---((HOMO-IP)^2+(LUMO+EA)^2)\nO2---((A_HOMO-LUMO)^2+(C_LUMO-HOMO)^2)\nor any other valid experssion, e.g. -c \"(A_HOMO-N_LUMO)**2\"\nIf multiple criterion is used, seperate them by comma. The order and number of criterions must be the same to that of structures. Default: J2",
                    type=str, default='J2')
parser.add_argument("-t", "--tuning", help="Set the tuning scheme. Available options are: w(omega),a(alpha),aw(alpha and omega),s(scan alpha or omege, must be used with -a or -w option)\nPut m in front of above options (mw,ma,maw,ms) means tuning with multiple structure,\n\
e.g. \"-t mw\" means optimize the omega to minimize sum of the criterion of all the input structure.\nIf different criterion is used for different structure, use -c option to designate criterion explicitly for each structure, \ne.g. \"-t ma -c Jh,Jl\" means optimize alpha to minimize sum of Jh of structure 1 and Jl of structure 2.  Default: w",
                    type=str, default='w')
parser.add_argument("-a", "--alpha",
                    help="Set the value of alpha. Can be a single value or multiple value seperated by space, e.g. -a '0.1 0.3 0.5 0.7' Default: None",
                    type=str)
parser.add_argument("-w", "--omega",
                    help="Set the value of omega. Can be a single value or multiple value seperated by space. Default: None",
                    type=str)
parser.add_argument("gaussianInput", nargs='*')
args = parser.parse_args()

# Set values based on parser input
program = args.gauss
level = args.level
memory = args.memory
nproc = args.nproc
bound = [float(i) for i in args.bound.split(',')]
alphaBound = [float(i) for i in args.alphaBound.split(',')]
keywords = args.keywords
tuning = args.tuning
criterions = args.criterion.split(',')

# Do check to see if multiple files are being requested to run and run them
# Tuning for single structures
if 'm' not in args.tuning:
    if len(criterions) > 1:
        print("ERROR: Multiple criterion have been requested, but this must be used with the '-t m; option.",
              file=sys.stderr)
        print("Exiting...", file=sys.stderr)
        sys.exit(1)
    for g09input in args.gaussianInput:
        tuningObj = tuningGaussian(g09input, nproc=nproc, mem=memory, level=level, addKeyWords=keywords,
                                   criterion=args.criterion, program=program)
        tuningObj.conver = args.precision
        if not args.alphaPrecision:
            args.alphaPrecision = args.precision
        tuningObj.alphaConver = args.alphaPrecision

        # Now run the requested tuning procedure and optimize as requested
        # Omega tuning
        if args.tuning == 'w':
            def wtune(omega, wtune=tuningObj):
                wtune.omega = omega
                return wtune.runGaussian()


            if args.alpha:
                for alpha in args.alpha.split():
                    tuningObj.alpha = float(alpha)
                    res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)
            else:
                res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)

        # Alpha tuning
        if args.tuning == 'a':
            def atune(alpha, atune=tuningObj):
                atune.alpha = alpha
                return atune.runGaussian()


            if args.omega:
                for omega in args.omega.split():
                    tuningObj.omega = float(omega)
                    res = optimize.fminbound(atune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)
            else:
                res = optimize.fminbound(atune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)

        # Alpha and Omega tuning
        if args.tuning == 'aw':
            OUTCYCLE = 1


            def awtune(alpha, awtune=tuningObj):
                global OUTCYCLE

                def wtune(omega, wtune=tuningObj):
                    global OUTCYCLE
                    wtune.omega = omega
                    return wtune.runGaussian()

                awtune.alpha = alpha
                res = optimize.fminbound(wtune, bound[0], bound[1], xtol=1e-05, disp=3)
                awtune.omega = res
                print("=====OUTCYCLE", str(OUTCYCLE).zfill(2), "=====", awtune.tuningNameSuffix, "  ", awtune.criterion,
                      ":", str(awtune.C), sep='', file=sys.stdout)
                OUTCYCLE = OUTCYCLE + 1
                return awtune.runGaussian()


            res = optimize.fminbound(awtune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)

        # Scan of Alpha or Omega
        if args.tuning == 's':
            if args.alpha:
                for alpha in args.alpha.split():
                    tuningObj.alpha = float(alpha)
                    if args.omega:
                        for omega in args.omega.split():
                            tuningObj.omega = float(omega)
                            tuningObj.runGaussian()
                    else:
                        tuningObj.runGaussian()
            elif args.omega:
                for omega in args.omega.split():
                    tuningObj.omega = float(omega)
                    tuningObj.runGaussian()
            else:
                print("You must set either alpha (-a) or omega (-w) to run a scan job.", file=sys.stderr)
                print("Exiting...", file=sys.stderr)
                sys.exit(1)

        shutil.copy(tuningObj.workingDIR + '/' + tuningObj.neutralFilename + '.gjf',
                    tuningObj.currentDIR + '/' + tuningObj.originFilename + '_' + re.sub('[+)(,]', '',
                                                                                         tuningObj.level).replace('/',
                                                                                                                  '_') + '_' + tuningObj.criterion + '_tuned.gjf')

# Tuning for multiple structures
if 'm' in args.tuning:
    allStruct = []
    if len(criterions) > 1:
        if len(args.gaussianInput) != len(criterions):
            print("ERROR: Number of structures should be equal to the number of criterion.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
        combinedCriterion = '+'.join(criterions)
    elif len(criterions) == 1:
        combinedCriterion = 'sum of ' + args.criterion
        criterions = criterions * len(args.gaussianInput)

    for i in range(len(args.gaussianInput)):
        tuningObj = tuningGaussian(args.gaussianInput[i], nproc=nproc, mem=memory, level=level, addKeyWords=keywords,
                                   criterion=criterions[i], program=program)
        tuningObj.conver = args.precision
        if not args.alphaPrecision:
            args.alphaPrecision = args.precision
        tuningObj.alphaConver = args.alphaPrecision
        allStruct.append(tuningObj)

    CinfoTable = ''
    listFilename = []
    for i in allStruct:
        listFilename.append(''.join(re.split('[/.]', i.filename)[-2:-1]) + '_' + i.criterion)
        CinfoTable = CinfoTable + ''.join(re.split('[/.]', i.filename)[-2:-1]) + '_'
    CinfoTable = CinfoTable + ''.join(args.criterion.split(',')) + '_cawDATA'
    infoTable = open(CinfoTable, 'w')
    lenC = len(combinedCriterion) + 7
    lenList_raw = [len(i) for i in listFilename]
    lenList = []
    for i in lenList_raw:
        if i < 8:
            i = 8
        lenList.append(str(i + 2))
    fnFormat = '}{:<'.join(lenList)

    infoTable.write(
        ('{:<10}{:<10}{:<' + str(lenC) + '}{:<' + fnFormat + '}\n').format('Alpha', 'Omega', 'opt:' + combinedCriterion,
                                                                           *listFilename))

    # Multiple structure Omega tuning
    if args.tuning == 'mw':
        def cwtune(omega, cwtune=allStruct):
            sumC = 0
            listDATA = []
            for i in cwtune:
                i.omega = omega
                result = i.runGaussian()
                listDATA.append(str(round(result, 6)))
                sumC = sumC + result
            print("=====COMBINED-CYCLE", str(i.cycle).zfill(2), "=====", i.tuningNameSuffix, "  ", combinedCriterion,
                  ":", str(sumC), sep='', file=sys.stdout)
            print('', file=sys.stdout)
            infoTable.write(
                ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega, sumC,
                                                                                      *listDATA))
            return sumC


        if args.alpha:
            for alpha in args.alpha.split():
                for i in allStruct:
                    i.alpha = float(alpha)
                res = optimize.fminbound(cwtune, bound[0], bound[1], xtol=1e-05, disp=3)
        else:
            res = optimize.fminbound(cwtune, bound[0], bound[1], xtol=1e-05, disp=3)

    # Multiple structure Alpha tuning
    if args.tuning == 'ma':
        def catune(alpha, catune=allStruct):
            sumC = 0
            listDATA = []
            for i in catune:
                i.alpha = alpha
                result = i.runGaussian()
                listDATA.append(str(round(result, 6)))
                sumC = sumC + result
            print("=====COMBINED-CYCLE", str(i.cycle).zfill(2), "=====", i.tuningNameSuffix, "  ", combinedCriterion,
                  ":", str(sumC), sep='', file=sys.stdout)
            print('', file=sys.stdout)
            infoTable.write(
                ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega, sumC,
                                                                                      *listDATA))
            return sumC


        if args.omega:
            for omega in args.omega.split():
                for i in allStruct:
                    i.omega = float(omega)
                res = optimize.fminbound(catune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)
        else:
            res = optimize.fminbound(catune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)

    # Multiple structure Alpha and Omega tuning
    if args.tuning == 'maw':
        OUTCYCLE = 1


        def cawtune(alpha, cawtune=allStruct):
            global OUTCYCLE

            def cwtune(omega, cwtune=allStruct):
                global OUTCYCLE
                sumC = 0
                listDATA = []
                for i in cwtune:
                    i.omega = omega
                    result = i.runGaussian()
                    listDATA.append(str(round(result, 6)))
                    sumC = sumC + result
                    print("=====COMBINED-CYCLE", str(OUTCYCLE), '_', str(i.cycle).zfill(2), "=====", i.tuningNameSuffix,
                          "  ", combinedCriterion, ":", str(sumC), sep='', file=sys.stdout)
                    print('', file=sys.stdout)
                    infoTable.write(
                        ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega, sumC,
                                                                                              *listDATA))
                    return sumC
                sumC = 0
                for i in cawtune:
                    i.alpha = alpha
                res = optimize.fminbound(cwtune, bound[0], bound[1], xtol=1e-05, disp=3)
                for i in cawtune:
                    i.omega = res
                    sumC = sumC + i.runGaussian()
                print("=====COMBINED-OUTERCYCLE", str(OUTCYCLE).zfill(2), "=====", i.tuningNameSuffix, '  ',
                      combinedCriterion, ':', str(sumC), sep='', file=sys.stdout)
                print('', file=sys.stdout)
                OUTCYCLE = OUTCYCLE + 1
                return sumC

            res = optimize.fminbound(cawtune, alphaBound[0], alphaBound[1], xtol=1e-05, disp=3)

    # Multiple structure scan
    if args.tuning == 'ms':
        if args.alpha:
            for alpha in args.alpha.split():
                sumC = 0
                listDATA = []
                for i in allStruct:
                    i.alpha = float(alpha)
                    if args.omega:
                        for omega in args.omega.split():
                            for i in allStruct:
                                i.omega = float(omega)
                                result = i.runGaussian()
                                listDATA.append(str(round(result, 6)))
                                sumC = sumC + result
                            infoTable.write(
                                ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega,
                                                                                                      sumC, *listDATA))
                    else:
                        result = i.runGaussian()
                        listDATA.append(str(round(result, 6)))
                        sumC = sumC + result
                print(listDATA, file=sys.stdout)
                infoTable.write(
                    ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega, sumC,
                                                                                          *listDATA))
        elif args.omega:
            for omega in args.omega.split():
                sumC = 0
                listDATA = []
                for i in allStruct:
                    i.omega = float(omega)
                    result = i.runGaussian()
                    listDATA.append(str(round(result, 6)))
                    sumC = sumC + result
                infoTable.write(
                    ('{:<10}{:<10}{:<' + str(lenC) + '.6f}{:<' + fnFormat + '}\n').format(i.alpha, i.omega, sumC,
                                                                                          *listDATA))
        else:
            print("Alpha (-a) or Omega (-o) must be set to initiate scan job.", file=sys.stderr)
            print("Exiting...", file=sys.stderr)
            sys.exit(1)
