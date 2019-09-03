import os
import glob
import warnings
import sys
import configuration, shellcall, vasp, gauss
import rawcifparser
import settings
import fileop
import datetime


class Rocket:

    def __init__(self, name, tankload):

        self.name = name
        self.cifname = self.name + '.cif'

        self.globalpath = settings.SYSPATH
        self.tankload = tankload

        # the path where calculations are done, eg. ./run/23/
        self.binpath = self.globalpath['BIN_PATH']
        self.runpath = self.globalpath['RUN_PATH'] + '/' + self.name + '/'
        fileop.createdir(self.runpath)

        f = open(self.globalpath['CIF_REPO'] + '/' + self.cifname, 'r')
        self.cifstring = f.read()
        f.close()

        # check path in settings
        if not fileop.check_path([self.globalpath[k] for k in self.globalpath.keys()]):
            warnings.warn('check system path!')
            sys.exit(1)

        self.logfile = self.runpath + '/' + self.name + '.rocketlog'
        with open(self.logfile, 'w') as f:
            f.write('# log file initialized at ' + datetime.datetime.now().strftime('%H:%M:%S on %Y-%m-%d') + '\r\n')

        self.log_append('# tanks are loaded as: ' + ' & '.join([str(i) for i in self.tankload]))

    @property
    def rocket_status(self):
        """
        :return: boolean for whether all tankload empty
        """
        status = 0  # clean launchpad
        emptyload = []
        for i in self.tankload:
            itank = Tank(self, i)
            if itank.status:
                emptyload.append(i)
        if tuple(emptyload) == tuple(self.tankload):
            status = 1
        return status

    def cleanup(self):
        os.system('rm -rf ' + self.runpath)

    def launch(self):
        fileop.createdir(self.runpath)
        if not self.rocket_status:
            print('rocket launch with tankload ', self.tankload)
            for i in self.tankload:
                current_tank = Tank(self, i)
                logstring = current_tank.ignite()
                if not current_tank.status:
                    warnings.warn('tank ' + str(i) + 'exploded at ' + os.getcwd())
                    sys.exit(1)
                self.log_append(self.name + logstring)

    def log_append(self, string):
        with open(self.logfile, 'a') as f:
            f.write(string)


class Tank:

    def __init__(self, rocket, stepnumber):
        self.runpath = rocket.runpath
        self.cifstring = rocket.cifstring
        self.stepnumber = stepnumber
        self.globalpath = rocket.globalpath

        self.wdir = self.runpath + '/' + 's' + str(self.stepnumber) + '/'

    def ignite(self):
        exhaust = ' \r\n step ' + str(self.stepnumber) + '\r\n'

        whereiam = os.getcwd()
        fileop.createdir(self.wdir)
        os.chdir(self.wdir)
        if self.status:
            exhaust = exhaust + 'already finished \r\n'
        else:
            dustbin = 'dustbin' + datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + '/'
            fileop.createdir(self.wdir + '/' + dustbin)
            os.system('mv * ./' + dustbin)
            self.workstep(self.stepnumber, self.cifstring, self.globalpath['BIN_PATH']+'/aflow',
                          self.wdir, self.globalpath['PAW_PATH'])
            exhaust = exhaust + 'finished at' + str(datetime.datetime.now()) + '\r\n'
        os.chdir(whereiam)

        return exhaust

    @property
    def status(self):
        fin = 0

        if not os.path.isdir(self.wdir):
            return fin

        else:
            whereiam = os.getcwd()
            os.chdir(self.wdir)
            if self.stepnumber == 0:
                files = ['r.cif', 'config0.cif', 'kpath_unwrap.poscar', 'mol0.poscar', 'mol0.xyz', 'molbone0.xyz',
                         'unwrap.poscar', 'std_unwrap.poscar', 'dimers.pkl']
                if all([os.path.isfile(f) for f in files]):
                    if all([os.stat(f).st_size > 0 for f in files]):
                        fin = 1
            elif self.stepnumber in (1, 2, 3, 4, 5, 6, 7):  # 1~7 are vasp, check outcar
                if fileop.check_outcar() and os.path.isfile('vasprun.xml'):
                    fin = 1
            elif self.stepnumber in (8, 9, 10, 11, 12):
                if fileop.check_allgausslog(len(glob.glob('*.com'))):
                    fin = 1
            os.chdir(whereiam)
            return fin

    @staticmethod
    def workstep(stepnumber, cifstring, aflowpath, wdir, pawpath):
        whereami = os.getcwd()
        os.chdir(wdir)

        if stepnumber == 0:  # config ana & aflow std
            with open('r.cif', 'w') as f:
                f.write(cifstring)
            rawcifparser.RawCifParser('r.cif')
            ca = configuration.Configuration('config0.cif')
            ca.unwrap_structure.to(fmt='poscar', filename='unwrap.poscar')
            for i in range(len(ca.boxmols)):
                ca.boxmols[i].to(fmt='poscar', filename='mol' + str(i) + '.poscar')
            for i in range(len(ca.mols)):
                ca.mols[i].to(fmt='xyz', filename='mol' + str(i) + '.xyz')
                ca.organic_mols[i].bone.mol.to(fmt='xyz', filename='molbone' + str(i) + '.xyz')
            fileop.write_obj(ca.dimers, 'dimers.pkl')
            shellcall.aflow(aflowpath, 'std', poscar='unwrap.poscar')

        elif stepnumber == 1:
            fileop.copyfile('../s0/mol0.poscar', './mol0.poscar')
            fileop.copyfile('./mol0.poscar', './POSCAR')
            vasp.WriteInp('fbrlx', pawpath)
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 2:
            fileop.copyfile('../s0/std_unwrap.poscar', './')
            fileop.copyfile('./std_unwrap.poscar', './POSCAR')
            vasp.WriteInp('fbrlx', pawpath)
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 3:
            fileop.copyfile('../s2/CONTCAR', './POSCAR')
            vasp.WriteInp('scf', pawpath)
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 4:
            fileop.copyfile('../s3/CONTCAR', './POSCAR')
            fileop.copyfile('../s3/CHGCAR', './CHGCAR')
            vasp.WriteInp('nscf', pawpath)
            fileop.copyfile('../s0/kpath_unwrap.poscar', './KPOINTS')
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 5:
            fileop.copyfile('../s0/std_unwrap.poscar', './')
            fileop.copyfile('./std_unwrap.poscar', './POSCAR')
            vasp.WriteInp('arlx', pawpath)
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 6:
            fileop.copyfile('../s5/CONTCAR', './POSCAR')
            vasp.WriteInp('scf', pawpath)
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 7:
            fileop.copyfile('../s6/CONTCAR', './POSCAR')
            fileop.copyfile('../s6/CHGCAR', './CHGCAR')
            vasp.WriteInp('nscf', pawpath)
            fileop.copyfile('../s0/kpath_unwrap.poscar', './KPOINTS')
            shellcall.runvasp()
            os.system('rm -rf CHG')

        elif stepnumber == 8:
            fileop.copyfile('../s0/dimers.pkl', './dimers.pkl')
            gauss.WriteInp(jobtype='dimers', rawfile='dimers.pkl', addkeywds='')
            for i in glob.glob('*.com'):
                shellcall.rungauss(i, 'g16')
            for i in glob.glob('*.chk'):
                os.system('formchk ' + i)

        elif stepnumber == 9:
            fileop.copyfile('../s0/mol0.xyz', './mol0.xyz')
            gauss.WriteInp(jobtype='reorg_opt', rawfile='mol0.xyz', addkeywds='')
            for i in glob.glob('*.com'):
                shellcall.rungauss(i, 'g16')

        elif stepnumber == 10:
            fileop.copyfile('../s9/nopt.chk', './')
            fileop.copyfile('../s9/copt.chk', './')
            fileop.copyfile('../s9/aopt.chk', './')
            gauss.WriteInp(jobtype='reorg_sp', rawfile='', addkeywds='')
            for i in glob.glob('*.com'):
                shellcall.rungauss(i, 'g16')
            for i in glob.glob('*.chk'):
                os.system('formchk ' + i)

        elif stepnumber == 11:
            fileop.copyfile('../s9/nopt.chk', './')
            gauss.WriteInp(jobtype='tdsinglet', rawfile='', addkeywds='')
            for i in glob.glob('*.com'):
                shellcall.rungauss(i, 'g16')
            for i in glob.glob('*.chk'):
                os.system('formchk ' + i)

        elif stepnumber == 12:

            fileop.copyfile('../s7/CONTCAR', './POSCAR')

            ca = configuration.Configuration('POSCAR')
            ca.unwrap_structure.to(fmt='poscar', filename='unwrap.poscar')
            fileop.write_obj(ca.dimers, 'dimers.pkl')
            gauss.WriteInp(jobtype='dimers', rawfile='dimers.pkl', addkeywds='')
            for i in glob.glob('*.com'):
                shellcall.rungauss(i, 'g16')
            for i in glob.glob('*.chk'):
                os.system('formchk ' + i)

        os.chdir(whereami)
