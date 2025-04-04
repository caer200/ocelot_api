"""
originally written by Sean M. Ryno, Cheng Zhong, Haitao Sun, see /legacy/aw_tuning.py
for a certain mol, get tuned w
default keywords for gaussian
'scf=(xqc,fermi,noincfock,ndamp=35,conver=6,vshift=500,novaracc)'
dev-ing
"""
import warnings
import datetime as dt
import os
import random
import re
from subprocess import Popen

from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput
from pymatgen.io.gaussian import GaussianOutput
from scipy import optimize

SUGGESTED_route_parameters = {
    'scf': {
        'xqc': '',
        'fermi': '',
        'noincfock': '',
        'novaracc': '',
        'ndamp': '35',
        'conver': '6',
        'vshift': '500'
    }
}


def gauss_in_gen(name, mol, type, func, basis, charge, spin, route_params, link0_params):
    mol_in = GaussianInput(mol, charge=charge, spin_multiplicity=spin, functional=func, basis_set=basis,
                           route_parameters=route_params, link0_parameters=link0_params, )
    fn = '{}_{}.gjf'.format(name, type)
    mol_in.write_file(filename=fn, cart_coords=True)
    return fn


class WtuningJob:

    def __init__(self, func='uLC-wHPBE', basis='def2tzvp', name='', nproc=16, mem=50,
                 n_charge=0, n_spin=1, wdir='./', rps=SUGGESTED_route_parameters, scheme='Jh', wbmin=0.05, wbmax=0.5,
                 abmin=0.05, abmax=0.5, gauss_cmd='g16'):

        self.name = name
        self.start_time = dt.datetime.now()
        self.func = func
        self.basis = basis
        self.nproc = nproc
        self.mem = mem
        # self.route_params = route_params
        self.link0_params = {'%nproc': nproc, '%mem': str(mem) + 'GB'}
        self.gauss_cmd = gauss_cmd

        # tuning params
        self.wbounds = (wbmin, wbmax)
        self.abounds = (abmin, abmax)
        self.conver = 5
        self.scheme = scheme  # Jl Jh Jn O2
        self.route_params = rps
        self.wdir = wdir
        self.n_charge = n_charge
        self.n_spin = n_spin
        self.c_charge = int(self.n_charge) + 1
        self.a_charge = int(self.n_charge) - 1
        if self.n_spin == 1:
            self.c_spin = 2
            self.a_spin = 2
        else:
            self.c_spin = self.n_spin - 1
            self.a_spin = self.n_spin - 1

        # tuning var
        self.omega = random.uniform(*self.wbounds)
        self.alpha = random.uniform(*self.abounds)
        self.beta = 1 - self.alpha
        # self.super_omega = self.omega
        self.cycle = 0
        self.ocycle = 0
        self.mol = None

    @staticmethod
    def gauss_run(fn, gcmd):
        job = Popen([gcmd, fn])
        job.wait()
        return fn.split('.')[0] + '.log'

    def omega_tune(self, dis=3, tol=1e-04, deltaomega=0.2):

        if self.ocycle > 0:
            bounds = (self.omega - deltaomega, self.omega + deltaomega)
        else:
            bounds = self.wbounds

        self.cycle = 0
        whereami = os.getcwd()
        os.chdir(self.wdir)
        new_dir = 'tuning_' + self.name + '_ocycle' + str(self.ocycle)
        new_dir = '{}/{}'.format(self.wdir, new_dir)
        # os.mkdir(new_dir)
        os.makedirs(new_dir)
        os.chdir(new_dir)
        omega_opt, C_opt, err, num = optimize.fminbound(self.omega_wtune, bounds[0],
                                                        bounds[1], disp=dis, xtol=tol,
                                                        full_output=True)
        os.chdir(whereami)
        return omega_opt, C_opt, num

    def alpha_tune(self, dis=3, tol=1e-04, deltaalpha=0.2):
        if self.ocycle > 0:
            bounds = (self.alpha - deltaalpha, self.alpha + deltaalpha)
        else:
            bounds = self.abounds

        self.cycle = 0
        whereami = os.getcwd()
        os.chdir(self.wdir)
        new_dir = 'tuning_' + self.name + '_acycle' + str(self.ocycle)
        new_dir = '{}/{}'.format(self.wdir, new_dir)
        # os.mkdir(new_dir)
        os.makedirs(new_dir)
        os.chdir(new_dir)
        alpha_opt, C_opt, err, num = optimize.fminbound(self.alpha_atune, bounds[0],
                                                        bounds[1], disp=dis, xtol=tol,
                                                        full_output=True)
        os.chdir(whereami)
        return alpha_opt, C_opt, num

    # iop_route_param = 'iop(3/107={}, 3/108={})'.format(self.omega_iopstr, self.omega_iopstr)

    @property
    def iop_route_params(self):
        return 'iop(3/107={}, 3/108={})'.format(self.omega_iopstr, self.omega_iopstr)

    def geo_opt(self, rps={}):
        rps[self.iop_route_params] = ''
        rps['opt'] = ''
        fnin = gauss_in_gen(
            name=self.name,
            mol=self.mol,
            route_params=rps,
            type='opt_{}'.format(self.ocycle),
            func=self.func,
            basis=self.basis,
            charge=self.n_charge,
            spin=self.n_spin,
            link0_params=self.link0_params
        )
        fnout = self.gauss_run(fnin, gcmd=self.gauss_cmd)
        return Molecule.from_file(fnout)  # Return mol

    def wtuning_cycle(self, eps=0.01,
                      dis=3, tol=1e-3, max_cycles=5):
        while True:
            oldomega = self.omega
            self.log_add({'Super Cycle': self.ocycle, 'cycle init with w': self.omega})
            omega = self.omega_tune(dis=dis, tol=tol)[0]
            self.log_add({'new': omega, 'old': oldomega})
            if abs(omega - oldomega) <= eps and self.ocycle > 0:
                self.omega = omega
                break
            elif self.ocycle >= max_cycles:
                warnings.warn('tuning cycle went over max cycles')
                break
            self.omega = omega
            self.mol = self.geo_opt()
            self.ocycle += 1

    def atuning_cycle(self, eps=0.01, dis=3, tol=1e-3, max_cycles=5):
        while True:
            oldalpha = self.alpha
            self.log_add(
                {'Super Cycle': self.ocycle, 'a': self.alpha, 'Elapsed Time': dt.datetime.now() - self.start_time})
            alpha = self.alpha_tune(dis=dis, tol=tol)[0]
            self.log_add({'new': alpha, 'old': oldalpha})
            if abs(alpha - oldalpha) <= eps and self.ocycle > 0:
                self.alpha = alpha
                break
            elif self.ocycle >= max_cycles:
                warnings.warn('tuning cycle went over max cycles')
                break
            self.alpha = alpha
            self.mol = self.geo_opt()
            self.ocycle += 1

    @staticmethod
    def omega_extract(fn):
        """
        Pull important data from Log Files after run
        :param fn: Filename of Log File to be read
        :return: Energies of Molecule, HOMO and LUMO in eV
        """
        with open(fn, 'r') as f:
            LogFile = f.readlines()
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
            return None
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

    @property
    def omega_iopstr(self):
        return str(int(float(round(self.omega * 1e4)))).zfill(5) + '00000'

    @property
    def alpha_iopstr(self):
        return str(int(float(round(self.alpha * 1e4)))).zfill(5) + '00000'

    @property
    def beta_iopstr(self):
        return str(int(float(round(self.beta * 1e4)))).zfill(5) + '00000'

    def omega_format(self):
        """
        Format IOp strings based off of omega value, then update route parameters dictionary and input filenames
        """
        # self.omega = round(self.omega, self.conver)
        try:
            self.route_params.pop(self.iop_route_paramw)
        except:
            pass
        self.iop_route_paramw = 'iop(3/107={}, 3/108={})'.format(self.omega_iopstr, self.omega_iopstr)
        self.route_params[self.iop_route_paramw] = ''

        self.n_input_fn = 'n_omega_{}_{}.com'.format(self.omega_iopstr, self.cycle)
        self.c_input_fn = 'c_omega_{}_{}.com'.format(self.omega_iopstr, self.cycle)
        self.a_input_fn = 'a_omega_{}_{}.com'.format(self.omega_iopstr, self.cycle)

    def alpha_format(self):
        try:
            self.route_params.pop(self.iop_route_parama)
        except:
            pass
        self.iop_route_parama = 'iop(3/130={}, 3/131={}, 3/119={}, 3/120={})'.format(self.alpha_iopstr,
                                                                                     self.alpha_iopstr,
                                                                                     self.beta_iopstr, self.beta_iopstr)
        self.route_params[self.iop_route_parama] = ''

        self.n_input_fn = 'n_alpha_{}_{}.com'.format(self.alpha_iopstr, self.cycle)
        self.c_input_fn = 'c_alpha_{}_{}.com'.format(self.alpha_iopstr, self.cycle)
        self.a_input_fn = 'a_alpha_{}_{}.com'.format(self.alpha_iopstr, self.cycle)

    def omega_gauss_do(self):
        """
        Run Gaussian in subprocess and wait for termination. Extract data from output when done
        """
        self.n_in_fn = gauss_in_gen(mol=self.mol, charge=self.n_charge, spin=self.n_spin,
                                    type='tune_n' + str(self.cycle), basis=self.basis, route_params=self.route_params,
                                    link0_params=self.link0_params, name=self.name, func=self.func)
        self.n_out_fn = self.gauss_run(self.n_in_fn, self.gauss_cmd)
        self.n_e, self.n_homo, self.n_lumo = self.omega_extract(self.n_out_fn)
        if self.scheme == 'Jh':
            self.c_in_fn = gauss_in_gen(mol=self.mol, charge=self.c_charge, spin=self.c_spin,
                                        type='tune_c' + str(self.cycle), basis=self.basis,
                                        route_params=self.route_params, link0_params=self.link0_params, name=self.name,
                                        func=self.func)
            self.c_out_fn = self.gauss_run(self.c_in_fn, self.gauss_cmd)
            self.c_e, self.c_homo, self.c_lumo = self.omega_extract(self.c_out_fn)
        elif self.scheme == 'Jl':
            self.a_in_fn = gauss_in_gen(mol=self.mol, charge=self.a_charge, spin=self.a_spin,
                                        type='tune_a' + str(self.cycle), basis=self.basis,
                                        route_params=self.route_params, link0_params=self.link0_params, name=self.name,
                                        func=self.func)
            self.a_out_fn = self.gauss_run(self.a_in_fn, self.gauss_cmd)
            self.a_e, self.a_homo, self.a_lumo = self.omega_extract(self.a_out_fn)
        else:

            self.c_in_fn = gauss_in_gen(mol=self.mol, charge=self.c_charge, spin=self.c_spin,
                                        type='tune_c' + str(self.cycle), basis=self.basis,
                                        route_params=self.route_params, link0_params=self.link0_params, name=self.name,
                                        func=self.func)
            self.a_in_fn = gauss_in_gen(mol=self.mol, charge=self.a_charge, spin=self.a_spin,
                                        type='tune_a' + str(self.cycle), basis=self.basis,
                                        route_params=self.route_params, link0_params=self.link0_params, name=self.name,
                                        func=self.func)

            self.c_out_fn = self.gauss_run(self.c_in_fn, self.gauss_cmd)
            self.a_out_fn = self.gauss_run(self.a_in_fn, self.gauss_cmd)

            self.c_e, self.c_homo, self.c_lumo = self.omega_extract(self.c_out_fn)
            self.a_e, self.a_homo, self.a_lumo = self.omega_extract(self.a_out_fn)
        self.log_add({'c_e':self.c_e,'c_homo':self.c_homo,'c_lumo':self.c_lumo,
                      'n_e':self.n_e,'n_homo':self.n_homo,'n_lumo':self.n_lumo,
                      })
        self.cycle += 1

    def omega_FindC(self):
        """
        Calculate scheme value from extracted data
        Set the optimization criterion (the value to minimize),
        available options are:
        J2---((HOMO-IP)^2+(A_HOMO+EA)^2)
        Jh---(HOMO-IP)
        Jl---(LUMO+EA)
        Jn2---((HOMO-IP)^2+(LUMO+EA)^2)
        O2---((A_HOMO-LUMO)^2+(C_LUMO-HOMO)^2)
        :return: Jn, Jl, Jh value (depending on scheme)
        """
        if self.scheme == 'Jh':
            IP = self.n_e - self.c_e
            Jh = abs(self.n_homo - IP)
        elif self.scheme == 'Jl':
            EA = self.n_e - self.a_e
            Jl = abs(self.n_lumo + EA)
        elif self.scheme == 'J2':
            IP = self.n_e - self.c_e
            EA = self.n_e - self.a_e
            J2 = (self.n_homo - IP) ** 2 + (self.a_homo + EA) ** 2
        elif self.scheme == 'Jn2':
            IP = self.n_e - self.c_e
            Jh = abs(self.n_homo - IP)
            EA = self.n_e - self.a_e
            Jl = abs(self.n_lumo + EA)
            Jn2 = Jh ** 2 + Jl ** 2
        elif self.scheme == 'O2':
            O2 = (self.n_homo - self.c_lumo) ** 2 + (self.n_lumo - self.a_homo) ** 2
        C = eval(self.scheme)
        self.log_add({'C':C})
        return C

    def omega_wtune(self, omega_in):
        """
        :param omega_in: Value for 'fminbound' function to pass scalars into
        :return: Scheme value from FindC
        """
        self.omega = omega_in
        self.omega_format()
        self.log_add({'Subcycle': self.cycle, 'w': self.omega})
        self.omega_gauss_do()
        return self.omega_FindC()

    def alpha_atune(self, alpha_in):
        self.alpha = alpha_in
        self.beta = 1 - self.alpha
        self.alpha_format()
        self.log_add({'Subcycle': self.cycle, 'a': self.alpha})
        self.omega_gauss_do()
        return self.omega_FindC()

    @staticmethod
    def extract_energy(fn):
        return GaussianOutput(fn).final_energy

    @staticmethod
    def extract_abs(fn):
        return GaussianOutput(fn).read_excitation_energies()

    # I/O
    def to_dict(self):
        pass

    def to_file(self):
        name = self.name + '.out'
        with open(name, 'w') as fn:
            fn.write(self.name + '\n')
            fn.write(self.func + '/' + self.basis + '\n')
            fn.write('omega: ' + str(self.omega) + '\n')

    def log_add(self,s):
        with open('output.log','a') as fn:
            for key in s:
                string = str(key)+': {}\n'.format(s[key])
                fn.write(string)

