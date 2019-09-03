import warnings
import sys
import os
import glob
import shellcall
import fileop
import pymatop, mathop
import settings, configuration
import omol
import reciprocal, coupling, reorgnization
from pymatgen.io.gaussian import GaussianInput
import gaussplot
from pymatgen.core import Molecule
import numpy as np
import core


class Analyzer:
    """
    perform analysis for a crystal, based on its unique name in the database, at the analysis folder
    """

    def __init__(self, name, tankload):

        print('------ analyzer launched')
        print('------ analyzer for ' + name)

        self.name = name
        self.cifname = name + '.cif'
        self.tankload = tankload

        self.globalpath = settings.SYSPATH

        self.runpath = self.globalpath['RUN_PATH'] + '/' + self.name + '/'  # the path where calculations are done
        self.anapath = self.globalpath['ANA_PATH'] + '/' + self.name + '/'  # where the analysis will be done

        fileop.createdir(self.anapath)

        # check path in settings
        if not fileop.check_path([self.globalpath[k] for k in self.globalpath.keys()]):
            warnings.warn('check system path!')
            sys.exit(1)

        rocket = core.Rocket(self.name, tankload=self.tankload)
        if not rocket.rocket_status:
            warnings.warn('rocket can still fly')
            sys.exit(1)

        self.data = dict(name=self.name)
        for i in self.tankload:
            self.task_analysis(int(i))
            print('------ task done', i)

        print(self.data)
        print('------ analyzer done for ' + self.name)

    def task_analysis(self, tanknumber):

        whereami = os.getcwd()

        os.chdir(self.anapath)
        if tanknumber == 0:
            self.data['config_per_cell'] = len(glob.glob(self.runpath + '/s0/' + 'config*.cif'))
            self.data['mol_per_cell'] = len(glob.glob(self.runpath + '/s0/' + 'mol*.poscar'))
            fileop.copyfile(self.runpath + '/s0/slip_vis.xyz', './')
            fileop.copyfile(self.runpath + '/s0/slip_data.txt', './')
            configs = glob.glob(self.runpath + '/s0/config*.cif')
            for i in configs:
                fileop.copyfile(i, self.anapath)
            fileop.copyfile(self.runpath + '/s0/mol0.xyz', self.anapath)
            self.data['molecule'] = omol.OMol.from_xyz('mol0.xyz').data
            # ca = configuration.Configuration('config0.cif')
            # fileop.write_obj(ca.dimers, 'dimers.pkl')

        elif tanknumber == 1:
            fileop.copyfile(self.runpath + '/s1/OSZICAR', './OSZICAR_s1')
            self.data['mol_e'] = shellcall.oszicar_energy('OSZICAR_s1')

        elif tanknumber == 3:
            fileop.copyfile(self.runpath + '/s3/OSZICAR', './OSZICAR_s3')
            self.data['fb_e'] = shellcall.oszicar_energy('OSZICAR_s3')

        elif tanknumber == 4:
            fileop.copyfile(self.runpath + '/s4/OUTCAR', self.anapath)
            fileop.copyfile(self.runpath + '/s4/EIGENVAL', self.anapath)
            fileop.copyfile(self.runpath + '/s4/KPOINTS', self.anapath)
            ra = reciprocal.ReciprocalAnalyzer('fb', usefitdata=True)
            os.remove('OUTCAR')
            os.remove('EIGENVAL')
            os.remove('KPOINTS')
            self.data['fb_gap'] = ra.gap
            self.data['fb_e_em'] = ra.e_em
            self.data['fb_h_em'] = ra.h_em

        elif tanknumber == 6:
            fileop.copyfile(self.runpath + '/s6/OSZICAR', './OSZICAR_s6')
            self.data['ar_e'] = shellcall.oszicar_energy('OSZICAR_s6')
            self.data['cohesive'] = abs((self.data['ar_e'] - self.data['mol_per_cell'] * self.data['mol_e']) /
                                        self.data['mol_per_cell'])

        elif tanknumber == 7:
            fileop.copyfile(self.runpath + '/s7/OUTCAR', self.anapath)
            fileop.copyfile(self.runpath + '/s7/EIGENVAL', self.anapath)
            fileop.copyfile(self.runpath + '/s7/KPOINTS', self.anapath)
            ra = reciprocal.ReciprocalAnalyzer('ar', usefitdata=True)
            os.remove('OUTCAR')
            os.remove('EIGENVAL')
            os.remove('KPOINTS')
            self.data['ar_gap'] = ra.gap
            self.data['ar_e_em'] = ra.e_em
            self.data['ar_h_em'] = ra.h_em

        elif tanknumber == 8:
            fns = glob.glob(self.runpath + '/s8/*.fchk') + glob.glob(self.runpath + '/s8/*_A_B.log')
            fns += glob.glob(self.runpath + '/s8/*_A_B.com')
            for fn in fns:
                fileop.copyfile(fn, self.anapath)

            # for aug 1st
            # dimers = fileop.read_obj('dimers.pkl')
            dimers = []
            for i in glob.glob('*_A_B.com'):
                label = i[5:-8]
                gi = GaussianInput.from_file(i)
                grps = pymatop.group_close_sites(gi.molecule.sites)
                refmol = configuration.OrganicMolecule(Molecule.from_sites(grps[0]))
                othmol = configuration.OrganicMolecule(Molecule.from_sites(grps[1]))
                pslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_p)
                oslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_o)
                qslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_q)
                pslipration = pslip / refmol.bone.lp
                qslipration = pslip / refmol.bone.lq
                h_angle = mathop.angle_btw(refmol.bone.uv_o, othmol.bone.uv_o, output='degree')
                p_angle = mathop.angle_btw(refmol.bone.uv_p, othmol.bone.uv_p, output='degree')
                jmol_command = '# draw arrow {' + ','.join([str(round(ii, 2)) for ii in refmol.bone.com]) + \
                               '} {' '.'.join([str(round(ii, 2)) for ii in othmol.bone.com]) + '} \r\n'
                s ={
                    'pslip':pslip,
                    'qslip':qslip,
                    'oslip':oslip,
                    'h_angle':h_angle,
                    'p_angle':p_angle,
                    'jmol':jmol_command,
                }
                d = configuration.Dimer(s, refmol, othmol, label)
                dimers.append(d)
            # end for aug 1st

            couplings = []
            for dimer in dimers:
                elec_coupling, hole_coupling = self.coupling('dimer' + dimer.label)
                couplings.append([dimer.label, elec_coupling, hole_coupling, dimer.slip])
            couplings = sorted(couplings, key=lambda x: x[1])
            self.data['max_ecoupling'] = couplings[-1][1]
            self.data['max_ecoupling_label'] = couplings[-1][0]
            self.data['max_e_dimer_pslip'] = couplings[-1][3]['pslip']
            self.data['max_e_dimer_qslip'] = couplings[-1][3]['qslip']
            self.data['max_e_dimer_oslip'] = couplings[-1][3]['oslip']
            self.data['max_e_dimer_h_angle'] = couplings[-1][3]['h_angle']
            self.data['max_e_dimer_p_angle'] = couplings[-1][3]['p_angle']
            self.data['max_e_dimer_jmol'] = couplings[-1][3]['jmol']
            couplings = sorted(couplings, key=lambda x: x[2])
            self.data['max_hcoupling'] = couplings[-1][2]
            self.data['max_hcoupling_label'] = couplings[-1][0]
            self.data['max_h_dimer_pslip'] = couplings[-1][3]['pslip']
            self.data['max_h_dimer_qslip'] = couplings[-1][3]['qslip']
            self.data['max_h_dimer_oslip'] = couplings[-1][3]['oslip']
            self.data['max_h_dimer_h_angle'] = couplings[-1][3]['h_angle']
            self.data['max_h_dimer_p_angle'] = couplings[-1][3]['p_angle']
            self.data['max_h_dimer_jmol'] = couplings[-1][3]['jmol']

        elif tanknumber == 9:
            fns = glob.glob(self.runpath + '/s9/*.chk') + glob.glob(self.runpath + '/s9/*.log')
            for fn in fns:
                fileop.copyfile(fn, self.anapath)
                if fn[-4:] == '.chk':
                    os.system('formchk ' + fn.split('/')[-1])
            with open('nopt.fchk', 'r') as f:
                fchklines = f.readlines()
            homo, lumo = coupling.ElectronicCoupling.g09_get_homo_lumo(fchklines)
            nbf = coupling.ElectronicCoupling.g09_get_basis_functions(fchklines)
            moeigens = coupling.ElectronicCoupling.g09_get_alpha_orbital_energies(fchklines, nbf)
            self.data['homo_m1'] = moeigens[0, homo-2] * 27.2114
            self.data['homo'] = moeigens[0, homo-1] * 27.2114
            self.data['lumo'] = moeigens[0, lumo-1] * 27.2114
            self.data['lumo_p1'] = moeigens[0, lumo] * 27.2114
            copt_energy = reorgnization.ReorganizationEnergy.extract_energy('copt.log')
            aopt_energy = reorgnization.ReorganizationEnergy.extract_energy('aopt.log')
            nopt_energy = reorgnization.ReorganizationEnergy.extract_energy('nopt.log')
            self.data['aip'] = abs(copt_energy - nopt_energy) * 27.2114
            self.data['aea'] = abs(aopt_energy - nopt_energy) * 27.2114

        elif tanknumber == 10:
            fns = glob.glob(self.runpath + '/s9/*.log') + glob.glob(self.runpath + '/s10/*.log')
            for fn in fns:
                fileop.copyfile(fn, self.anapath)
            re = reorgnization.ReorganizationEnergy()
            self.data['reorg_h'] = re.hole_reorg
            self.data['reorg_e'] = re.elec_reorg

        elif tanknumber == 11:
            fileop.copyfile(self.runpath + '/s11/tdsinglet.log', self.anapath)
            gaussplot.pltuvvis(fn='tdsinglet.log')
            excites = self.td_parser('tdsinglet.log')
            self.data['lambda1'] = excites[0][0]
            self.data['lambda2'] = excites[1][0]
            self.data['lambda3'] = excites[2][0]
            self.data['os1'] = excites[0][1]
            self.data['os2'] = excites[1][1]
            self.data['os3'] = excites[2][1]
            self.data['dom_trans1'] = excites[0][2]
            self.data['dom_trans2'] = excites[1][2]
            self.data['dom_trans3'] = excites[2][2]
            self.data['dom_trans1_percent'] = excites[0][3]
            self.data['dom_trans2_percent'] = excites[1][3]
            self.data['dom_trans3_percent'] = excites[2][3]

        elif tanknumber == 12:
            fns = glob.glob(self.runpath + '/s12/*.fchk') + glob.glob(self.runpath + '/s12/*_A_B.log')
            fns += glob.glob(self.runpath + '/s12/*_A_B.com')
            for fn in fns:
                fileop.copyfile(fn, self.anapath)

            # for aug 1st
            # dimers = fileop.read_obj('dimers.pkl')
            dimers = []
            for i in glob.glob('*_A_B.com'):
                label = i[5:-8]
                gi = GaussianInput.from_file(i)
                grps = pymatop.group_close_sites(gi.molecule.sites)
                refmol = configuration.OrganicMolecule(Molecule.from_sites(grps[0]))
                othmol = configuration.OrganicMolecule(Molecule.from_sites(grps[1]))
                pslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_p)
                oslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_o)
                qslip = np.dot(refmol.bone.com - othmol.bone.com, refmol.bone.uv_q)
                pslipration = pslip / refmol.bone.lp
                qslipration = pslip / refmol.bone.lq
                h_angle = mathop.angle_btw(refmol.bone.uv_o, othmol.bone.uv_o, output='degree')
                p_angle = mathop.angle_btw(refmol.bone.uv_p, othmol.bone.uv_p, output='degree')
                jmol_command = '# draw arrow {' + ','.join([str(round(ii, 2)) for ii in refmol.bone.com]) + \
                               '} {' '.'.join([str(round(ii, 2)) for ii in othmol.bone.com]) + '} \r\n'
                s ={
                    'pslip':pslip,
                    'qslip':qslip,
                    'oslip':oslip,
                    'h_angle':h_angle,
                    'p_angle':p_angle,
                    'jmol':jmol_command,
                }
                d = configuration.Dimer(s, refmol, othmol, label)
                dimers.append(d)
                print(oslip)
            # end for aug 1st

            couplings = []
            for dimer in dimers:
                elec_coupling, hole_coupling = self.coupling('dimer' + dimer.label)
                couplings.append([dimer.label, elec_coupling, hole_coupling, dimer.slip])
            couplings = sorted(couplings, key=lambda x: x[1])
            self.data['max_ecoupling'] = couplings[-1][1]
            self.data['max_ecoupling_label'] = couplings[-1][0]
            self.data['max_e_dimer_pslip'] = couplings[-1][3]['pslip']
            self.data['max_e_dimer_qslip'] = couplings[-1][3]['qslip']
            self.data['max_e_dimer_oslip'] = couplings[-1][3]['oslip']
            self.data['max_e_dimer_h_angle'] = couplings[-1][3]['h_angle']
            self.data['max_e_dimer_p_angle'] = couplings[-1][3]['p_angle']
            self.data['max_e_dimer_jmol'] = couplings[-1][3]['jmol']
            couplings = sorted(couplings, key=lambda x: x[2])
            self.data['max_hcoupling'] = couplings[-1][2]
            self.data['max_hcoupling_label'] = couplings[-1][0]
            self.data['max_h_dimer_pslip'] = couplings[-1][3]['pslip']
            self.data['max_h_dimer_qslip'] = couplings[-1][3]['qslip']
            self.data['max_h_dimer_oslip'] = couplings[-1][3]['oslip']
            self.data['max_h_dimer_h_angle'] = couplings[-1][3]['h_angle']
            self.data['max_h_dimer_p_angle'] = couplings[-1][3]['p_angle']
            self.data['max_h_dimer_jmol'] = couplings[-1][3]['jmol']

        os.chdir(whereami)

    @staticmethod
    def coupling(dimername):
        if not all([os.path.isfile(dimername + suffix) for suffix in ['_A.fchk', '_B.fchk',
                                                                      '_A_B.fchk', '_A_B.log']]):
            sys.exit(1)

        ec = coupling.ElectronicCoupling(dimername)

        je12_elec, je12_hole = [0.0, 0.0]
        for entry in ec.output:
            if entry[0] == ec.ahomo and entry[1] == ec.bhomo:
                je12_hole = entry[8]
            elif entry[0] == ec.alumo and entry[1] == ec.blumo:
                je12_elec = entry[8]
        return abs(je12_elec), abs(je12_hole)

    @staticmethod
    def td_parser(log):
        with open(log, 'r') as f:
            lines = f.readlines()

        excites = np.empty(3, dtype=[('lambdanm', 'f8'), ('os', 'f8'), ('dom.trans.', 'U16'),
                                     ('dom.trans.percent', 'f8')])

        homo, lumo = (0, 0)
        for line in lines:
            if 'alpha electrons' in line:
                num_alpha_electrons = int(line.strip().split()[0])
                homo = num_alpha_electrons  # starts from 1
                lumo = homo + 1
            elif 'Excited State   1' in line:
                items = line.strip().split()
                lambdanm = float(items[6])
                osci = float(items[8].split('=')[1])
                line_num = lines.index(line)
                next_line = 1
                trans = []
                while len(lines[line_num + next_line].strip().split()) in [3, 4] and \
                        '->' in lines[line_num + next_line]:
                    line_items = lines[line_num + next_line].strip().split()
                    line_items = [i.strip('->') for i in line_items]
                    line_items = [i for i in line_items if i != '']
                    frommo = int(line_items[0])
                    frommo_diff = homo - frommo
                    tomo = int(line_items[1])
                    tomo_diff = tomo - lumo
                    transstr = 'HOMO-' + str(frommo_diff) + ' -- ' + 'LUMO+' + str(tomo_diff)
                    percent = float(line_items[2])
                    trans.append([transstr, percent])
                    next_line += 1
                domtrans = sorted(trans, key=lambda x: x[1], reverse=True)[0]
                excites[0] = (lambdanm, osci, domtrans[0], domtrans[1])
            elif 'Excited State   2' in line:
                items = line.strip().split()
                lambdanm = float(items[6])
                osci = float(items[8].split('=')[1])
                line_num = lines.index(line)
                next_line = 1
                trans = []
                while len(lines[line_num + next_line].strip().split()) in [3, 4] and \
                        '->' in lines[line_num + next_line]:
                    line_items = lines[line_num + next_line].strip().split()
                    line_items = [i.strip('->') for i in line_items]
                    line_items = [i for i in line_items if i != '']
                    frommo = int(line_items[0])
                    frommo_diff = homo - frommo
                    tomo = int(line_items[1])
                    tomo_diff = tomo - lumo
                    transstr = 'HOMO-' + str(frommo_diff) + ' -- ' + 'LUMO+' + str(tomo_diff)
                    percent = float(line_items[2])
                    trans.append([transstr, percent])
                    next_line += 1
                domtrans = sorted(trans, key=lambda x: x[1], reverse=True)[0]
                excites[1] = (lambdanm, osci, domtrans[0], domtrans[1])
            elif 'Excited State   3' in line:
                items = line.strip().split()
                lambdanm = float(items[6])
                osci = float(items[8].split('=')[1])
                line_num = lines.index(line)
                next_line = 1
                trans = []
                while len(lines[line_num + next_line].strip().split()) in [3, 4] and \
                        '->' in lines[line_num + next_line]:
                    line_items = lines[line_num + next_line].strip().split()
                    line_items = [i.strip('->') for i in line_items]
                    line_items = [i for i in line_items if i != '']
                    frommo = int(line_items[0])
                    frommo_diff = homo - frommo
                    tomo = int(line_items[1])
                    tomo_diff = tomo - lumo
                    transstr = 'HOMO-' + str(frommo_diff) + ' -- ' + 'LUMO+' + str(tomo_diff)
                    percent = float(line_items[2])
                    trans.append([transstr, percent])
                    next_line += 1
                domtrans = sorted(trans, key=lambda x: x[1], reverse=True)[0]
                excites[2] = (lambdanm, osci, domtrans[0], domtrans[1])
        return excites

    # def config_analysis(self):
    #     whereami = os.getcwd()
    #
    #     os.chdir(self.anapath)
    #     self.data['num_of_config'] = len(glob.glob(self.runpath + '/s0/' + 'config*.cif'))
    #     self.data['num_of_mol'] = len(glob.glob(self.runpath + '/s0/' + 'mol*.poscar'))
    #     sysrwc.copyfile(self.runpath + '/s0/config0.cif', self.anapath)
    #     config = Analysis.schema.Configuration('config0.cif')
    #
    #     config.write_slip()
    #     dimer = sorted(config.dimers, key=lambda x: x.slip['pslip_ratio'])[0]
    #     self.data['dimer'] = dimer
    #
    #     omols = config.organic_mols
    #     if len(list(set([i.composition for i in omols]))) == 1 and len(list(set([i.composition for i
    #                                                                              in omols[0].side_chains]))) == 1:
    #         self.data['side_chain_sum_v_norm'] = omols[0].side_chains[0].sum_v_norm
    #         self.data['side_chain_shortest'] = omols[0].side_chains[0].shortest
    #         self.data['side_chain_longest'] = omols[0].side_chains[0].longest
    #         self.data['electron_density'] = omols[0].side_chains[0].electron_density
    #         self.data['dis_list_std'] = omols[0].side_chains[0].dis_list_std
    #         self.data['dis_list_mean'] = omols[0].side_chains[0].dis_list_mean
    #     else:
    #         warnings.warn('not single crystal or side chains are different')
    #         sys.exit(1)
    #
    #     os.chdir(whereami)
    #
    #
    # def electronic_analysis(self, jtype):
    #     whereami = os.getcwd()
    #
    #     os.chdir(self.anapath)
    #     if jtype == 'fb':
    #         sysrwc.copyfile(self.runpath + '/s4/OUTCAR', self.anapath)
    #         sysrwc.copyfile(self.runpath + '/s4/EIGENVAL', self.anapath)
    #         sysrwc.copyfile(self.runpath + '/s4/KPOINTS', self.anapath)
    #     elif jtype == 'ar':
    #         sysrwc.copyfile(self.runpath + '/s7/OUTCAR', self.anapath)
    #         sysrwc.copyfile(self.runpath + '/s7/EIGENVAL', self.anapath)
    #         sysrwc.copyfile(self.runpath + '/s7/KPOINTS', self.anapath)
    #
    #     ra = Analysis.reciprocal.ReciprocalAnalyzer(jtype)
    #
    #     os.remove('OUTCAR')
    #     os.remove('EIGENVAL')
    #     os.remove('KPOINTS')
    #     self.data[jtype + '_e_em'] = ra.e_em
    #     self.data[jtype + '_h_em'] = ra.h_em
    #
    #     os.chdir(whereami)
