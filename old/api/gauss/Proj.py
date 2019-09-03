import datetime
import multiprocessing
import sys

<<<<<<< HEAD
from api.schema.omol import Mol
=======
from api.schema.Mol import Mol
>>>>>>> 88b3e96e2d37f403133c7972e5e02f735a27c716
from api.schema.msite import Site


class Gaussproj:

    def __init__(self, mol, jobname,
                 hashline='#p wb97xd/def2svp scf=(xqc,fermi,noincfock,ndamp=35,conver=8,vshift=500,novaracc) int=(acc2e=16)',
                 mem='20GB', nproc=multiprocessing.cpu_count()):
        """

        :param mol:
        :param jobtype:
        :param jobname:
        :param hashline:
        :param mem:
        :param nproc:
        """
        self.mol = mol
        self.jobname = jobname
        self.hashline = hashline
        self.addkwds = []
        self.level = ''
        self.jobtype = []
        for kw in self.hashline.strip().split():
            if '/' in kw and 'iop' not in kw.lower():
                self.level = kw
            else:
                self.addkwds.append(kw)
            if 'opt' in kw.lower():
                self.jobtype.append(kw)
            if 'freq' in kw.lower():
                self.jobtype.append(kw)
            if 'nmr' in kw.lower():
                self.jobtype.append(kw)
            if 'td' in kw.lower():
                self.jobtype.append(kw)
        self.mem = mem
        self.nproc = nproc

    @classmethod
    def from_gjf(cls, gjf, mem='20GB', nproc=multiprocessing.cpu_count()):
        mol = Mol.from_gjf(gjf)
        return cls(mol, jobname=gjf[:-4], hashline=mol.comment['hashline'], mem=mem, nproc=nproc)


    @staticmethod
    def writecom(jobname, nproc, mem, oldchk, hashline, charge, spin, sites):
        """

        :param jobname:
        :param nproc:
        :param oldchk:
        :param mem:
        :param charge:
        :param spin:
        :param sites:
        :return:
        """
        filename = jobname + '.com'
        content = []
        content.append('%rwf=' + jobname + '.rwf')
        content.append('%nosave')
        if oldchk != '':
            content.append('%oldchk=' + oldchk)
        content.append('%chk=' + jobname + '.chk')
        content.append('%nprocshared=' + str(nproc))
        content.append('%mem=' + str(mem))
        content.append(hashline)
        content.append('')
        content.append('Input Generated at ' + str(datetime.datetime.now()) + ' for ' + jobname)
        content.append('')
        content.append(str(charge) + ' ' + str(spin))
        if len(sites) != 0:
            for s in sites:
                x = str(round(s.x, 5))
                y = str(round(s.y, 5))
                z = str(round(s.z, 5))
                content.append(s.element.name + '\t\t' + x + '  ' + y + '  ' + z)
        else:
            content.append('')
        content.append('')
        with open(filename, 'w') as f:
            f.write('\r\n'.join(content))

    # @classmethod
    # def from_log(cls, log, mem='20GB', nproc=multiprocessing.cpu_count()):
    #     mol =
    #
    #
    #
    # @staticmethod
    # def parselog(logname, parsertype):
    #
    #     if parsertype == 'nicsxy':
    #         with open(logname, 'r') as f:
    #             lines = f.readlines()
    #
    #         ns = 0
    #         for i in range(len(lines)):
    #             if 'Charge' in lines[i] and 'Multiplicity' in lines[i]:
    #                 cursor = i + 1
    #                 while len(lines[cursor].strip().split()) == 4:
    #                     ns += 1
    #                     cursor += 1
    #                 break
    #
    #         tensor_zzs = []
    #         for i in range(len(lines)):
    #             if 'GIAO Magnetic shielding tensor' in lines[i]:
    #                 for j in range(i + 1, i + 1 + ns * 5):
    #                     if 'Bq' in lines[j]:
    #                         tensor_zzs.append(-float(lines[j + 3].strip().split()[-1]))
    #                 break
    #         with open(logname[:-4] + '.zz', 'w') as f:
    #             for i in range(len(tensor_zzs)):
    #                 f.write(str(i) + '\t' + str(tensor_zzs[i]) + '\r\n')
    #
    #
    # def write_inp(self, instructions):
    #     if self.jobtype == 'nicsxy-backbone-total':
    #         if any([k not in instructions.keys() for k in ['step_size', 'height', 'maxnbq', 'normaldirection']]):
    #             sys.exit('value not found in instructions')
    #
    #         pts, pt_idx, xnumbers, xticks = self.mol.nics_line_scan_path(instructions['step_size'],
    #                                                                      instructions['nrings'], instructions['height'],
    #                                                                      instructions['normaldirection'])
    #         bqsites = []
    #         for p in pts:
    #             bqsites.append(Site('bq', p))
    #         chunksize = instructions['maxnbq']
    #         if len(bqsites) > instructions['maxnbq']:
    #             bqsites_list = [bqsites[i * chunksize:(i + 1) * chunksize] for i in
    #                             range((len(bqsites) + chunksize - 1) // chunksize)]
    #             for idx in range(len(bqsites_list)):
    #                 self.writecom(self.jobname + '-' + str(idx) + '-nicsxy-backbone-total', self.nproc, self.mem, '',
    #                               self.level, self.addkwds + ' NMR=GIAO SCF(MAXCYCLE=512)', self.mol.charge,
    #                               self.mol.multiplicity, self.mol.msites + bqsites_list[idx])
    #         else:
    #             self.writecom(self.jobname + '-0-nicsxy-backbone-total', self.nproc, self.mem, '', self.level,
    #                           self.addkwds + ' NMR=GIAO SCF=(qc,MAXCYCLE=512)', self.mol.charge,
    #                           self.mol.multiplicity, self.mol.msites + bqsites)
    #
    #         with open(self.jobname + '-nicsxy-plot.data', 'w') as f:
    #             f.write('\t'.join([str(round(n, 4)) for n in xticks]) + '\r\n')
    #             for i in range(len(pt_idx)):
    #                 f.write(str(round(xnumbers[i], 4)) + '\t' + str(pt_idx[i]) + '\r\n')
    #
    #     if self.jobtype == 'nicsxy-backbone-sigma':
    #         if any([k not in instructions.keys() for k in ['step_size', 'height', 'maxnbq', 'normaldirection']]):
    #             sys.exit('value not found in instructions')
    #
    #         pts, pt_idx, xnumbers, xticks = self.mol.nics_line_scan_path(instructions['step_size'],
    #                                                                      instructions['nrings'], instructions['height'],
    #                                                                      instructions['normaldirection'])
    #         bqsites = []
    #         for p in pts:
    #             bqsites.append(Site('bq', p))
    #         chunksize = instructions['maxnbq']
    #         sigma_structure = self.mol.nics_sigma_structure(normaldirection=1 - instructions['normaldirection'])
    #         if len(bqsites) > instructions['maxnbq']:
    #             bqsites_list = [bqsites[i * chunksize:(i + 1) * chunksize] for i in
    #                             range((len(bqsites) + chunksize - 1) // chunksize)]
    #             for idx in range(len(bqsites_list)):
    #                 self.writecom(self.jobname + '-' + str(idx) + '-nicsxy-backbone-sigma', self.nproc, self.mem, '',
    #                               self.level, self.addkwds + ' NMR=GIAO SCF(MAXCYCLE=512)', self.mol.charge,
    #                               self.mol.multiplicity, sigma_structure.msites + bqsites_list[idx])
    #         else:
    #             self.writecom(self.jobname + '-0-nicsxy-backbone-sigma', self.nproc, self.mem, '', self.level,
    #                           self.addkwds + ' NMR=GIAO SCF=(qc,MAXCYCLE=512)', self.mol.charge,
    #                           self.mol.multiplicity, sigma_structure.msites + bqsites)
