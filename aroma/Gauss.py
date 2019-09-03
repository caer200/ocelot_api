import multiprocessing
import datetime
import sys
import numpy as np
from Schema import Mol, Site, Element

class GaussProj:

    def __init__(self, mol, jobtype, jobname,
                 addkwds='scf=(xqc,fermi,noincfock,ndamp=35,conver=8,vshift=500,novaracc) int=(acc2e=16)',
                 level='wb97xd/def2svp', mem='20GB',
                 nproc=multiprocessing.cpu_count()):

        self.mol = mol
        self.jobtype = jobtype
        self.jobname = jobname
        self.addkwds = addkwds
        self.level = level
        self.mem = mem
        self.nproc = nproc

    @classmethod
    def from_gjf(cls, gjf):
        mol = Mol.from_gjf(gjf)
        keywords = mol.comment.split()[1:]  # remove #p
        addkwds = []
        jobtype = []
        for kw in keywords:
            if '/' in kw and 'iop' not in kw.lower():
                level = kw
            else:
                addkwds.append(kw)
            if 'opt' in kw.lower():
                jobtype += 'opt'
            if 'freq' in kw.lower():
                jobtype += 'freq'
            if 'nmr' in kw.lower():
                jobtype += 'nmr'
        return cls(mol, '-'.join(jobtype), gjf[:-4],
                 addkwds=' '.join(addkwds),
                 level=level, mem='20GB',
                 nproc=multiprocessing.cpu_count())

    @staticmethod
    def parselog(logname, parsertype):

        if parsertype == 'nicsxy':
            with open(logname, 'r') as f:
                lines = f.readlines()

            ns = 0
            for i in range(len(lines)):
                if 'Charge' in lines[i] and 'Multiplicity' in lines[i]:
                    cursor = i+1
                    while len(lines[cursor].strip().split()) == 4:
                        ns += 1
                        cursor += 1
                    break

            tensor_zzs = []
            for i in range(len(lines)):
                if 'GIAO Magnetic shielding tensor' in lines[i]:
                    for j in range(i+1, i+1+ns*5):
                        if 'Bq' in lines[j]:
                            tensor_zzs.append(-float(lines[j+3].strip().split()[-1]))
                    break
            with open(logname[:-4]+'.zz', 'w') as f:
                for i in range(len(tensor_zzs)):
                    f.write(str(i) + '\t' + str(tensor_zzs[i]) + '\r\n')

    @staticmethod
    def writecom(jobname, nproc, mem, oldchk, level, keywords, charge, spin, sites):
        """

        :param jobname:
        :param nproc:
        :param oldchk:
        :param mem:
        :param level:
        :param keywords:
        :param charge:
        :param spin:
        :param sites:
        :return:
        """
        filename = jobname + '.com'
        with open(filename, 'w') as f:
            f.write('%rwf=' + jobname + '.rwf\r\n')
            f.write('%nosave\r\n')
            if oldchk != '':
                f.write('%oldchk=' + oldchk + '\r\n')
            f.write('%chk=' + jobname + '.chk\r\n')
            f.write('%nprocshared=' + str(nproc) + '\r\n')
            f.write('%mem=' + str(mem) + '\r\n')
            f.write('#p ' + level + ' ' + keywords + ' \r\n')
            f.write('\r\n')
            f.write('Input Generated at ' + str(datetime.datetime.now()) + ' for ' + jobname + '\r\n')
            f.write('\r\n')
            f.write(str(charge) + ' ' + str(spin) + ' \r\n')
            if len(sites) != 0:
                for s in sites:
                    x = str(round(s.x, 5))
                    y = str(round(s.y, 5))
                    z = str(round(s.z, 5))
                    f.write(s.element.name + '\t\t' + x + '  ' + y + '  ' + z + '  \r\n')
            else:
                f.write('\r\n')
            f.write('\r\n')

    def write_inp(self, instructions):
        if self.jobtype == 'nicsxy-backbone-total':
            if any([k not in instructions.keys() for k in ['step_size', 'height', 'maxnbq', 'normaldirection']]):
                sys.exit('value not found in instructions')

            pts, pt_idx, xnumbers, xticks = self.mol.nics_line_scan_path(instructions['step_size'], instructions['nrings'], instructions['height'], instructions['normaldirection'])
            bqsites = []
            for p in pts:
                bqsites.append(Site('bq', p))
            chunksize = instructions['maxnbq']
            if len(bqsites) > instructions['maxnbq']:
                bqsites_list = [bqsites[i * chunksize:(i + 1) * chunksize] for i in range((len(bqsites) + chunksize - 1) // chunksize )]
                for idx in range(len(bqsites_list)):
                    self.writecom(self.jobname+'-'+str(idx)+'-nicsxy-backbone-total', self.nproc, self.mem, '', self.level, self.addkwds+' NMR=GIAO SCF(MAXCYCLE=512)', self.mol.charge,
                                  self.mol.multiplicity, self.mol.sites+bqsites_list[idx])
            else:
                self.writecom(self.jobname+'-0-nicsxy-backbone-total', self.nproc, self.mem, '', self.level, self.addkwds+' NMR=GIAO SCF=(qc,MAXCYCLE=512)', self.mol.charge,
                          self.mol.multiplicity, self.mol.sites+bqsites)

            with open(self.jobname + '-nicsxy-plot.data', 'w') as f:
                f.write('\t'.join([str(round(n, 4)) for n in xticks]) + '\r\n')
                for i in range(len(pt_idx)):
                    f.write(str(round(xnumbers[i], 4)) + '\t' + str(pt_idx[i]) + '\r\n')


        if self.jobtype == 'nicsxy-backbone-sigma':
            if any([k not in instructions.keys() for k in ['step_size', 'height', 'maxnbq', 'normaldirection']]):
                sys.exit('value not found in instructions')

            pts, pt_idx, xnumbers, xticks = self.mol.nics_line_scan_path(instructions['step_size'], instructions['nrings'], instructions['height'], instructions['normaldirection'])
            bqsites = []
            for p in pts:
                bqsites.append(Site('bq', p))
            chunksize = instructions['maxnbq']
            sigma_structure = self.mol.nics_sigma_structure(normaldirection=1-instructions['normaldirection'])
            if len(bqsites) > instructions['maxnbq']:
                bqsites_list = [bqsites[i * chunksize:(i + 1) * chunksize] for i in range((len(bqsites) + chunksize - 1) // chunksize )]
                for idx in range(len(bqsites_list)):
                    self.writecom(self.jobname+'-'+str(idx)+'-nicsxy-backbone-sigma', self.nproc, self.mem, '', self.level, self.addkwds+' NMR=GIAO SCF(MAXCYCLE=512)', self.mol.charge,
                                  self.mol.multiplicity, sigma_structure.sites+bqsites_list[idx])
            else:
                self.writecom(self.jobname+'-0-nicsxy-backbone-sigma', self.nproc, self.mem, '', self.level, self.addkwds+' NMR=GIAO SCF=(qc,MAXCYCLE=512)', self.mol.charge,
                              self.mol.multiplicity, sigma_structure.sites+bqsites)










